/**
 * PRIVATE HEADER
 *
 * Data structures for main memory storage.
 *
 * Copyright (c) 2015 Eugene Lazin <4lazin@gmail.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#pragma once
#include "page.h"
#include "cursor.h"
#include "queryprocessor_framework.h"
#include "datetime.h"

#include <tuple>
#include <vector>
#include <algorithm>
#include <memory>
#include <mutex>

namespace Akumuli {

struct TimeSeriesValue {

    enum ValueType {
        BLOB,
        DOUBLE,
    };

    // Data members
    aku_Timestamp                           key_ts_;  // Key value (time)
    aku_ParamId                             key_id_;  // Key value (id)
    double                                  value;    // Numeric payload

    TimeSeriesValue();

    TimeSeriesValue(aku_Timestamp ts, aku_ParamId id, double value);

    aku_Timestamp get_timestamp() const;

    aku_ParamId get_paramid() const;

    aku_Sample to_result() const;

    void add_to_header(UncompressedChunk *chunk_header) const;

    friend bool operator < (TimeSeriesValue const& lhs, TimeSeriesValue const& rhs);

    //! Chunk order less then operator (id goes first, then goes timestamp)
    friend bool chunk_order_LT (TimeSeriesValue const& lhs, TimeSeriesValue const& rhs);

};


enum class SlidingWindowMatch {
    FUTURE_WRITE,
    PRESENT,
    PRESENT_SW_MOVED,
    LATE_WRITE,
};

struct SlidingWindow {
    virtual ~SlidingWindow() = default;
    virtual SlidingWindowMatch classify(aku_Timestamp ts) = 0;
};


/** Implements sequencer sliding window alg. (should be used for testing)
  * Doesn't uses system timer. Assumes there is no erroneus future writes.
  */
struct SimpleSlidingWindow {
    const aku_Duration           window_size_;
    aku_Timestamp                top_timestamp_;    //< Largest timestamp ever seen
    aku_Timestamp                checkpoint_;       //< Last checkpoint timestamp

    SimpleSlidingWindow(aku_Duration window_size)
        : window_size_(window_size)
        , top_timestamp_(0ull)
        , checkpoint_(0ull)
    {
    }

    virtual SlidingWindowMatch classify(aku_Timestamp ts) {
        if (ts < top_timestamp_) {
            auto delta = top_timestamp_ - ts;
            if (delta > window_size_) {
                return SlidingWindowMatch::LATE_WRITE;
            }
            return SlidingWindowMatch::PRESENT;
        }
        auto point = ts / window_size_;
        if (point > checkpoint_) {
            // Create new checkpoint
            checkpoint_ = point;
            top_timestamp_ = ts;
            return SlidingWindowMatch::PRESENT_SW_MOVED;
        }
        top_timestamp_ = ts;
        return SlidingWindowMatch::PRESENT;
    }
};

struct SimpleSlidingWindowV2 {
    const aku_Duration           window_size_;
    const aku_Duration           step_size_;
    aku_Timestamp                top_timestamp_;    //< Largest timestamp ever seen
    aku_Timestamp                checkpoint_;       //< Last checkpoint timestamp

    SimpleSlidingWindowV2(aku_Duration window_size)
        : window_size_(window_size)
        , step_size_(window_size/32ull)
        , top_timestamp_(0ull)
        , checkpoint_(0ull)
    {
        assert(step_size_);
    }

    virtual SlidingWindowMatch classify(aku_Timestamp ts) {
        if (ts < top_timestamp_) {
            auto delta = top_timestamp_ - ts;
            if (delta > window_size_) {
                return SlidingWindowMatch::LATE_WRITE;
            }
            return SlidingWindowMatch::PRESENT;
        }
        auto point = ts / step_size_;
        if (point > checkpoint_) {
            // Create new checkpoint
            checkpoint_ = point;
            top_timestamp_ = ts;
            return SlidingWindowMatch::PRESENT_SW_MOVED;
        }
        top_timestamp_ = ts;
        return SlidingWindowMatch::PRESENT;
    }
};


struct TimerBasedSlidingWindow {
    const aku_Duration window_size_;
    aku_Timestamp      top_timestamp_;    //< Largest timestamp ever seen
    aku_Timestamp      checkpoint_;       //< Last checkpoint timestamp

    TimerBasedSlidingWindow(aku_Duration window_size)
        : window_size_(window_size)
        , top_timestamp_(0ull)
        , checkpoint_(get_time() + window_size_/2)
    {
    }

    static aku_Timestamp get_time() {
        auto tp = std::chrono::high_resolution_clock::now();
        return DateTimeUtil::from_std_chrono(tp);
    }

    virtual SlidingWindowMatch classify(aku_Timestamp ts) {
        if (ts < (checkpoint_ - window_size_/2)) {
            // Late write
            return SlidingWindowMatch::LATE_WRITE;
        } else if (ts > checkpoint_) {
            // Check time
            auto newtp = get_time() + window_size_/2;
            if (ts > newtp) {
                return SlidingWindowMatch::FUTURE_WRITE;
            } else {
                // Success
                checkpoint_ = newtp;
                top_timestamp_ = ts;
                return SlidingWindowMatch::PRESENT_SW_MOVED;
            }
        }
        // Success
        top_timestamp_ = ts;
        return SlidingWindowMatch::PRESENT;
    }
};

/** Time-series sequencer.
  * @brief Akumuli can accept unordered time-series (this is the case when
  * clocks of the different time-series sources are slightly out of sync).
  * This component accepts all of them, filter out late writes and reorder
  * all the remaining samples by timestamp and parameter id.
  */
struct Sequencer {
    typedef std::vector<TimeSeriesValue> SortedRun;
    typedef std::shared_ptr<SortedRun>   PSortedRun;
    typedef std::mutex                   Mutex;
    typedef std::unique_lock<Mutex>      Lock;

    static const int RUN_LOCK_MAX_BACKOFF = 0x100;
    static const int RUN_LOCK_BUSY_COUNT = 0xFFF;
    static const int RUN_LOCK_FLAGS_MASK = 0x0FF;
    static const int RUN_LOCK_FLAGS_SIZE = 0x100;

    // TODO: space usage should be limited

    std::vector<PSortedRun>      runs_;             //< Active sorted runs
    std::vector<PSortedRun>      ready_;            //< Ready to merge
    PSortedRun                   key_;
    //
    TimerBasedSlidingWindow      sliding_window_;
    //
    mutable std::atomic_int      sequence_number_;  //< Flag indicates that merge operation is in progress and
                                                    //< search will return inaccurate results.
                                                    //< If progress_flag_ is odd - merge is in progress if it is
                                                    //< even - there is no merge and search will work correctly.
    mutable Mutex                runs_resize_lock_;
    mutable std::vector<RWLock>  run_locks_;
    const size_t                 c_threshold_;      //< Compression threshold

    Sequencer(aku_FineTuneParams const& config);

    /** Add new sample to sequence.
      * @brief Timestamp of the sample can be out of order.
      * @returns error code and flag that indicates whether of not new checkpoint is createf
      */
    std::tuple<aku_Status, int> add(TimeSeriesValue const& value);

    //! Simple merge and sync without compression. (depricated)
    void merge(Caller& caller, InternalCursor* cur);

    /** Merge all values (ts, id, offset, length)
      * and write it to target page.
      * caller and cur parameters used for communication with storage (error reporting).
      */
    aku_Status merge_and_compress(PageHeader* target, bool enforce_write=false);

    //! Close cache for writing, merge everything to page header.
    aku_Status close(PageHeader* target);

    /** Reset sequencer.
      * All runs are ready for merging.
      * @returns new sequence number.
      */
    int reset();

    /** Search in sequencer data.
      * @param query represents search query
      * @param sequence_number sequence number obtained with get_window function
      * @note search method follows common pattern used by all methods except sequence_number
      * parameter. This parameter is used to organize optimistic concurrency control. User must
      * call get_window fn and get current window and seq-number. This seq-number then passed to
      * search method. If seq-number is changed between calls to get_window and search - search
      * will be aborted and AKU_EBUSY.error code will be returned If merge occures during search -
      * search will be aborted and AKU_EBUSY error code will be returned.
      */
    void search(std::shared_ptr<QP::IQueryProcessor> query, int sequence_number) const;

    std::tuple<aku_Timestamp, int> get_window() const;

private:

    // move sorted runs to ready_ collection
    int make_checkpoint_(aku_Timestamp new_checkpoint);

    /** Check timestamp and make checkpoint if timestamp is large enough.
      * @returns error code and flag that indicates whether or not new checkpoint is created
      */
    std::tuple<aku_Status, int> check_timestamp_(aku_Timestamp ts);

    void filter(PSortedRun run, std::shared_ptr<QP::IQueryProcessor> query, std::vector<PSortedRun>* results) const;
};
}
