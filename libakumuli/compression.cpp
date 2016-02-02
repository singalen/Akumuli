#include "compression.h"
#include "util.h"

#include <unordered_map>
#include <algorithm>
#include <iostream> // TODO: remove me

namespace Akumuli {

StreamOutOfBounds::StreamOutOfBounds(const char* msg) : std::runtime_error(msg)
{
}

struct PrevValPredictor {
    uint64_t last_value;
    PrevValPredictor(int) : last_value(0)
    {
    }
    uint64_t predict_next() const {
        return last_value;
    }
    void update(uint64_t value) {
        last_value = value;
    }
};

struct FcmPredictor {
    std::vector<uint64_t> table;
    uint64_t last_hash;
    const uint64_t MASK_;

    FcmPredictor(size_t table_size)
        : last_hash(0ull)
        , MASK_(table_size - 1)
    {
        assert((table_size & MASK_) == 0);
        table.resize(table_size);
    }

    uint64_t predict_next() const {
        return table[last_hash];
    }

    void update(uint64_t value) {
        table[last_hash] = value;
        last_hash = ((last_hash << 6) ^ (value >> 48)) & MASK_;
    }
};

struct DfcmPredictor {
    std::vector<uint64_t> table;
    uint64_t last_hash;
    uint64_t last_value;
    const uint64_t MASK_;

    //! C-tor. `table_size` should be a power of two.
    DfcmPredictor(int table_size)
        : last_hash (0ul)
        , last_value(0ul)
        , MASK_(table_size - 1)
    {
       assert((table_size & MASK_) == 0);
       table.resize(table_size);
    }

    uint64_t predict_next() const {
        return table.at(last_hash) + last_value;
    }

    void update(uint64_t value) {
        table[last_hash] = value - last_value;
        last_hash = ((last_hash << 2) ^ ((value - last_value) >> 40)) & MASK_;
        last_value = value;
    }
};

typedef FcmPredictor PredictorT;

static const int PREDICTOR_N = 1 << 10;

static inline void encode_value(Base128StreamWriter& wstream, uint64_t diff, unsigned char flag) {
    int nbytes = (flag & 7) + 1;
    int nshift = (64 - nbytes*8)*(flag >> 3);
    diff >>= nshift;
    switch(nbytes) {
    case 8:
        wstream.put_raw(diff);
        break;
    case 7:
        wstream.put_raw(static_cast<unsigned char>(diff & 0xFF));
        diff >>= 8;
    case 6:
        wstream.put_raw(static_cast<unsigned char>(diff & 0xFF));
        diff >>= 8;
    case 5:
        wstream.put_raw(static_cast<unsigned char>(diff & 0xFF));
        diff >>= 8;
    case 4:
        wstream.put_raw(static_cast<uint32_t>(diff & 0xFFFFFFFF));
        diff >>= 32;
        break;
    case 3:
        wstream.put_raw(static_cast<unsigned char>(diff & 0xFF));
        diff >>= 8;
    case 2:
        wstream.put_raw(static_cast<unsigned char>(diff & 0xFF));
        diff >>= 8;
    case 1:
        wstream.put_raw(static_cast<unsigned char>(diff & 0xFF));
    }
}

static inline uint64_t decode_value(Base128StreamReader& rstream, unsigned char flag) {
    uint64_t diff = 0ul;
    int nbytes = (flag & 7) + 1;
    for (int i = 0; i < nbytes; i++) {
        uint64_t delta = rstream.read_raw<unsigned char>();
        diff |= delta << (i*8);
    }
    int shift_width = (64 - nbytes*8)*(flag >> 3);
    diff <<= shift_width;
    return diff;
}

size_t CompressionUtil::compress_doubles(std::vector<double> const& input,
                                         Base128StreamWriter&       wstream,
                                         size_t begin, size_t end)
{
    if (begin == 0 && end == 0) {
        end = input.size();
    }
    PredictorT predictor(PREDICTOR_N);
    uint64_t prev_diff = 0;
    unsigned char prev_flag = 0;
    const size_t num_iter = end - begin;
    for (size_t ix = 0; ix < num_iter; ix++) {
        union {
            double real;
            uint64_t bits;
        } curr = {};
        curr.real = input.at(ix + begin);
        uint64_t predicted = predictor.predict_next();
        predictor.update(curr.bits);
        uint64_t diff = curr.bits ^ predicted;

        int leading_zeros = 64;
        int trailing_zeros = 64;

        if (diff != 0) {
            trailing_zeros = __builtin_ctzl(diff);
        }
        if (diff != 0) {
            leading_zeros = __builtin_clzl(diff);
        }

        int nbytes;
        unsigned char flag;

        if (trailing_zeros > leading_zeros) {
            // this would be the case with low precision values
            nbytes = 8 - trailing_zeros / 8;
            if (nbytes > 0) {
                nbytes--;
            }
            // 4th bit indicates that only leading bytes are stored
            flag = 8 | (nbytes&7);
        } else {
            nbytes = 8 - leading_zeros / 8;
            if (nbytes > 0) {
                nbytes--;
            }
            // zeroed 4th bit indicates that only trailing bytes are stored
            flag = nbytes&7;
        }

        if (ix % 2 == 0) {
            prev_diff = diff;
            prev_flag = flag;
        } else {
            // we're storing values by pairs to save space
            unsigned char flags = (prev_flag << 4) | flag;
            wstream.put_raw(flags);
            encode_value(wstream, prev_diff, prev_flag);
            encode_value(wstream, diff, flag);
        }
    }
    if ((end - begin) % 2 != 0) {
        // `input` contains odd number of values so we should use
        // empty second value that will take one byte in output
        unsigned char flags = prev_flag << 4;
        wstream.put_raw(flags);
        encode_value(wstream, prev_diff, prev_flag);
        //encode_value(wstream, 0, 0);
    }
    return end - begin;
}

void CompressionUtil::decompress_doubles(Base128StreamReader&     rstream,
                                         size_t                   numvalues,
                                         std::vector<double>     *output)
{
    PredictorT predictor(PREDICTOR_N);
    auto it = std::back_inserter(*output);
    int flags = 0;
    for (auto i = 0u; i < numvalues; i++) {
        unsigned char flag = 0;
        if (i % 2 == 0) {
            flags = (int)rstream.read_raw<unsigned char>();
            flag = static_cast<unsigned char>(flags >> 4);
        } else {
            flag = static_cast<unsigned char>(flags & 0xF);
        }
        uint64_t diff = decode_value(rstream, flag);
        union {
            uint64_t bits;
            double real;
        } curr = {};
        uint64_t predicted = predictor.predict_next();
        curr.bits = predicted ^ diff;
        predictor.update(curr.bits);
        *it++ = curr.real;
    }
}


static std::vector<std::pair<size_t, size_t>> split_chunk(const UncompressedChunk& data, int N) {
    std::vector<std::pair<size_t, size_t>> indexes;
    size_t sz = data.paramids.size();
    size_t base = 0;
    for (int i = 0; i < N; i++) {
        // from base to sz/N + something
        auto begin = base;
        auto end = base + sz/N;
        if (begin == end) {
            break;
        }
        auto lastid = data.paramids.at(end - 1);
        auto it = std::find_if(data.paramids.cbegin() + end, data.paramids.cend(), [lastid](aku_ParamId id) { return id != lastid; });
        end = it - data.paramids.cbegin();
        indexes.push_back(std::make_pair(begin, end));
        base = end;
        if (base == sz) {
            break;
        }
    }
    return indexes;
}


using ZDeltaRleWriter = DeltaStreamWriter<ZigZagStreamWriter<RLEStreamWriter<int64_t>, int64_t>, int64_t>;

using DeltaRleWriter = DeltaStreamWriter<RLEStreamWriter<uint64_t>, uint64_t>;

using ZDeltaRleReader = DeltaStreamReader<ZigZagStreamReader<RLEStreamReader<int64_t>, int64_t>, int64_t>;

using DeltaRleReader = DeltaStreamReader<RLEStreamReader<uint64_t>, uint64_t>;

static const int BLOCK_HEADER_SIZE = 4*sizeof(uint32_t) + 2*sizeof(aku_ParamId);

static std::vector<uint8_t> encode_block(const UncompressedChunk& data,
                                         size_t ixbegin, size_t ixend,
                                         size_t bytes_per_el,
                                         aku_Timestamp* out_min,
                                         aku_Timestamp* out_max) {
    if (bytes_per_el > 128) {
        // something goes really wrong
        throw StreamOutOfBounds("compression alg. failure");
    }
    std::vector<uint8_t> result;
    size_t buffer_size = (ixend - ixbegin) * bytes_per_el + BLOCK_HEADER_SIZE;
    result.resize(buffer_size);

    // Reserve space for metadata
    uint8_t* begin = result.data();
    const uint8_t* end = result.data() + buffer_size;
    aku_ParamId* first_id = reinterpret_cast<aku_ParamId*>(begin);
    begin += sizeof(aku_ParamId);
    aku_ParamId* last_id = reinterpret_cast<aku_ParamId*>(begin);
    begin += sizeof(aku_ParamId);
    uint32_t* num_elements = reinterpret_cast<uint32_t*>(begin);
    begin += sizeof(uint32_t);
    uint32_t* paramid_stream_sz = reinterpret_cast<uint32_t*>(begin);
    begin += sizeof(uint32_t);
    uint32_t* tmstamp_stream_sz = reinterpret_cast<uint32_t*>(begin);
    begin += sizeof(uint32_t);
    uint32_t* payload_stream_sz = reinterpret_cast<uint32_t*>(begin);
    begin += sizeof(uint32_t);

    try {
        // Data should contain at least one element
        assert(data.paramids.size() > 0);
        *num_elements = ixend - ixbegin;
        *first_id = data.paramids.at(ixbegin);
        *last_id = data.paramids.at(ixend-1);

        Base128StreamWriter wstream(begin, end);

        // Write ParamId stream
        DeltaRleWriter idwriter(wstream);
        for (size_t ix = ixbegin; ix < ixend; ix++) {
            idwriter.put(data.paramids[ix]);
        }
        idwriter.commit();
        *paramid_stream_sz = (uint32_t)idwriter.size();

        // Write Timestamp stream
        // Timestamps can decrease, because of that we should use zig-zag encoding. To be able to do this
        // we should cast each timestamp to int64_t. This is safe because timestamps will fit into int64_t
        // (and BTW they will fit into double).
        // TODO: calculate max timestamp that can be stored
        ZDeltaRleWriter tswriter(wstream);
        aku_Timestamp mints = ~0, maxts = 0;
        for (size_t ix = ixbegin; ix < ixend; ix++) {
            auto ts = data.timestamps[ix];
            mints = std::min(ts, mints);
            maxts = std::max(ts, maxts);
            tswriter.put(ts);
        }
        tswriter.commit();
        *tmstamp_stream_sz = (uint32_t)tswriter.size();
        *out_min = mints;
        *out_max = maxts;

        // Write Payload stream
        size_t payload_begin = wstream.size();
        CompressionUtil::compress_doubles(data.values, wstream, ixbegin, ixend);
        size_t payload_end = wstream.size();
        *payload_stream_sz = payload_end - payload_begin;
        result.resize(wstream.pos() - result.data());
        assert(static_cast<size_t>(*payload_stream_sz + *paramid_stream_sz + *tmstamp_stream_sz + BLOCK_HEADER_SIZE) == result.size());
    } catch(StreamOutOfBounds const&) {
        // TODO: add this information to self-monitoring
        std::vector<uint8_t> tmp = std::move(encode_block(data, ixbegin, ixend, bytes_per_el * 2, out_min, out_max));
        std::swap(result, tmp);
    }
    return std::move(result);
}

aku_Status CompressionUtil::encode_chunk( uint32_t           *n_elements
                                        , aku_Timestamp      *ts_begin
                                        , aku_Timestamp      *ts_end
                                        , ChunkWriter        *writer
                                        , const UncompressedChunk&  data)
{
    aku_MemRange available_space = writer->allocate();
    uint8_t* const origin = (uint8_t*)available_space.address;
    uint8_t* begin = (uint8_t*)available_space.address;
    uint8_t* end = begin + (available_space.length - 2*sizeof(uint32_t));  // 2*sizeof(aku_EntryOffset)

    if (available_space.length < 8) {
        return AKU_EOVERFLOW;
    }

    auto indexes = split_chunk(data, 4);
    *reinterpret_cast<uint32_t*>(begin) = (uint32_t)indexes.size();
    begin += sizeof(uint32_t);
    *reinterpret_cast<uint32_t*>(begin) = (uint32_t)data.values.size();
    begin += sizeof(uint32_t);

    aku_Timestamp mints = ~0, maxts = 0;
    for (auto ix: indexes) {
        aku_Timestamp itmin, itmax;

        std::vector<uint8_t> buffer = encode_block(data, ix.first, ix.second, 10, &itmin, &itmax);
        if (static_cast<size_t>(end - begin) > buffer.size()) {
            memcpy(begin, buffer.data(), buffer.size());
            begin += buffer.size();
            mints = std::min(mints, itmin);
            maxts = std::max(maxts, itmax);
        } else {
            return AKU_EOVERFLOW;
        }
    }
    *ts_begin = mints;
    *ts_end = maxts;
    *n_elements = (uint32_t)data.paramids.size();
    return writer->commit(begin - origin);
}

static aku_Status decode_block(UncompressedChunk *header, const uint8_t* begin, const uint8_t* end, size_t* out_size_bytes) {
    try {
        auto origin = begin;
        const auto BUFFER_SPACE = static_cast<uint32_t>(end - begin);

        if (end - begin < 32) {
            throw StreamOutOfBounds("chunk too small");
        }

        //aku_ParamId first_id = *reinterpret_cast<const aku_ParamId*>(begin);
        begin += sizeof(aku_ParamId);
        //aku_ParamId last_id = *reinterpret_cast<const aku_ParamId*>(begin);
        begin += sizeof(aku_ParamId);
        uint32_t num_elements = *reinterpret_cast<const uint32_t*>(begin);
        begin += sizeof(uint32_t);
        const uint32_t paramid_stream_sz = *reinterpret_cast<const uint32_t*>(begin);
        begin += sizeof(uint32_t);
        const uint32_t tmstamp_stream_sz = *reinterpret_cast<const uint32_t*>(begin);
        begin += sizeof(uint32_t);
        const uint32_t payload_stream_sz = *reinterpret_cast<const uint32_t*>(begin);
        begin += sizeof(uint32_t);

        auto bytes_to_read = payload_stream_sz + paramid_stream_sz + tmstamp_stream_sz + BLOCK_HEADER_SIZE;
        if (bytes_to_read > BUFFER_SPACE) {
            return AKU_EOVERFLOW;
        }

        Base128StreamReader rstream(begin, end);

        // Read ids
        DeltaRleReader idreader(rstream);
        for (uint32_t i = 0; i < num_elements; i++) {
            auto id = idreader.next();
            header->paramids.push_back(id);
        }

        // Read timestamps
        ZDeltaRleReader tsreader(rstream);
        for (uint32_t i = 0; i < num_elements; i++) {
            auto ts = tsreader.next();
            header->timestamps.push_back(ts);
        }

        CompressionUtil::decompress_doubles(rstream, num_elements, &header->values);

        *out_size_bytes = rstream.pos() - origin;
        assert(*out_size_bytes == bytes_to_read);

    } catch (StreamOutOfBounds const&) {
        return AKU_EBAD_DATA;
    }
    return AKU_SUCCESS;
}


// obsolete
aku_Status CompressionUtil::decode_chunk( UncompressedChunk   *header
                                        , const unsigned char *pbegin
                                        , const unsigned char *pend
                                        , uint32_t             nelements)
{
    auto num_blocks = *reinterpret_cast<const uint32_t*>(pbegin);
    pbegin += sizeof(uint32_t);
    auto total_elements = *reinterpret_cast<const uint32_t*>(pbegin);
    pbegin += sizeof(uint32_t);

    header->values.reserve(total_elements);
    header->paramids.reserve(total_elements);
    header->timestamps.reserve(total_elements);

    for (uint32_t i = 0; i < num_blocks; i++) {
        size_t block_size = 0;
        auto status = decode_block(header, pbegin, pend, &block_size);
        if (status != AKU_SUCCESS) {
            return status;
        }
        pbegin += block_size;
        if (pbegin > pend) {
            return AKU_EOVERFLOW;
        }
    }
    return AKU_SUCCESS;
}

template<class Fn>
bool reorder_chunk_header(UncompressedChunk const& header, UncompressedChunk* out, Fn const& f) {
    auto len = header.timestamps.size();
    if (len != header.values.size() || len != header.paramids.size()) {
        return false;
    }
    // prepare indexes
    std::vector<int> index;
    for (auto i = 0u; i < header.timestamps.size(); i++) {
        index.push_back(i);
    }
    std::stable_sort(index.begin(), index.end(), f);
    out->paramids.reserve(index.size());
    out->timestamps.reserve(index.size());
    out->values.reserve(index.size());
    for(auto ix: index) {
        out->paramids.push_back(header.paramids.at(ix));
        out->timestamps.push_back(header.timestamps.at(ix));
        out->values.push_back(header.values.at(ix));
    }
    return true;
}

bool CompressionUtil::convert_from_chunk_order(UncompressedChunk const& header, UncompressedChunk* out) {
    auto fn = [&header](int lhs, int rhs) {
        auto lhstup = header.timestamps[lhs];
        auto rhstup = header.timestamps[rhs];
        return lhstup < rhstup;
    };
    return reorder_chunk_header(header, out, fn);
}

bool CompressionUtil::convert_from_time_order(UncompressedChunk const& header, UncompressedChunk* out) {
    auto fn = [&header](int lhs, int rhs) {
        auto lhstup = header.paramids[lhs];
        auto rhstup = header.paramids[rhs];
        return lhstup < rhstup;
    };
    return reorder_chunk_header(header, out, fn);
}

}
