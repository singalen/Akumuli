#pragma once

#include <atomic>

#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/lockfree/spsc_queue.hpp>

namespace Akumuli {

template<class Job, class Res, size_t S>
struct ThreadPool
{
    std::atomic<uint64_t> rdcounter_;
    std::atomic<uint64_t> wrcounter_;
    boost::thread_group tg_;
    boost::asio::io_service io_;
    boost::asio::io_service::work work_;
    typedef boost::lockfree::spsc_queue<Job, boost::lockfree::capacity<100>> JobQueue;
    typedef boost::lockfree::spsc_queue<Res, boost::lockfree::capacity<100>> ResQueue;
    std::array<JobQueue, S> jobs_;
    std::array<ResQueue, S> results_;

    ThreadPool()
        : rdcounter_{0ull}
        , wrcounter_{0ull}
        , work_(io_)
    {
        for (auto i = 0ull; i < S; i++) {
            auto &jobqueue = jobs_.at(i);
            auto &resqueue = results_.at(i);
            tg_.create_thread([&jobqueue, &resqueue]() {
                while(1) {
                    Job job;
                    if (jobqueue.pop(job)) {
                        Res result = job();
                        resqueue.push(result);
                    } else {
                        boost::this_thread::yield();
                    }
                }
            });
        }
    }

    void push(const Job& job) {
        auto val = rdcounter_++;
        auto &queue = jobs_.at(val % S);
        queue.push(job);
    }

    Res pop() {
        auto val = wrcounter_++;
        auto &queue = results_.at(val % S);
        while(1) {
            Res res;
            if (queue.pop(res)) {
                return res;
            } else {
                boost::this_thread::yield();
            }
        }
    }

};

}

