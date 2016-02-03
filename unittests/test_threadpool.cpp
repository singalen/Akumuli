#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main
#include <boost/test/unit_test.hpp>
#include <apr.h>
#include <string>
#include <functional>

#include "thread_pool.h"

using namespace Akumuli;


BOOST_AUTO_TEST_CASE(test_threadpool_1)
{
    ThreadPool<std::function<std::string()>, std::string, 4> tpool;

    std::function<std::string()> fn1 = []() {
        return std::string("hello");
    };
    std::function<std::string()> fn2 = []() {
        return std::string("world");
    };

    tpool.push(fn1);
    tpool.push(fn2);

    std::string val1 = tpool.pop();
    std::string val2 = tpool.pop();
    BOOST_REQUIRE_EQUAL(val1, "hello");
    BOOST_REQUIRE_EQUAL(val2, "world");
}
