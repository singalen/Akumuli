#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "resp.h"

using namespace Akumuli;

// Test inegers

BOOST_AUTO_TEST_CASE(Test_respstream_read_integer) {

    const char* buffer = ":1234567890\r\n";
    MemStreamReader stream(buffer, 14);
    RESPStream resp(&stream);
    BOOST_REQUIRE_EQUAL(resp.next_type(), RESPStream::INTEGER);
    uint64_t result = -1;
    BOOST_REQUIRE(resp.read_int(&result));
    BOOST_REQUIRE_EQUAL(result, 1234567890);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_integer_wrong_type) {

    const char* buffer = "+1234567890\r\n";
    MemStreamReader stream(buffer, 14);
    RESPStream resp(&stream);
    BOOST_REQUIRE_EQUAL(resp.next_type(), RESPStream::STRING);
    uint64_t result = -1;
    BOOST_REQUIRE(!resp.read_int(&result));
    BOOST_REQUIRE_EQUAL(result, -1);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_integer_bad_value) {

    const char* buffer = ":123fl\r\n";
    MemStreamReader stream(buffer, 14);
    RESPStream resp(&stream);
    uint64_t result = -1;
    BOOST_REQUIRE(!resp.read_int(&result));
    BOOST_REQUIRE_EQUAL(result, -1);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_integer_bad_end_seq) {

    const char* buffer = ":1234567890\r00";
    MemStreamReader stream(buffer, 14);
    RESPStream resp(&stream);
    uint64_t result = -1;
    BOOST_REQUIRE(!resp.read_int(&result));
    BOOST_REQUIRE_EQUAL(result, -1);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_integer_too_long) {

    const char* buffer = ":"
            "11111111111111111111"
            "22222222222222222222"
            "11111111111111111111"
            "22222222222222222222"
            "11110000000000000000"
            "\r\n";
    // Integer is too long
    MemStreamReader stream(buffer, 104);
    RESPStream resp(&stream);
    uint64_t result = -1;
    BOOST_REQUIRE(!resp.read_int(&result));
    BOOST_REQUIRE_EQUAL(result, -1);
}

// Test strings

BOOST_AUTO_TEST_CASE(Test_respstream_read_string) {

    const char* orig = "+foobar\r\n";
    MemStreamReader stream(orig, 10);
    RESPStream resp(&stream);
    BOOST_REQUIRE_EQUAL(resp.next_type(), RESPStream::STRING);
    const size_t buffer_size = RESPStream::STRING_LENGTH_MAX;
    Byte buffer[buffer_size];
    int bytes = resp.read_string(buffer, buffer_size);
    BOOST_REQUIRE(bytes > 0);
    BOOST_REQUIRE_EQUAL(bytes, 6);
    BOOST_REQUIRE_EQUAL(std::string(buffer, buffer + bytes), "foobar");
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_string_wrong_type) {

    const char* orig = ":foobar\r\n";
    MemStreamReader stream(orig, 10);
    RESPStream resp(&stream);
    const size_t buffer_size = RESPStream::STRING_LENGTH_MAX;
    Byte buffer[buffer_size];
    int bytes = resp.read_string(buffer, buffer_size);
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_string_small_buffer) {

    const char* orig = "+foobar\r\n";
    MemStreamReader stream(orig, 10);
    RESPStream resp(&stream);
    const size_t buffer_size = 4;
    Byte buffer[buffer_size];
    int bytes = resp.read_string(buffer, buffer_size);
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_string_large_string) {

    std::string orig = "+";
    for (int i = 0; i < RESPStream::STRING_LENGTH_MAX + 1; i++) {
        orig.push_back('X');
    }
    orig.push_back('\r');
    orig.push_back('\n');
    MemStreamReader stream(orig.data(), orig.size());
    RESPStream resp(&stream);
    const size_t buffer_size = RESPStream::STRING_LENGTH_MAX;
    Byte buffer[buffer_size];
    int bytes = resp.read_string(buffer, buffer_size);
    BOOST_REQUIRE(bytes < 0);
}

// Test bulk strings

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring) {

    const char* orig = "$6\r\nfoobar\r\n";
    MemStreamReader stream(orig, 13);
    RESPStream resp(&stream);
    BOOST_REQUIRE_EQUAL(resp.next_type(), RESPStream::BULK_STR);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes > 0);
    BOOST_REQUIRE_EQUAL(bytes, 6);
    BOOST_REQUIRE_EQUAL(std::string(buffer.begin(), buffer.begin() + bytes), "foobar");
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring_bad_type) {

    const char* orig = ":6\r\nfoobar\r\n";
    MemStreamReader stream(orig, 13);
    RESPStream resp(&stream);
    BOOST_REQUIRE_NE(resp.next_type(), RESPStream::BULK_STR);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring_bad_header_1) {

    const char* orig = "$f\r\nfoobar\r\n";
    MemStreamReader stream(orig, 13);
    RESPStream resp(&stream);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring_bad_header_2) {

    const char* orig = "$\r\nfoobar\r\n";
    MemStreamReader stream(orig, 13);
    RESPStream resp(&stream);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring_bad_header_3) {

    const char* orig = "$6r\nfoobar\r\n";
    MemStreamReader stream(orig, 13);
    RESPStream resp(&stream);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring_bad_len_1) {

    const char* orig = "$1\r\nfoobar\r\n";
    MemStreamReader stream(orig, 13);
    RESPStream resp(&stream);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring_bad_len_2) {

    const char* orig = "$7\r\nfoobar\r\n";
    MemStreamReader stream(orig, 13);
    RESPStream resp(&stream);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring_bad_tail) {

    const char* orig = "$6\r\nfoobar\n";
    MemStreamReader stream(orig, 12);
    RESPStream resp(&stream);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes < 0);
}

BOOST_AUTO_TEST_CASE(Test_respstream_read_bulkstring_too_large_to_handle) {

    std::string orig = "$10000000\r\n";
    for (int i = 10000000; i-->0;) {
        orig.push_back('x');
    }
    orig.push_back('\r');
    orig.push_back('\n');
    MemStreamReader stream(orig.data(), orig.size());
    RESPStream resp(&stream);
    std::vector<Byte> buffer;
    buffer.resize(RESPStream::BULK_LENGTH_MAX);
    int bytes = resp.read_bulkstr(buffer.data(), buffer.size());
    BOOST_REQUIRE(bytes < 0);
}
