#include <iostream>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main
#include <boost/test/unit_test.hpp>

#include <apr.h>
#include "akumuli.h"
#include "storage_engine/blockstore.h"
#include "storage_engine/volume.h"
#include "log_iface.h"

void test_logger(aku_LogLevel tag, const char* msg) {
    BOOST_MESSAGE(msg);
}

struct AkumuliInitializer {
    AkumuliInitializer() {
        apr_initialize();
        Akumuli::Logger::set_logger(&test_logger);
    }
};

AkumuliInitializer initializer;

using namespace Akumuli;
using namespace Akumuli::StorageEngine;


static const std::vector<u32> CAPACITIES = { 8, 8 };  // two 64KB volumes
static const std::vector<std::string> VOLPATH = { "volume0", "volume1" };
static const std::string METAPATH = "metavolume";


static void create_blockstore() {
    Volume::create_new(VOLPATH[0].c_str(), CAPACITIES[0]);
    Volume::create_new(VOLPATH[1].c_str(), CAPACITIES[1]);
    MetaVolume::create_new(METAPATH.c_str(), 2, CAPACITIES.data());
}

static std::shared_ptr<FixedSizeFileStorage> open_blockstore() {
    auto bstore = FixedSizeFileStorage::open(METAPATH, VOLPATH);
    return bstore;
}


static void delete_blockstore() {
    apr_pool_t* pool;
    apr_pool_create(&pool, nullptr);
    apr_file_remove(METAPATH.c_str(), pool);
    apr_file_remove(VOLPATH[0].c_str(), pool);
    apr_file_remove(VOLPATH[1].c_str(), pool);
    apr_pool_destroy(pool);
}


BOOST_AUTO_TEST_CASE(Test_blockstore_0) {
    delete_blockstore();
    create_blockstore();
    auto bstore = open_blockstore();
    std::shared_ptr<Block> block;
    aku_Status status;
    // Should be unreadable
    std::tie(status, block) = bstore->read_block(0);
    BOOST_REQUIRE_NE(status, AKU_SUCCESS);

    // Append first block
    auto buffer = std::make_shared<Block>();
    buffer->get_data()[0] = 1;
    LogicAddr addr;
    std::tie(status, addr) = bstore->append_block(buffer);

    BOOST_REQUIRE_EQUAL(status, AKU_SUCCESS);

    BOOST_REQUIRE_EQUAL(addr, 0);

    // Should be readable now
    std::tie(status, block) = bstore->read_block(0);
    BOOST_REQUIRE_EQUAL(status, AKU_SUCCESS);

    const u8* block_data = block->get_data();
    size_t block_size = block->get_size();

    BOOST_REQUIRE_EQUAL(block_size, 4096);
    BOOST_REQUIRE_EQUAL(block_data[0], 1);

    delete_blockstore();
}

BOOST_AUTO_TEST_CASE(Test_blockstore_1) {
    delete_blockstore();
    create_blockstore();
    auto bstore = open_blockstore();


    // Fill data in
    auto buffer = std::make_shared<Block>();


    LogicAddr addr;
    aku_Status status;

    for (int i = 0; i < 17; i++) {
        buffer->get_data()[0] = static_cast<u8>(i);
        std::tie(status, addr) = bstore->append_block(buffer);
        BOOST_REQUIRE_EQUAL(status, AKU_SUCCESS);
    }
    BOOST_REQUIRE_EQUAL(addr, (2ull << 32));

    std::shared_ptr<Block> block;

    // Should be unreadable now
    std::tie(status, block) = bstore->read_block(0);
    BOOST_REQUIRE_EQUAL(status, AKU_EBAD_ARG);

    // Reada this block
    std::tie(status, block) = bstore->read_block(2ull << 32);
    BOOST_REQUIRE_EQUAL(status, AKU_SUCCESS);

    const u8* block_data = block->get_data();
    size_t block_size = block->get_size();

    BOOST_REQUIRE_EQUAL(block_size, 4096);
    BOOST_REQUIRE_EQUAL(block_data[0], 16);

    delete_blockstore();
}

