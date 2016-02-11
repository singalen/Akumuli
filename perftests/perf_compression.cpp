#include "compression.h"
#include "perftest_tools.h"
#include "datetime.h"

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <zlib.h>
#include <cstring>

#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <omp.h>

using namespace Akumuli;
namespace fs = boost::filesystem;

//! Generate time-series from random walk
struct RandomWalk {
    std::random_device                  randdev;
    std::mt19937                        generator;
    std::normal_distribution<double>    distribution;
    size_t                              N;
    std::vector<double>                 values;

    RandomWalk(double start, double mean, double stddev, size_t N)
        : generator(randdev())
        , distribution(mean, stddev)
        , N(N)
    {
        values.resize(N, start);
    }

    double generate(aku_ParamId id) {
        values.at(id) += distribution(generator);
        return (float)values.at(id);
    }

    void add_anomaly(aku_ParamId id, double value) {
        values.at(id) += value;
    }
};


UncompressedChunk read_data(fs::path path) {
    UncompressedChunk res;
    std::fstream in(path.c_str());
    std::string line;
    aku_ParamId base_pid = 1;
    std::map<std::string, aku_ParamId> pid_map;
    while(std::getline(in, line)) {
        std::istringstream lstr(line);
        std::string series, timestamp, value;
        std::getline(lstr, series, ',');
        std::getline(lstr, timestamp, ',');
        std::getline(lstr, value, ',');

        aku_ParamId id;
        auto pid_it = pid_map.find(series);
        if (pid_it == pid_map.end()) {
            pid_map[series] = base_pid;
            id = base_pid;
            base_pid++;
        } else {
            id = pid_it->second;
        }
        res.paramids.push_back(id);

        res.timestamps.push_back(DateTimeUtil::from_iso_string(timestamp.c_str()));

        res.values.push_back(std::stod(value));

    }
    return res;
}

struct Writer : ChunkWriter {
    ByteVector *out;
    Writer(ByteVector *out) : out(out) {}

    virtual aku_MemRange allocate() {
        aku_MemRange range = {
            out->data(),
            static_cast<uint32_t>(out->size())
        };
        return range;
    }

    //! Commit changes
    virtual aku_Status commit(size_t bytes_written) {
        out->resize(bytes_written);
        return AKU_SUCCESS;
    }
};

struct TestRunResults {
    // Akumuli stats
    std::string file_name;
    size_t uncompressed;
    size_t compressed;
    size_t nelements;

    double bytes_per_element;
    double compression_ratio;

    // Gzip stats
    double gz_bytes_per_element;
    double gz_compression_ratio;
    double gz_compressed;

    // Performance
    std::vector<double> perf;
    std::vector<double> gz_perf;
};

TestRunResults run_tests(fs::path path, bool bench) {
    TestRunResults runresults;
    runresults.file_name = fs::basename(path);

    auto header = read_data(path);

    const size_t UNCOMPRESSED_SIZE = header.paramids.size()*8    // Didn't count lengths and offsets
                                   + header.timestamps.size()*8  // because because this arrays contains
                                   + header.values.size()*8;     // no information and should be compressed
                                                                 // to a few bytes
    ByteVector out;
    out.resize(UNCOMPRESSED_SIZE * 4);

    Writer writer(&out);

    aku_Timestamp tsbegin, tsend;
    uint32_t n;
    auto status = CompressionUtil::encode_chunk(&n, &tsbegin, &tsend, &writer, header);
    if (status != AKU_SUCCESS) {
        std::cout << "Encoding error" << std::endl;
        exit(1);
    }

    // Compress using zlib

    // Ids copy (zlib need all input data to be aligned because it uses SSE2 internally)
    Bytef* pgz_ids = (Bytef*)aligned_alloc(64, header.paramids.size()*8);
    memcpy(pgz_ids, header.paramids.data(), header.paramids.size()*8);
    // Timestamps copy
    Bytef* pgz_ts = (Bytef*)aligned_alloc(64, header.timestamps.size()*8);
    memcpy(pgz_ts, header.timestamps.data(), header.timestamps.size()*8);
    // Values copy
    Bytef* pgz_val = (Bytef*)aligned_alloc(64, header.values.size()*8);
    memcpy(pgz_val, header.values.data(), header.values.size()*8);

    const auto gz_max_size = UNCOMPRESSED_SIZE*10;
    Bytef* pgzout = (Bytef*)aligned_alloc(64, gz_max_size);
    uLongf gzoutlen = gz_max_size;
    size_t total_gz_size = 0, id_gz_size = 0, ts_gz_size = 0, float_gz_size = 0;
    // compress param ids
    auto zstatus = compress(pgzout, &gzoutlen, pgz_ids, header.paramids.size()*8);
    if (zstatus != Z_OK) {
        std::cout << "GZip error" << std::endl;
        exit(zstatus);
    }
    total_gz_size += gzoutlen;
    id_gz_size = gzoutlen;
    gzoutlen = gz_max_size;
    // compress timestamps
    zstatus = compress(pgzout, &gzoutlen, pgz_ts, header.timestamps.size()*8);
    if (zstatus != Z_OK) {
        std::cout << "GZip error" << std::endl;
        exit(zstatus);
    }
    total_gz_size += gzoutlen;
    ts_gz_size = gzoutlen;
    gzoutlen = gz_max_size;
    // compress floats
    zstatus = compress(pgzout, &gzoutlen, pgz_val, header.values.size()*8);
    if (zstatus != Z_OK) {
        std::cout << "GZip error" << std::endl;
        exit(zstatus);
    }
    total_gz_size += gzoutlen;
    float_gz_size = gzoutlen;

    const float GZ_BPE = float(total_gz_size)/header.paramids.size();
    const float GZ_RATIO = float(UNCOMPRESSED_SIZE)/float(total_gz_size);


    const size_t COMPRESSED_SIZE = out.size();
    const float BYTES_PER_EL = float(COMPRESSED_SIZE)/header.paramids.size();
    const float COMPRESSION_RATIO = float(UNCOMPRESSED_SIZE)/COMPRESSED_SIZE;

    // Save compression stats

    runresults.uncompressed         = UNCOMPRESSED_SIZE;
    runresults.compressed           = COMPRESSED_SIZE;
    runresults.nelements            = header.timestamps.size();
    runresults.bytes_per_element    = BYTES_PER_EL;
    runresults.compression_ratio    = COMPRESSION_RATIO;

    runresults.gz_compressed        = id_gz_size + ts_gz_size + float_gz_size;
    runresults.gz_compression_ratio = GZ_RATIO;
    runresults.gz_bytes_per_element = GZ_BPE;

    // Try to decompress
    UncompressedChunk decomp;
    const unsigned char* pbegin = out.data();
    const unsigned char* pend = pbegin + out.size();
    CompressionUtil::decode_chunk(&decomp, pbegin, pend, header.timestamps.size());
    bool first_error = true;
    for (auto i = 0u; i < header.timestamps.size(); i++) {
        if (header.timestamps.at(i) != decomp.timestamps.at(i) && first_error) {
            std::cout << "Error, bad timestamp at " << i << std::endl;
            first_error = false;
        }
        if (header.paramids.at(i) != decomp.paramids.at(i) && first_error) {
            std::cout << "Error, bad paramid at " << i << std::endl;
            first_error = false;
        }
        double origvalue = header.values.at(i);
        double decvalue = decomp.values.at(i);
        if (origvalue != decvalue && first_error) {
            std::cout << "Error, bad value at " << i << std::endl;
            std::cout << "Expected: " << origvalue << std::endl;
            std::cout << "Actual:   " << decvalue << std::endl;
            first_error = false;
        }
    }

    if (bench) {
        // Bench compression process
        const int NRUNS = 100;
        aku_Status tstatus;
        ByteVector vec;
        runresults.perf.reserve(NRUNS);
        runresults.gz_perf.reserve(NRUNS);
        for (int i = 0; i < NRUNS; i++) {
            vec.resize(UNCOMPRESSED_SIZE*4);
            PerfTimer tm;
            Writer w(&vec);
            aku_Timestamp ts;
            uint32_t n;
            tstatus = CompressionUtil::encode_chunk(&n, &ts, &ts, &w, header);
            if (tstatus != AKU_SUCCESS) {
                std::cout << "Encoding error" << std::endl;
                exit(1);
            }
            runresults.perf.push_back(tm.elapsed());
        }

        for (int i = 0; i < NRUNS; i++) {
            PerfTimer tm;
            uLongf offset = 0;
            // compress param ids
            auto zstatus = compress(pgzout, &gzoutlen, pgz_ids, header.paramids.size()*8);
            if (zstatus != Z_OK) {
                std::cout << "GZip error" << std::endl;
                exit(zstatus);
            }
            offset += gzoutlen;
            gzoutlen = gz_max_size - offset;
            // compress timestamps
            zstatus = compress(pgzout + offset, &gzoutlen, pgz_ts, header.timestamps.size()*8);
            if (zstatus != Z_OK) {
                std::cout << "GZip error" << std::endl;
                exit(zstatus);
            }
            offset += gzoutlen;
            gzoutlen = gz_max_size - offset;
            // compress floats
            zstatus = compress(pgzout + offset, &gzoutlen, pgz_val, header.values.size()*8);
            if (zstatus != Z_OK) {
                std::cout << "GZip error" << std::endl;
                exit(zstatus);
            }
            runresults.gz_perf.push_back(tm.elapsed());
        }
    }
    return runresults;
}

int main(int argc, char** argv) {
    auto nthreads = omp_get_max_threads();
    omp_set_num_threads(nthreads);

    if (argc < 2) {
        std::cout << "Path to dataset required" << std::endl;
        exit(1);
    }

    bool benchmark = false;

    if (argc == 3 && std::string(argv[2]) == "bench") {
        benchmark = true;
    }

    // Iter directory
    fs::path dir{argv[1]};
    fs::directory_iterator begin(dir), end;
    std::vector<fs::path> files;
    std::list<TestRunResults> results;
    for (auto it = begin; it != end; it++) {
        files.push_back(*it);
    }
    std::sort(files.begin(), files.end());
    for (auto fname: files) {
        std::cout << "Run tests for " << fs::basename(fname) << std::endl;
        results.push_back(run_tests(fname, benchmark));
    }

    // Write table
    std::cout << "| File name | num elements | uncompressed | compressed | gzip compressed | ratio " <<
                 "| gzip ratio | bytes/el | gzip bytes/el | time (s) | gzip time (s) | speedup | compression factor (aku/gzip) " << std::endl;
    std::cout << "| ----- | ---- | ----- | ---- | ----- | ---- | ----- | ----- | ----- | ----- | ----- | ----- | -----| " << std::endl;
    for (auto const& run: results) {
        double maxdouble = std::numeric_limits<double>::max();
        auto perf = std::accumulate(run.perf.begin(), run.perf.end(), maxdouble, [](double a, double b) { return std::min(a, b); });
        auto gz_perf = std::accumulate(run.gz_perf.begin(), run.gz_perf.end(), maxdouble, [](double a, double b) { return std::min(a, b); });
        std::cout << //run.file_name << " | " <<
                     //run.nelements << " | " <<
                     //run.uncompressed << " | " <<
                     //run.compressed << " | " <<
                     //run.gz_compressed << " | " <<
                     //run.compression_ratio << " | " <<
                     //run.gz_compression_ratio << " | " <<
                     //run.bytes_per_element << " | " <<
                     //run.gz_bytes_per_element << " | " <<
                     perf << " | " <<
                     gz_perf << " | " <<
                     //(gz_perf/perf) << " | " <<
                     //(double(run.compressed)/run.gz_compressed) << " | " <<
                     std::endl;

    }
}
