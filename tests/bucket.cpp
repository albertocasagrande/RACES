/**
 * @file bucket.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Bucket tests
 * @version 1.0
 * @date 2025-09-13
 *
 * @copyright Copyright (c) 2023-2025
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE bucket

#include <stdint.h>
#include <type_traits>
#include <functional>
#include <random>
#include <set>

#include <boost/test/unit_test.hpp>

#include "common.hpp"

#include "bucket.hpp"

#define DEFAULT_DATASET_SIZE 10000
#define DEFAULT_WRITE_CACHE_SIZE 700
#define DEFAULT_READ_CACHE_SIZE 900

class TestData
{
public:
    size_t data_1;
    int16_t data_2;

    CHECK_CONSTANT_SPACE_ON_DISK(data_1, data_2)

    TestData():
        data_1{0}, data_2{static_cast<int16_t>(0)}
    {}

    TestData(const size_t data_1, const int16_t data_2):
        data_1{data_1}, data_2{data_2}
    {}

    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<RACES::Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        archive & data_1
                & data_2;
    }

    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<RACES::Archive::Basic::In, ARCHIVE>, bool> = true>
    inline static TestData load(ARCHIVE& archive)
    {
        TestData data;

        archive & data.data_1
                & data.data_2;

        return data;
    }

    inline bool operator==(const TestData& data) const
    {
        return (data.data_1 == data_1) && (data.data_2 == data_2);
    }

    inline bool operator!=(const TestData& data) const
    {
        return !(data==*this);
    }
};

std::ostream& operator<<(std::ostream& os, const TestData& data)
{
    os << "(" << data.data_1 << "," << data.data_2 << ")";

    return os;
}

template<>
struct std::less<TestData>
{
    bool operator()(const TestData& a, const TestData& b) const
    {
        return (a.data_1<b.data_1) || ((a.data_1==b.data_1) && (a.data_2<b.data_2));
    }
};


template<typename TYPE>
TYPE create_data(const size_t& i);

template<>
size_t create_data(const size_t& i)
{
    return i;
}

template<>
TestData create_data(const size_t& i)
{
    return TestData(i, static_cast<uint16_t>(i));
}

template<typename TYPE>
struct BucketFixture
{
    RACES::Archive::Bucket<TYPE> bucket;
    std::set<TYPE> dataset;

    BucketFixture():
        bucket{get_a_temporary_path(), DEFAULT_WRITE_CACHE_SIZE}
    {
        for (size_t i=0; i<DEFAULT_DATASET_SIZE; ++i) {
            TYPE data = create_data<TYPE>(i);
            BOOST_CHECK_NO_THROW(bucket.push_back(data));
            dataset.emplace(std::move(data));
        }

        bucket.flush();
    }

    ~BucketFixture()
    {
        bucket.flush();

        std::filesystem::remove(bucket.path());
    }
};

struct BucketFixtureSizeT : public BucketFixture<size_t>
{
};

template<typename TYPE>
void create_bucket()
{
    using namespace RACES::Archive;
    namespace fs = std::filesystem;

    fs::path bucket_filepath = get_a_temporary_path();

    BOOST_CHECK(!fs::exists(bucket_filepath));

    Bucket<TYPE> bucket(bucket_filepath);

    BOOST_CHECK(fs::exists(bucket.path()));

    bucket.flush();

    fs::remove(bucket.path());

    BOOST_CHECK_THROW(Bucket<TYPE>("/"),
                      std::domain_error);

    BOOST_CHECK_THROW(Bucket<TYPE>(get_a_temporary_path(), 0),
                      std::domain_error);
}

BOOST_AUTO_TEST_CASE(create_bucket_size)
{
    create_bucket<size_t>();
}

BOOST_AUTO_TEST_CASE(create_bucket_data)
{
    create_bucket<TestData>();
}

template<typename TYPE>
void load_bucket(const RACES::Archive::Bucket<TYPE>& bucket)
{
    using namespace RACES::Archive;
    {
        Bucket<TYPE> load_bucket(bucket.path());
        BOOST_CHECK(load_bucket.size()==DEFAULT_DATASET_SIZE);
    }
}

template<typename TYPE>
void sequential_bucket(const RACES::Archive::Bucket<TYPE>& bucket)
{
    using namespace RACES::Archive;

    {
        Bucket<TYPE> reading_bucket(bucket.path(), DEFAULT_READ_CACHE_SIZE);

        BOOST_CHECK(reading_bucket.size()==DEFAULT_DATASET_SIZE);

        size_t i=0;
        for (const auto& value: reading_bucket) {
            BOOST_CHECK(create_data<TYPE>(i)==value);
            ++i;
        }
    }
}


template<typename TYPE>
void random_io_bucket(const RACES::Archive::Bucket<TYPE>& bucket)
{
    using namespace RACES::Archive;

    // create a list of indices and randomly shuffle them
    std::vector<size_t> indices(bucket.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937_64 gen(0);

    std::shuffle(indices.begin(), indices.end(), gen);

    // check whether the values in the dataset correspond to
    // those in the same positions in the bucket
    for (const auto& index: indices) {
        BOOST_CHECK(create_data<TYPE>(index)==bucket[index]);
    }
}

template<typename TYPE, typename RANDOM_GENERATOR>
TYPE
test_random_tour_on(const RACES::Archive::BucketRandomTour<TYPE, RANDOM_GENERATOR>& tour,
                                  const std::set<TYPE>& dataset)
{
    std::set<TYPE> local_dataset(dataset.begin(), dataset.end());

    TYPE last_value;
    for (TYPE value : tour) {
        auto found = local_dataset.find(value);
        BOOST_CHECK(found != local_dataset.end());

        last_value = value;
        if (found != local_dataset.end()) {
            local_dataset.erase(found);
        }
    }

    BOOST_CHECK(local_dataset.empty());

    return last_value;
}

template<typename TYPE>
void random_tour(const RACES::Archive::Bucket<TYPE>& bucket,
                 const std::set<TYPE>& dataset)
{
    using namespace RACES::Archive;

    std::mt19937_64 gen(0);

    // testing different tours with different generators
    std::set<TYPE> last_values;
    for (size_t i=0; i<5; ++i) {
        std::mt19937_64 gen(i);

        auto tour = bucket.random_tour(gen, DEFAULT_READ_CACHE_SIZE);

        last_values.insert(test_random_tour_on(tour, dataset));
    }
    BOOST_CHECK(last_values.size()>1);

    // testing different tours with the same generator
    last_values = std::set<TYPE>();
    for (size_t i=0; i<5; ++i) {
        auto tour = bucket.random_tour(gen, DEFAULT_READ_CACHE_SIZE);

        last_values.insert(test_random_tour_on(tour, dataset));
    }
    BOOST_CHECK(last_values.size()==1);

    // testing multiple time the same generator
    last_values = std::set<TYPE>();
    auto tour = bucket.random_tour(gen, DEFAULT_READ_CACHE_SIZE);
    for (size_t i=0; i<5; ++i) {
        last_values.insert(test_random_tour_on(tour, dataset));
    }

    BOOST_CHECK(last_values.size()==1);
}

template<typename TYPE>
void copy_bucket(const RACES::Archive::Bucket<TYPE>& bucket)
{
    using namespace RACES::Archive;

    BOOST_CHECK_NO_THROW(Bucket<TYPE> bucket2(bucket));

    try {
        Bucket<TYPE> bucket2(bucket);

        BOOST_CHECK(bucket.size() == bucket2.size());
        BOOST_CHECK(bucket.get_cache_size() == bucket2.get_cache_size());

        auto it = bucket2.begin();
        for (const TYPE& value: bucket) {
            BOOST_CHECK(value == *it);
            ++it;
        }
    } catch (std::exception& ex)
    {}
}

template<typename TYPE>
void shuffle_bucket(RACES::Archive::Bucket<TYPE>& bucket,
                    const std::set<TYPE>& dataset)
{
    using namespace RACES::Archive;
    namespace fs = std::filesystem;

    std::mt19937_64 gen(0);

    {
        BOOST_CHECK_NO_THROW(bucket.shuffle(gen, fs::temp_directory_path()));
    }

    {
        Bucket<TYPE> bucket2(bucket.path(), DEFAULT_READ_CACHE_SIZE);
        std::set<TYPE> local_dataset(dataset.begin(), dataset.end());

        BOOST_CHECK(bucket2.size()==dataset.size());

        for (const TYPE& data: bucket2) {
            auto found = local_dataset.find(data);
            BOOST_CHECK(found != dataset.end());

            if (found != local_dataset.end()) {
                local_dataset.erase(found);
            }
        }
        BOOST_CHECK(local_dataset.empty());
    }
}

BOOST_FIXTURE_TEST_SUITE( BucketSizeT, BucketFixture<size_t>)

BOOST_AUTO_TEST_CASE(load_bucket_size_t)
{
    load_bucket<size_t>(bucket);
}

BOOST_AUTO_TEST_CASE(sequential_bucket_size_t)
{
    sequential_bucket<size_t>(bucket);
}

BOOST_AUTO_TEST_CASE(random_io_bucket_size_t)
{
    random_io_bucket<size_t>(bucket);
}

BOOST_AUTO_TEST_CASE(random_tour_size_t)
{
    random_tour<size_t>(bucket, dataset);
}

BOOST_AUTO_TEST_CASE(copy_bucket_size_t)
{
    copy_bucket<size_t>(bucket);
}

BOOST_AUTO_TEST_CASE(shuffle_bucket_size_t)
{
    shuffle_bucket<size_t>(bucket, dataset);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE( BucketTestData, BucketFixture<TestData> )

BOOST_AUTO_TEST_CASE(load_bucket_data)
{
    load_bucket<TestData>(bucket);
}

BOOST_AUTO_TEST_CASE(sequential_bucket_data)
{
    sequential_bucket<TestData>(bucket);
}

BOOST_AUTO_TEST_CASE(random_io_bucket_data)
{
    random_io_bucket<TestData>(bucket);
}

BOOST_AUTO_TEST_CASE(random_tour_data)
{
    random_tour<TestData>(bucket, dataset);
}

BOOST_AUTO_TEST_CASE(copy_bucket_data)
{
    copy_bucket<TestData>(bucket);
}

BOOST_AUTO_TEST_CASE(shuffle_bucket_data)
{
    shuffle_bucket<TestData>(bucket, dataset);
}

BOOST_AUTO_TEST_SUITE_END()
