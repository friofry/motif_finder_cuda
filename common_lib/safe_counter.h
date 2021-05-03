#pragma once
#include <mutex>
#include <utility>
#include <vector>

class SafeCounter {
public:
    SafeCounter() = default;
    SafeCounter(unsigned int max_value)
        : max_value_(max_value)
    {}
    SafeCounter(const std::vector<uint32_t> &values)
        : values_(values)
    {}

    unsigned int get_and_increment()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        auto prev = value_;
        value_++;
        return prev;
    }

    uint32_t get_next_value(bool &ok, uint32_t &id)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        auto prev = value_;
        value_++;

        ok = false;
        if (!values_.empty() && prev < values_.size()) {
            ok = true;
            id = prev;
            return values_[prev];
        }
        return 0;
    }

    std::pair<unsigned int, unsigned int> get_and_increment_range(unsigned int count)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        auto ret = std::make_pair(value_, value_);
        value_ += count;
        if (value_ > max_value_) {
            value_ = max_value_;
        }
        ret.second = value_;
        return ret;
    }
    struct RangeInfo {
        uint32_t start;
        uint32_t end;
        uint32_t idx;
    };
    RangeInfo get_and_increment_range_info(unsigned int count)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        RangeInfo ret;
        ret.start = value_;
        ret.idx = value_/count;

        value_ += count;
        if (value_ > max_value_) {
            value_ = max_value_;
        }
        ret.end  = value_;
        return ret;
    }

private:
    mutable std::mutex mutex_;
    unsigned int value_ = 0;
    unsigned int max_value_ = 0;
    std::vector<uint32_t> values_;
};