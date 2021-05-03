#pragma once
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cstdint>
#include "config.h"

class Timer
{
public:
    using Time = std::chrono::steady_clock;
    using float_sec = std::chrono::duration<float>;
    using float_time_point = std::chrono::time_point<Time, float_sec>;

public:
	Timer(const char *message = 0, unsigned int total=0)
		: start(now())
		, message(message)
		, total(total)
	{}

	~Timer()
	{
	    if (!_silence) {
	        std::cout << "TIMER: ";
	        print_elapsed();
        }
	}

	void print_gmps()
	{
		if (total)
		    std::cout << message << " GMPS: " << gmps() << std::endl;
	}

	double gmps() const
	{
        return (2.562890625 * total) / elapsed();
	}

	void silence()
	{
	    _silence = true;
	}

	void stop()
	{
	    if (!_silence)
	    print_elapsed();
	}

	void add_total(unsigned int count)
	{
	    total += count;
	}

    float_time_point now() const {
        return Time::now();
    }

	float elapsed() const
	{
        return (now() - start).count();
    }

    void print_elapsed()
    {
        if (message)
			std::cout << message << ": ";
		std::cout << elapsed() << "s. " << std::endl;
    }

	double estimate(unsigned current, unsigned total) const
	{
		unsigned estimate = 0u;
		if (current > 0)
		{
			auto time_elapsed = elapsed();
			estimate = (time_elapsed * double(total)) / double(current) - time_elapsed;
		}
		return estimate;
	}

	void estimate_print(unsigned current, unsigned all) const
	{
	    auto e = estimate(current, all);
        unsigned int now = elapsed();

        if (now - last_print > 0) {
            double part_done = (double)current / all;
            int percentage = 100.0 * part_done;
            std::cout << percentage << "% (" <<current << " of " << all << "): " << elapsed() << "s. /" << e << "s. left";
            if (total)
                std::cout << "; GMPS=" << gmps() * part_done;
            std::cout << std::endl;
            last_print = now;
        }
	}

	unsigned int get_total() const
	{
	    return total;
	}

	bool silent() const
	{
	    return _silence;
	}

private:
	float_time_point start;
	mutable int last_print { 0 };
	const char *message;
    unsigned int total;
	bool _silence { false };
};
