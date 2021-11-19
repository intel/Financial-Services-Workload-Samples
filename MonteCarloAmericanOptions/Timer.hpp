/*

Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iostream>

class ScopedTimer
{
public:
	ScopedTimer(double& res) : res_(res)
	{
		t1_ = std::chrono::system_clock::now();
	}

	~ScopedTimer()
	{
		t2_ = std::chrono::system_clock::now();
		res_ = std::chrono::duration_cast<std::chrono::milliseconds>(t2_ - t1_).count();
		std::cout << "res_" << res_ << std::endl;
	}
private:
	double &res_;
	std::chrono::time_point<std::chrono::system_clock> t1_;
	std::chrono::time_point<std::chrono::system_clock> t2_;
};

class Timer
{
public:
	Timer() { start(); }
	void start()
	{
		t1_ = std::chrono::system_clock::now();
	}
	void stop()
	{
		t2_ = std::chrono::system_clock::now();
	}
	double getTime() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(t2_ - t1_).count();
	}
private:
	std::chrono::time_point<std::chrono::system_clock> t1_;
	std::chrono::time_point<std::chrono::system_clock> t2_;
};


#endif
