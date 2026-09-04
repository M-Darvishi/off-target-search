#pragma once
#include <iostream>
#include <vector>
#include "Hit.h"

SlidingwindowHits slidingWindow_serial(std::string genome, std::string gRNA, int max_mismatch);
SlidingwindowHits slidingWindow_parallel(std::string genome, std::string gRNA, int max_mismatch , int threads_number);