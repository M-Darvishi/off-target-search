#include <iostream>
#include "Validation.h"

bool validateSearchInput(const std::string& genome, const std::string& guide, int max_mismatch)
{
    if (genome.empty()) {
        std::cerr << "Error: genome is empty\n";
        return false;
    }

    if (guide.empty()) {
        std::cerr << "Error: guide is empty\n";
        return false;
    }

    if (max_mismatch < 0) {
        std::cerr << "Error: max_mismatch cannot be negative\n";
        return false;
    }

    return true;
}

bool validateIndexInput(const std::string& genome, int k)
{
    if (genome.empty()) {
        std::cerr << "Error: genome is empty\n";
        return false;
    }

    if (k <+ 1) {
        std::cerr << "Error: k must be positive\n";
        return false;
    }

    if (k > genome.length()) {
        std::cerr << "Error: k > genome length\n";
        return false;
    }

    return true;
}


bool validateKmerInput(const std::string& genome, const std::string& guide, int max_mismatch, int k)
{
    if (!validateSearchInput(genome, guide, max_mismatch))  return false;
    if (!validateIndexInput(genome, k))  return false;

    if (k > guide.length()) {
        std::cerr << "Error: k > guide length\n";
        return false;
    }

    return true;
}