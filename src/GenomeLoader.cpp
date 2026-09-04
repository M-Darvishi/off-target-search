#include <fstream>
#include <iostream>
#include <string>
#include "GenomeLoader.h"

std::string loadFasta(const std::string& file_path) {
	std::ifstream file(file_path);
	if (!file.is_open()) {
		std::cerr << "Cannot open file " << file_path << "!" << std::endl;
		return "";
	}

    std::string genome;
    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        if (line[0] == '>')
        {
            continue;
        }

        genome += line;
    }

    return genome;
}