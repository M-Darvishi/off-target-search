#include <iostream>
#include "GenomeLoader.h"

void simpleLoadTest();
void realDataTest();

int main() {
	simpleLoadTest();
	realDataTest();
}

void simpleLoadTest() {
	auto genome = loadFasta("data/TestRNA.fna");
	std::cout <<"Genome Length: "<< genome.size() << std::endl;
}

void realDataTest() {
	auto genome = loadFasta("../../../data/GCF_000005845.2_ASM584v2_genomic.fna");
	std::cout << "Genome Length: " << genome.size() << std::endl;
}
