#include <iostream>
#include <vector>
#include <random>
#include "bitmatrix.h"
#include "rules.h"
#include "axioms.h"

template<typename T>
void printVector(const std::vector<T> & v) { 
    for(const auto & e:v) {
        std::cout << e << ", ";
    }
    std::cout << std::endl;
}

void printDelegationAndApprovalMatrix(size_t k, const std::vector<size_t> & edges) {

    BitMatrix matrix = BitMatrix(k);

    for(size_t i = 0; i < k; i++){
        matrix.set(i, edges[i], true);
    }
    std::cout << "----- Delegation matrix -----\n";
    matrix.print();
    std::cout << "-----  Approval matrix  -----\n";
    (matrix.transitiveClosure()).print();
    std::cout << "\n";
}

// test for isLiquidProfile method. Works by generating all liquid profiles of given size. Then generates all possible matrixes of given size and prints all graphs that are LP but weren't recognized or vise versa.
void checkIsLiquidProfileMethod(size_t k){

    std::set<BitMatrix> mSet;

    uint64_t kk = 1; //number of possible delegations
    for(size_t i = 0; i < k; i++) {
        kk *= k;
    }

    std::vector<size_t> edges(k, 0);

    for(uint64_t num = 0; num < kk; num++) {
        uint64_t currentNum = num;

        //conversion from number i base k to k numbers representing edges
        for(size_t i = 0; i < k; i++) {
            edges[i] = currentNum%k;
            currentNum /= k;
        }

        mSet.insert(BitMatrix(edges).transitiveClosure());
    }

    uint64_t max = (1ULL << (k * k));

    for (uint64_t mask = 0; mask < max; ++mask) {

        if (mask%1000000 == 0) {
            std::cout << 100.0*mask/max << "percent done\n";
        }

        BitMatrix m(mask, k);
        BitMatrix out(k);

        if (m.isLiquidProfileFast(out) != (mSet.count(m) == 1)) {
            m.print();
            std::cout << "\n";
        }
    }
}

// simple generator 
void printAllMatrixesOfGivenSize(size_t k) {

    uint64_t kk = 1; //number of possible delegations
    for(size_t i = 0; i < k; i++) {
        kk *= k;
    }

    std::vector<size_t> edges(k, 0);

    for(uint64_t num = 0; num < kk; num++) {
        uint64_t currentNum = num;

        //conversion from number i base k to k numbers representing edges
        for(size_t i = 0; i < k; i++) {
            edges[i] = currentNum%k;
            currentNum /= k;
        }

        printDelegationAndApprovalMatrix(k, edges);
    }
}

int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <k> (size of matrixes to be generated/tested)\n";
        return 1;
    }

    size_t k = std::stoi(argv[1]);

    // generating 
    printAllMatrixesOfGivenSize(k);

    // testing of IsLiquidProfile Method (for development and debug)
    // checkIsLiquidProfileMethod(k);

    // here you can write your own matrix to check if it is LP and deconstruct possible delegations. 
    auto m = BitMatrix( "0 0 0 1 0 0 0 0 0\n"
                        "0 0 0 1 0 0 0 0 0\n"
                        "0 1 0 1 0 0 0 0 0\n"
                        "0 0 0 1 0 0 0 0 0\n"
                        "0 0 0 0 0 0 0 1 1\n"
                        "0 0 0 0 0 0 0 1 1\n"
                        "0 0 0 0 0 0 0 1 1\n"
                        "0 0 0 0 0 0 0 1 1\n"
                        "0 0 0 0 0 0 0 1 1");

    BitMatrix out(m.size());

    if (m.isLiquidProfileFast(out)) {
        std::cout << "----approval profile----\n";
        m.print();
        std::cout << "--possible delegations--\n";
        out.print();
        std::cout << "\n";
    }
    else {
        std::cout << "Matrix is not a Liquid Profile\n";
    }

    size_t numVoters = 20;
    std::vector<size_t> v(numVoters);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, numVoters - 1);
    for (size_t i = 0; i < numVoters; ++i) {
        v[i] = dist(gen);
    }

    std::cout << "------------------------\n";

    BitMatrix out2(numVoters);
    BitMatrix in(v);
     std::cout << "generated\n";
    BitMatrix in2 = in.transitiveClosure();
    std::cout << "transitiveClosure done\n";

    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "check LP "<<in2.isLiquidProfileFast(out2) << std::endl;

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time of calculation: " << elapsed.count() << " s\n";

    std::cout << "check LP done\n";
    //in.print();
    std::cout << "\n";
    //out2.print();
    std::cout << "check" << (in.transitiveClosure() == out2.transitiveClosure()) << std::endl;

    Rules r(in2);
    size_t commiteeSize = 2;
    printVector(r.approvalVoting(commiteeSize));
    printVector(r.sequentialPhragmen(commiteeSize));
    //printVector(r.thieleRule(commiteeSize, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
    printVector(r.thieleGreedy(commiteeSize, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
    printVector(r.thieleGreedy(commiteeSize, {1, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10}));
    printVector(r.MES(commiteeSize));
    printVector(r.GJCR(commiteeSize));

    AxiomChecker a(in2);
    std::cout << a.isJR(r.approvalVoting(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isJR(r.sequentialPhragmen(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isJR(r.MES(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isJR(r.GJCR(commiteeSize), commiteeSize) << std::endl;

    std::cout << a.isPJR(r.approvalVoting(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isPJR(r.sequentialPhragmen(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isPJR(r.MES(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isPJR(r.GJCR(commiteeSize), commiteeSize) << std::endl;

    std::cout << a.isEJR(r.approvalVoting(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isEJR(r.sequentialPhragmen(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isEJR(r.MES(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isEJR(r.GJCR(commiteeSize), commiteeSize) << std::endl;

    std::cout << a.isEJRplus(r.approvalVoting(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isEJRplus(r.sequentialPhragmen(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isEJRplus(r.MES(commiteeSize), commiteeSize) << std::endl;
    std::cout << a.isEJRplus(r.GJCR(commiteeSize), commiteeSize) << std::endl;
    return 0;
}
