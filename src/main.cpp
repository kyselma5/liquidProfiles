#include <iostream>
#include <vector>
#include <random>
#include "bitmatrix.h"
#include "rules.h"
#include "axioms.h"



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
    //printAllMatrixesOfGivenSize(k);

    // testing of IsLiquidProfile Method (for development and debug)
    // checkIsLiquidProfileMethod(k);

    // here you can write your own matrix to check if it is LP and deconstruct possible delegations.
    /* 
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
    */
    size_t numVoters = k;
    std::vector<size_t> v(numVoters);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, numVoters - 1);
    for (size_t i = 0; i < numVoters; ++i) {
        v[i] = dist(gen);
    }
    
    BitMatrix out2(numVoters);
    BitMatrix in(v);
    BitMatrix in2 = in.transitiveClosure();
    std::cout << "transitiveClosure done\n";

    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "check LP "<<in2.isLiquidProfileFast(out2) << std::endl;

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time of calculation: " << elapsed.count() << " s\n";

    std::cout << "Sanity check " << (in.transitiveClosure() == out2.transitiveClosure()) << std::endl;

    std::vector<double> CC(k, 0);
    CC[0] = 1;
    std::vector<double> PAV(k, 0);
    for(size_t i = 0; i < k; i++){
        PAV[i] = 1/i;
    }

    Rules r(in2);
    size_t committeeSize = 5;

    auto resAV = r.approvalVoting(committeeSize);
    auto resSP = r.sequentialPhragmen(committeeSize);
    auto resCC = r.thieleGreedy(committeeSize, CC);
    auto resPAV = r.thieleGreedy(committeeSize, PAV);
    auto resMES = r.MES(committeeSize);
    auto resGJCR = r.GJCR(committeeSize);

    printVector(resAV);
    printVector(resSP);
    printVector(resCC);
    printVector(resPAV);
    printVector(resMES);
    printVector(resGJCR);

    AxiomChecker a(in2);

    std::cout << "     AV SP CC PAV MES GJCR\nJR   :";

    std::cout << a.isJRFast(resAV, committeeSize) << "  ";
    std::cout << a.isJRFast(resSP, committeeSize) << "  ";
    std::cout << a.isJRFast(resCC, committeeSize) << "  ";
    std::cout << a.isJRFast(resPAV, committeeSize) << "  ";
    std::cout << a.isJRFast(resMES, committeeSize) << "  ";
    std::cout << a.isJRFast(resGJCR, committeeSize) << "  ";

    std::cout << "\nPJR  :";

    std::cout << a.isPJRFast(resAV, committeeSize) << "  ";
    std::cout << a.isPJRFast(resSP, committeeSize) << "  ";
    std::cout << a.isPJRFast(resCC, committeeSize) << "  ";
    std::cout << a.isPJRFast(resPAV, committeeSize) << "  ";
    std::cout << a.isPJRFast(resMES, committeeSize) << "  ";
    std::cout << a.isPJRFast(resGJCR, committeeSize) << "  ";

    std::cout << "\nEJR  :";

    std::cout << a.isEJRFast(resAV, committeeSize) << "  ";
    std::cout << a.isEJRFast(resSP, committeeSize) << "  ";
    std::cout << a.isEJRFast(resCC, committeeSize) << "  ";
    std::cout << a.isEJRFast(resPAV, committeeSize) << "  ";
    std::cout << a.isEJRFast(resMES, committeeSize) << "  ";
    std::cout << a.isEJRFast(resGJCR, committeeSize) << "  ";

    std::cout << "\nEJR+ :";

    std::cout << a.isEJRplusFast(resAV, committeeSize) << "  ";
    std::cout << a.isEJRplusFast(resSP, committeeSize) << "  ";
    std::cout << a.isEJRplusFast(resCC, committeeSize) << "  ";
    std::cout << a.isEJRplusFast(resPAV, committeeSize) << "  ";
    std::cout << a.isEJRplusFast(resMES, committeeSize) << "  ";
    std::cout << a.isEJRplusFast(resGJCR, committeeSize) << "  ";

    std::cout << "\nCS   :";

    std::cout << a.isCS(resAV, committeeSize) << "  ";
    std::cout << a.isCS(resSP, committeeSize) << "  ";
    std::cout << a.isCS(resCC, committeeSize) << "  ";
    std::cout << a.isCS(resPAV, committeeSize) << "  ";
    std::cout << a.isCS(resMES, committeeSize) << "  ";
    std::cout << a.isCS(resGJCR, committeeSize) << "  ";

    std::cout << "\nSJR  :";

    std::cout << a.isSJRFast(resAV, committeeSize) << "  ";
    std::cout << a.isSJRFast(resSP, committeeSize) << "  ";
    std::cout << a.isSJRFast(resCC, committeeSize) << "  ";
    std::cout << a.isSJRFast(resPAV, committeeSize) << "  ";
    std::cout << a.isSJRFast(resMES, committeeSize) << "  ";
    std::cout << a.isSJRFast(resGJCR, committeeSize) << "  ";

    std::cout << "\nLR   :";

    std::cout << a.isLRFast(resAV, committeeSize) << "  ";
    std::cout << a.isLRFast(resSP, committeeSize) << "  ";
    std::cout << a.isLRFast(resCC, committeeSize) << "  ";
    std::cout << a.isLRFast(resPAV, committeeSize) << "  ";
    std::cout << a.isLRFast(resMES, committeeSize) << "  ";
    std::cout << a.isLRFast(resGJCR, committeeSize) << "  ";

    std::cout << "\nSEJR :";

    std::cout << a.isSEJRFast(resAV, committeeSize) << "  ";
    std::cout << a.isSEJRFast(resSP, committeeSize) << "  ";
    std::cout << a.isSEJRFast(resCC, committeeSize) << "  ";
    std::cout << a.isSEJRFast(resPAV, committeeSize) << "  ";
    std::cout << a.isSEJRFast(resMES, committeeSize) << "  ";
    std::cout << a.isSEJRFast(resGJCR, committeeSize) << "  ";

    std::cout << "\nSEJR+:";

    std::cout << a.isSEJRPlusFast(resAV, committeeSize) << "  ";
    std::cout << a.isSEJRPlusFast(resSP, committeeSize) << "  ";
    std::cout << a.isSEJRPlusFast(resCC, committeeSize) << "  ";
    std::cout << a.isSEJRPlusFast(resPAV, committeeSize) << "  ";
    std::cout << a.isSEJRPlusFast(resMES, committeeSize) << "  ";
    std::cout << a.isSEJRPlusFast(resGJCR, committeeSize) << "  ";

    std::cout << std::endl;
    return 0;
}
