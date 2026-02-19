#include <iostream>
#include <cstdint>
#include <stdexcept>
#include <set>

class BitMatrix {
    std::vector<std::vector<bool>> bits{};
    //uint64_t bits{}; // can speed things up, but would limit us to 8x8 matrixes
    size_t k{};

public:
    BitMatrix(size_t k_) 
        : bits(k_, std::vector<bool>(k_, false)), k(k_) {}

    bool at(size_t row, size_t col) const {
        if (row >= k || col >= k) throw std::out_of_range("Index out of range");
        return bits[row][col];
    }

    void set(size_t row, size_t col, bool value) {
        if (row >= k || col >= k) throw std::out_of_range("Index out of range");
        bits[row][col] = value;
    }

    bool operator == (const BitMatrix& other) const {
        return k == other.k && bits == other.bits;
    }

    void print() const {
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < k; ++j) {
                std::cout << at(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    bool isTransitive() const {
        return transitiveClosure() == *this;
    }

    BitMatrix transitiveClosure() const {
        BitMatrix res(*this);

        for (size_t kIdx = 0; kIdx < k; ++kIdx) {
            for (size_t i = 0; i < k; ++i) {
                if (!res.at(i, kIdx)) continue;
                for (size_t j = 0; j < k; ++j) {
                    if (res.at(kIdx, j)){
                        res.set(i, j, true);
                    }
                }
            }
        }

        return res;
    }

    bool hasWitnessingGraph() const {
        if (k == 0) return true;

        std::vector<size_t> edges(k);
        return dfsCandidateGraph(0, edges);
    }

    bool dfsCandidateGraph(size_t v, std::vector<size_t>& edges) const {
        if (v == k) {
            BitMatrix candidate(k);
            for (size_t i = 0; i < k; ++i) {
                candidate.set(i, edges[i], true);
            }

            BitMatrix closure = candidate.transitiveClosure();
            return closure == *this;
        }

        for (size_t to = 0; to < k; ++to) {
            edges[v] = to;
            if (dfsCandidateGraph(v + 1, edges)) {
                return true;
            }
        }

        return false;
    }
};

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


void printAllMatrixiesOfgivenSize(size_t k) {

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
        std::cerr << "Usage: " << argv[0] << " <k>\n";
        return 1;
    }

    size_t k = std::stoi(argv[1]);

    printAllMatrixiesOfgivenSize(k);

    return 0;
}
