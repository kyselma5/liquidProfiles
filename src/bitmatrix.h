#pragma once

#include <iostream>
#include <cstdint>
#include <stdexcept>
#include <set>
#include <sstream>
#include <numeric>

class BitMatrix {
    std::vector<std::vector<bool>> bits{};
    //uint64_t bits{}; // can speed things up, but would limit us to 8x8 matrixes
    size_t k{};

    static constexpr size_t NOT_DEFINED = -1;

    class CondensedGraph {
    public:
        struct Node {
            std::set<size_t> outEdges;
            std::set<size_t> reverseEdges;
            std::vector<size_t> condensedNodes;
            bool sink = false;
            size_t maxDistanceToSink = NOT_DEFINED;
            size_t finalDelegation = NOT_DEFINED;
        };
        std::vector<Node> nodes;
        size_t maxDistanceToCycle = 0;
        size_t maxCycleSize = 0;
        size_t sinkCount = 0;
        size_t sourceCount = 0;
        size_t isolatedCount = 0;
        size_t maxTreeSize = 0; // TODO

        bool identifySinks();
        size_t calculateMaxDistanceToCycle(size_t idx);
        bool reconstructDelegations();
        bool hasSingleDelegations();
        std::vector<size_t> reconstructDelegationEdges();
        void calculateStats();
    };

public:
    BitMatrix(size_t k_) : bits(k_, std::vector<bool>(k_, false)), k(k_) {}
    BitMatrix(size_t k_, const std::vector<std::vector<bool>> & bits_) : bits(bits_), k(k_) {}
    BitMatrix(const std::string &input);
    BitMatrix(const std::vector<size_t>& edges);
    BitMatrix(uint64_t mask, size_t k_);

    CondensedGraph m_condensed;
    
    bool at(size_t row, size_t col) const {
        if (row >= k || col >= k) throw std::out_of_range("Index out of range");
        return bits[row][col];
    }

    void set(size_t row, size_t col, bool value) {
        if (row >= k || col >= k) throw std::out_of_range("Index out of range");
        bits[row][col] = value;
    }

    bool operator<(const BitMatrix &other) const {
        if (k != other.k) return k < other.k;

        for (size_t i = 0; i < k; ++i) {
            if (bits[i] < other.bits[i]) return true;
            if (other.bits[i] < bits[i]) return false;
        }

        return false; 
    }

    bool operator == (const BitMatrix& other) const {
        return k == other.k && bits == other.bits;
    }

    size_t size() const {return k;}

    void print() const;

    bool isTransitive() const {
        return transitiveClosure() == *this;
    }

    BitMatrix transitiveClosure() const;

    // DFS on original graph for Kosaraju
    void dfsOriginal(size_t v, std::vector<bool>& visited, std::stack<size_t>& st) const {
        visited[v] = true;
        for (size_t u = 0; u < k; ++u) {
            if (at(v, u) && !visited[u]) {
                dfsOriginal(u, visited, st);
            }
        }
        st.push(v);
    }

    // DFS on reversed graph for Kosaraju
    void dfsReversed(size_t v, const std::vector<std::vector<size_t>>& rev,
                     std::vector<bool>& visited, std::vector<size_t>& component, size_t compId) const {
        visited[v] = true;
        component[v] = compId;
        for (size_t u : rev[v]) {
            if (!visited[u]) {
                dfsReversed(u, rev, visited, component, compId);
            }
        }
    }

    // Find Strongly Connected Components (SCCs)
    // Returns a vector "component" where component[i] = SCC id of vertex i
    std::vector<size_t> findSCCs() const ;
    void buildCondensedGraph () ;
    bool everyoneVoted() const ;
    // checks is this matrix can represent LP is it does, returns true and in out parameter is one of the possible delegation matrixes.
    bool isLiquidProfileFast(BitMatrix & out);
    bool isLiquidProfileSlow() const;
    bool dfsCandidateGraph(size_t v, std::vector<size_t>& edges) const;
};


