#include <iostream>
#include <sstream>
#include "bitmatrix.h"

bool BitMatrix::CondensedGraph::identifySinks() {
    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i].outEdges.count(i) == 1){
            nodes[i].sink = true;
            nodes[i].maxDistanceToSink = 0;
            nodes[i].finalDelegation = i;
            if (nodes[i].outEdges.size() != 1){
                return false;
            }
        }
    }
    return true;
}

size_t BitMatrix::CondensedGraph::calculateMaxDistanceToCycle(size_t idx) {
    if (nodes[idx].maxDistanceToSink != NOT_DEFINED) {
        return nodes[idx].maxDistanceToSink;
    }
    if(nodes[idx].sink) {
        return 0;
    }
    size_t acc = 0;
    for(auto & i:nodes[idx].outEdges){
        acc = std::max(acc,calculateMaxDistanceToCycle(i));
    }
    return nodes[idx].maxDistanceToSink = acc+1;
}

bool BitMatrix::CondensedGraph::reconstructDelegations() {
    std::queue<size_t> q;
    for(size_t i = 0; i < nodes.size(); i++){
        if(nodes[i].sink) {
            q.push(i);
        }
    }
    while(!q.empty()) {
        size_t curr = q.front();
        q.pop();
        nodes[curr].outEdges.insert(curr);
        for(auto & i:nodes[curr].reverseEdges){
            if(nodes[i].maxDistanceToSink == nodes[curr].maxDistanceToSink + 1) {
                if (nodes[i].outEdges == nodes[curr].outEdges) {
                    nodes[i].finalDelegation = curr;
                    q.push(i);
                }
                else {
                    return false;
                }
            }
        }
    }
    return true;
}

bool BitMatrix::CondensedGraph::hasSingleDelegations() {
    if (!identifySinks()) {
        return false;
    }
    for(size_t i = 0; i < nodes.size(); i++){
        calculateMaxDistanceToCycle(i);
    }
    return reconstructDelegations();
}

std::vector<size_t> BitMatrix::CondensedGraph::reconstructDelegationEdges(){
    size_t count = 0;
    for (const auto & n:nodes) {
        count += n.condensedNodes.size();
    }
    std::vector<size_t> edges(count, BitMatrix::NOT_DEFINED);
    for (const auto & n:nodes) {
        if(n.sink) {
            if (n.condensedNodes.size() == 1) {
                // sink with one node meaning one self delegating voter
                edges[n.condensedNodes[0]] = n.condensedNodes[0];
            }
            else {
                // more then one voter in this component -> we need to make cycle using those voter delegations
                for(size_t i = 0; i < n.condensedNodes.size() - 1; i++) {
                    edges[n.condensedNodes[i]] = n.condensedNodes[i+1];
                }
                edges[n.condensedNodes[n.condensedNodes.size() - 1]] = n.condensedNodes[0];
            }
        }
        else {
            // node is not self delegating cycle, meaning it delegates to some other component -> we delegate the vote to some voter in that component
            edges[n.condensedNodes[0]] = nodes[n.finalDelegation].condensedNodes[0];
        }
    }
    return edges;
}

void BitMatrix::CondensedGraph::calculateStats() {
    for (const auto & n:nodes) {
        if (n.reverseEdges.size() == 0){
            sourceCount++;
            if (n.sink) {
                isolatedCount++;
            }
        }
        if (n.sink) {
            maxTreeSize = std::max(maxTreeSize, n.reverseEdges.size());
            maxCycleSize = std::max(maxCycleSize, n.condensedNodes.size());
            sinkCount++;
        }
        maxDistanceToCycle = std::max(maxDistanceToCycle, n.maxDistanceToSink);
    }
}

BitMatrix::BitMatrix(const std::string &input) {
    std::istringstream iss(input);
    std::string line;
    std::vector<std::vector<bool>> temp;
    while (std::getline(iss, line)) {
        if (line.empty()) continue;
        std::vector<bool> row;
        for (char c : line) {
            if (c == '0') row.push_back(false);
            else if (c == '1') row.push_back(true);
            else if (isspace(c)) continue;
            else throw std::invalid_argument("Invalid character in input");
        }
        if (!temp.empty() && row.size() != temp[0].size()) {
            throw std::invalid_argument("Matrix is not rectangular");
        }
        temp.push_back(row);
    }
    if (temp.empty()) {
        throw std::invalid_argument("Empty matrix");
    }
    if (temp.size() != temp[0].size()) {
        throw std::invalid_argument("Matrix must be square");
    }
    k = temp.size();
    bits = std::move(temp);
}
    
BitMatrix::BitMatrix(const std::vector<size_t>& edges) {
    k = edges.size();
    bits.assign(k, std::vector<bool>(k, false));
    for (size_t i = 0; i < k; ++i) {
        if (edges[i] >= k) {
            throw std::out_of_range("Edge index out of range");
        }
        bits[i][edges[i]] = true;
    }
}

BitMatrix::BitMatrix(uint64_t mask, size_t k_) : bits(k_, std::vector<bool>(k_, false)), k(k_) {
    if (k > 8) {
        throw std::invalid_argument("Matrix too large for uint64_t");
    }
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            size_t bitIndex = i * k + j;
            bool value = (mask >> bitIndex) & 1ULL;
            bits[i][j] = value;
        }
    }
}

void BitMatrix::print() const {
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            std::cout << at(i, j) << " ";
        }
        std::cout << "\n";
    }
}

BitMatrix BitMatrix::transitiveClosure() const {
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

std::vector<size_t> BitMatrix::findSCCs() const {
    std::vector<size_t> component(k, -1);
    std::vector<bool> visited(k, false);
    std::stack<size_t> st;

    // First pass: DFS on original graph
    for (size_t i = 0; i < k; ++i) {
        if (!visited[i]) {
            dfsOriginal(i, visited, st);
        }
    }
    // Build reversed graph
    std::vector<std::vector<size_t>> rev(k);
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            if (bits[i][j]) {
                rev[j].push_back(i);
            }
        }
    }
    // Second pass: DFS on reversed graph
    visited.assign(k, false);
    size_t compId = 0;
    while (!st.empty()) {
        size_t v = st.top();
        st.pop();
        if (!visited[v]) {
            dfsReversed(v, rev, visited, component, compId);
            compId++;
        }
    }
    return component;
}

void BitMatrix::buildCondensedGraph () {
    std::vector<size_t> component = findSCCs();
    size_t componentcount = 1 + *std::max_element(component.begin(), component.end());
    CondensedGraph cg;
    cg.nodes.resize(componentcount);
    
    for (size_t i = 0; i < component.size(); i++){
        cg.nodes[component[i]].condensedNodes.push_back(i);
        for (size_t j = 0; j < component.size(); j++){
            // we also want to add self delegation edges within the same component
            if (at(i,j)) {
                cg.nodes[component[i]].outEdges.insert(component[j]);
                cg.nodes[component[j]].reverseEdges.insert(component[i]);
            }
        }
    }
    m_condensed = cg;
}

bool BitMatrix::everyoneVoted() const {
    for(size_t i = 0; i < k; i++){
        bool res = false;
        for(size_t j = 0; j < k; j++){
            res = at(i, j) || res;
        }
        if (!res){
            return false;
        }
    }
    return true;
}

bool BitMatrix::isLiquidProfileFast(BitMatrix & out) {
    if (isTransitive() && everyoneVoted()) {
        buildCondensedGraph();

        if (m_condensed.hasSingleDelegations()) {
            out = BitMatrix(m_condensed.reconstructDelegationEdges());
            m_condensed.calculateStats();
            return true;
        }
    }
    return false;
}

bool BitMatrix::isLiquidProfileSlow() const {
    if (k == 0) return true;

    std::vector<size_t> edges(k);
    return dfsCandidateGraph(0, edges);
}

bool BitMatrix::dfsCandidateGraph(size_t v, std::vector<size_t>& edges) const {
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