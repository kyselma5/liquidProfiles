#include <iostream>
#include <cstdint>
#include <stdexcept>
#include <set>
#include <sstream>

class CondesedGraph {
public:
    struct Node {
        std::set<size_t> outEdges;
        std::set<size_t> reverseEdges;
        std::set<size_t> sinks;
        std::vector<size_t> condensedNodes;
        bool closed = false;
        bool cycle = false;
        size_t distanceToCycle = -1;
        size_t finalDelegation = -1;
    };
    std::vector<Node> nodes;

    bool identifySinks() {
        for (size_t i = 0; i < nodes.size(); i++) {
            if (nodes[i].outEdges.count(i) == 1){
                nodes[i].cycle = true;
                nodes[i].distanceToCycle = 0;
                nodes[i].finalDelegation = i;
                if (nodes[i].outEdges.size() != 1){
                    return false;
                }
            }
        }
        return true;
    }

    std::set<size_t> findSinks(size_t idx, std::vector<Node> & nodes) {
        if (!nodes[idx].sinks.empty()) {
            return nodes[idx].sinks;
        }

        if (nodes[idx].cycle) {
            return nodes[idx].sinks = {idx};
        }

        for (auto &i : nodes[idx].outEdges) {
            auto otherSinks = findSinks(i, nodes);
            nodes[idx].sinks.insert(otherSinks.begin(), otherSinks.end());
        }

        return nodes[idx].sinks;
    }

    bool checkEveryoneHasExactlyOneSink() {

        for(const auto & n:nodes) {
            if (n.sinks.size() != 1) {
                return false;
            }
        }
        return true;
    }

    size_t maxDistanceToCycle(size_t idx, std::vector<Node> & nodes) {
        if(nodes[idx].cycle) {
            return 0;
        }

        size_t acc = 0;
        for(auto & i:nodes[idx].outEdges){
            acc = std::max(acc,maxDistanceToCycle(i, nodes));
        }

        return nodes[idx].distanceToCycle = acc+1;
    }

    bool reconstructDelegations() {
        std::queue <size_t> q;
        for(size_t i = 0; i < nodes.size(); i++){
            if (nodes[i].distanceToCycle == 0){
                q.push(i);
                nodes[i].closed = true;
            }
        }

        while(!q.empty()){
            size_t idx = q.front();
            q.pop();

            for (auto & n:nodes[idx].reverseEdges) {
                if ((nodes[n].distanceToCycle > nodes[idx].distanceToCycle) && nodes[n].closed){
                    return false;
                }
                if (nodes[n].distanceToCycle == nodes[idx].distanceToCycle + 1){
                    nodes[n].closed = true;
                    nodes[n].finalDelegation = idx;
                    q.push(n);
                }
            }
        }
        
        bool res = true;
        for (const auto & n:nodes){
            res = res && n.closed;
        }
        return res;
    }

    bool hasSingleDelegations() {
        if (!identifySinks()) {
            return false;
        }
        
        for(size_t i = 0; i < nodes.size(); i++){
            findSinks(i, nodes);
        }

        if (!checkEveryoneHasExactlyOneSink()) {
            return false;
        }

        for(size_t i = 0; i < nodes.size(); i++){
            maxDistanceToCycle(i, nodes);
        }
        return reconstructDelegations();
    }

    std::vector<std::vector<bool>> reconstructDelegationMartix(){
        size_t count = 0;
        for (const auto & n:nodes) {
            count += n.condensedNodes.size();
        }
        std::vector<std::vector<bool>> bits;
        bits.assign(count, std::vector<bool>(count, false));

        for (const auto & n:nodes) {
            if(n.cycle) {
                if (n.condensedNodes.size() == 1) {
                    // cycle with one node meaning one self delegating voter
                    bits[n.condensedNodes[0]][n.condensedNodes[0]] = true;
                }
                else {
                    // more then one voter in this component -> we need to make cycle using those voter delegations
                    for(size_t i = 0; i < n.condensedNodes.size() - 1; i++) {
                        bits[n.condensedNodes[i]][n.condensedNodes[i+1]] = true;
                    }
                    bits[n.condensedNodes[n.condensedNodes.size() - 1]][ n.condensedNodes[0]] = true;
                }
            }
            else {
                // node is not self delegating cycle, meaning it delegates to some other component -> we delegate the vote to some voter in that component
                bits[n.condensedNodes[0]][nodes[n.finalDelegation].condensedNodes[0]] = true;
            }
        }
        return bits;
    }

    void print(std::ostream &os = std::cout) const {
        os << "CondensedGraph:\n";
        os << "Number of nodes: " << nodes.size() << "\n\n";

        for (size_t i = 0; i < nodes.size(); ++i) {
            const Node &n = nodes[i];

            os << "Node " << i << ":\n";

            // condensed nodes
            os << "  condensedNodes: ";
            for (auto x : n.condensedNodes)
                os << x << " ";
            os << "\n";

            // out edges
            os << "  outEdges: ";
            for (auto x : n.outEdges)
                os << x << " ";
            os << "\n";

            // reverse edges
            os << "  reverseEdges: ";
            for (auto x : n.reverseEdges)
                os << x << " ";
            os << "\n";

            // flags
            os << "  closed: " << (n.closed ? "true" : "false") << "\n";
            os << "  cycle: " << (n.cycle ? "true" : "false") << "\n";

            // distances
            os << "  distanceToCycle: ";
            if (n.distanceToCycle == (size_t)-1)
                os << "undefined";
            else
                os << n.distanceToCycle;
            os << "\n";

            os << "  finalDelegation: ";
            if (n.finalDelegation == (size_t)-1)
                os << "undefined";
            else
                os << n.finalDelegation;
            os << "\n";

            os << "\n";
        }
    }
};

class BitMatrix {
    std::vector<std::vector<bool>> bits{};
    //uint64_t bits{}; // can speed things up, but would limit us to 8x8 matrixes
    size_t k{};

public:
    BitMatrix(size_t k_) : bits(k_, std::vector<bool>(k_, false)), k(k_) {}

    BitMatrix(size_t k_, const std::vector<std::vector<bool>> & bits_) : bits(bits_), k(k_) {}

    BitMatrix(const std::string &input) {
        std::istringstream iss(input);
        std::string line;
        std::vector<std::vector<bool>> temp;

        while (std::getline(iss, line)) {
            if (line.empty()) continue;

            std::vector<bool> row;
            for (char c : line) {
                if (c == '0') row.push_back(false);
                else if (c == '1') row.push_back(true);
                else if (isspace(c)) continue; // ignoruj mezery
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
    
    BitMatrix(const std::vector<size_t>& edges)
    {
        k = edges.size();
        bits.assign(k, std::vector<bool>(k, false));

        for (size_t i = 0; i < k; ++i) {
            if (edges[i] >= k) {
                throw std::out_of_range("Edge index out of range");
            }

            bits[i][edges[i]] = true;
        }
    }

    BitMatrix(uint64_t mask, size_t k_) : bits(k_, std::vector<bool>(k_, false)), k(k_) {
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

    bool at(size_t row, size_t col) const {
        if (row >= k || col >= k) throw std::out_of_range("Index out of range");
        return bits[row][col];
    }

    void set(size_t row, size_t col, bool value) {
        if (row >= k || col >= k) throw std::out_of_range("Index out of range");
        bits[row][col] = value;
    }

    bool operator<(const BitMatrix &other) const {
        if (k != other.k) return k < other.k; // nejdřív porovnáme velikost

        for (size_t i = 0; i < k; ++i) {
            if (bits[i] < other.bits[i]) return true;
            if (other.bits[i] < bits[i]) return false;
        }

        return false; 
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
    std::vector<size_t> findSCCs() const {
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

    std::vector<std::set<size_t>> buildComponentGraph() const {

        std::vector<size_t> component = findSCCs();

        size_t componentcount = 1 + *std::max_element(component.begin(), component.end());
        std::vector<std::set<size_t>> compAdj(componentcount);

        for (size_t i = 0; i < component.size(); i++){
            for (size_t j = 0; j < component.size(); j++){
                if (component[i] != component[j] && at(i,j)) {
                    compAdj[component[i]].insert(component[j]);
                }
            }
        }
        return compAdj;
    }

    BitMatrix buildComponentCondensationMartix() const {

        std::vector<size_t> component = findSCCs();

        size_t componentcount = 1 + *std::max_element(component.begin(), component.end());
        BitMatrix matrix = BitMatrix(componentcount);

        for (size_t i = 0; i < component.size(); i++){
            for (size_t j = 0; j < component.size(); j++){
                if (component[i] != component[j] && at(i,j)) {
                    matrix.set(component[i], component[j], true);
                }
            }
        }
        return matrix;
    }

    CondesedGraph buildCondesedGraph () const {
        std::vector<size_t> component = findSCCs();

        size_t componentcount = 1 + *std::max_element(component.begin(), component.end());
        CondesedGraph cg;

        cg.nodes.resize(componentcount);
        
        for (size_t i = 0; i < component.size(); i++){
            cg.nodes[component[i]].condensedNodes.push_back(i);
            for (size_t j = 0; j < component.size(); j++){

                if (at(i,j)) {
                    cg.nodes[component[i]].outEdges.insert(component[j]);
                    cg.nodes[component[j]].reverseEdges.insert(component[i]);
                }
            }
            cg.nodes[component[i]].closed = false;
        }
        return cg;
    }

    bool everyoneVoted() const {
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

    bool hasWitnessingGraphOld() const {
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

    BitMatrix matrix = BitMatrix(edges);

    //std::cout << "----- Delegation matrix -----\n";
    //matrix.print();
    //std::cout << "-----  Approval matrix  -----\n";
    auto transitiveClosure = matrix.transitiveClosure();
    std::cout << "martix - \n";
    matrix.print();
    std::cout << "closure - \n";
    transitiveClosure.print();
    //std::cout << "\n";
    auto ssc = transitiveClosure.findSCCs();
    std::cout << "scc - \n";
    for(const auto & i :ssc){
        std::cout << i << " ";
    }
    std::cout << std::endl;

    auto cg = transitiveClosure.buildComponentGraph();
    std::cout << "cg - \n";
    for(size_t i = 0; i < cg.size(); i++){
        std::cout << i << " - ";
        for(const auto & e:cg[i]){
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    auto condensedGraph = transitiveClosure.buildCondesedGraph();
    condensedGraph.print();


    std::cout << "hasWitnessingGraph " << condensedGraph.hasSingleDelegations() << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
}

void checkWitnessingGraphMethod(size_t k){

    for (uint64_t mask = 0; mask < (1ULL << (k * k)); ++mask) {

        BitMatrix m(mask, k);

        if (!m.isTransitive()) {
            continue;
        }

        auto cg = m.buildCondesedGraph();

        bool newmethod = cg.hasSingleDelegations();
        bool oldmethod = m.hasWitnessingGraphOld();
        if ( newmethod != oldmethod){
            m.print();
            std::cout << "old - " << newmethod << std::endl;
            std::cout << "new - " << oldmethod << std::endl;
            std::cout << std::endl;
        }
    }
}

void checkWitnessingGraphMethod2(size_t k){

    for (uint64_t mask = 0; mask < (1ULL << (k * k)); ++mask) {

        BitMatrix m(mask, k);

        if (m.isTransitive()) {
            auto cg = m.buildCondesedGraph();

            bool newmethod = cg.hasSingleDelegations();
            if (newmethod) {
                m.print();
                std::cout << "--------------\n";
                auto bits = cg.reconstructDelegationMartix();
                BitMatrix(bits.size(), bits).transitiveClosure().print();
                std::cout << "--------------\n";
                BitMatrix(bits.size(), bits).print();
                std::cout << "\n\n";
            }
        }
    }
}

void checkWitnessingGraphMethod3(size_t k){

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

    for (uint64_t mask = 0; mask < (1ULL << (k * k)); ++mask) {

        BitMatrix m(mask, k);

        if (m.isTransitive()) {
            auto cg = m.buildCondesedGraph();

            if (cg.hasSingleDelegations() != (mSet.count(m) == 1)) {
                m.print();
                std::cout << "\n";
            }
        }
    }
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

    //printAllMatrixiesOfgivenSize(k);
    checkWitnessingGraphMethod3(k);

    /*
    std::cout << "----------------------\n";
    auto m = BitMatrix( "1 0 0 0\n"
                        "1 0 1 1\n"
                        "1 0 0 0\n"
                        "0 0 0 1");

    auto cg = m.buildCondesedGraph();
    cg.print();
    std::cout << "\n\nhas single delegations? " << cg.hasSingleDelegations() << "\n\n";
    cg.print();
    auto bits = cg.reconstructDelegationMartix();
    BitMatrix(bits.size(), bits).print();
    BitMatrix(bits.size(), bits).transitiveClosure().print();
    */

    return 0;
}
