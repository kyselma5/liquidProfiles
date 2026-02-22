#include <iostream>
#include <cstdint>
#include <stdexcept>
#include <set>
#include <sstream>
#include <numeric>

constexpr size_t NOT_DEFINED = -1;

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

    bool identifySinks() {
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

    size_t maxDistanceToCycle(size_t idx) {
        if (nodes[idx].maxDistanceToSink != NOT_DEFINED) {
            return nodes[idx].maxDistanceToSink;
        }
        if(nodes[idx].sink) {
            return 0;
        }

        size_t acc = 0;
        for(auto & i:nodes[idx].outEdges){
            acc = std::max(acc,maxDistanceToCycle(i));
        }
        return nodes[idx].maxDistanceToSink = acc+1;
    }

    bool reconstructDelegations() {
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

    bool hasSingleDelegations() {
        if (!identifySinks()) {
            return false;
        }

        for(size_t i = 0; i < nodes.size(); i++){
            maxDistanceToCycle(i);
        }
        return reconstructDelegations();
    }

    std::vector<size_t> reconstructDelegationEdges(){
        size_t count = 0;
        for (const auto & n:nodes) {
            count += n.condensedNodes.size();
        }

        std::vector<size_t> edges(count, NOT_DEFINED);

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
            os << "  cycle: " << (n.sink ? "true" : "false") << "\n";

            // distances
            os << "  distanceToCycle: ";
            if (n.maxDistanceToSink == NOT_DEFINED)
                os << "undefined";
            else
                os << n.maxDistanceToSink;
            os << "\n";

            os << "  finalDelegation: ";
            if (n.finalDelegation == NOT_DEFINED)
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

    void print() const {
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < k; ++j) {
                std::cout << at(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    size_t size() {return k;}

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

    CondensedGraph buildCondensedGraph () const {
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

    // checks is this matrix can represent LP is it does, returns true and in out parameter is one of the possible delegation matrixes.
    bool isLiquidProfileFast(BitMatrix & out) const {
        if (isTransitive() && everyoneVoted()) {
            auto cg = buildCondensedGraph();

            if (cg.hasSingleDelegations()) {
                out = BitMatrix(cg.reconstructDelegationEdges());
                return true;
            }
        }
        return false;
    }

    bool isLiquidProfileSlow() const {
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
    return 0;
}
