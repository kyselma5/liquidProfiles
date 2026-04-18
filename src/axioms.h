#pragma once

#include <iostream>
#include <vector>
#include <set>
#include "bitmatrix.h"
#include "other/maxflow2.h"

template<typename T>
void printContainer(const T & v) { 
    for(const auto & e:v) {
        std::cout << e << ", ";
    }
    std::cout << std::endl;
}

class AxiomChecker {
private:
    const BitMatrix& m;
    size_t V;
    size_t C;

public:
    AxiomChecker(const BitMatrix& matrix) : m(matrix), V(matrix.size()), C(matrix.size()) {}

    bool isJR(const std::set<size_t>& W, size_t k) {
        for(size_t c = 0; c < C; c++) {
            if(std::find(W.begin(), W.end(), c) == W.end()){
                continue;
            }

            // find group of voters supporting c
            // and count supporters, who are not represented in committee
            size_t countNotRepresented = 0;
            for(size_t v = 0; v < V; v++) {
                if(m.at(v, c)) {
                    bool represented = false;
                    for(size_t w:W){
                        if(m.at(v, w)) {
                            represented = true;
                            break;
                        }
                    }
                    if(!represented){
                        countNotRepresented++;
                    }
                }
            }

            // check if the number of not represented voters is enough to form 1 cohesive group
            if (countNotRepresented * k > V /*(*1)*/){
                return false;
            }
        }
        return true;
    }

    bool PJRhelp(size_t l, size_t k, 
        const std::set<size_t> & elected,
        const std::set<size_t> & supporters, 
        const std::set<size_t> & barrier, 
        const std::vector<std::set<size_t>> & V_v, 
        const std::vector<std::set<size_t>> & N_v, 
        const std::vector<size_t> & delegations) {
        
        if(l * V <= k * supporters.size()) {
            if (elected.size() < l) {
                return false;
            }
            else {
                return true;
            }
        }
        if(elected.size() >= l) {
            return true;
        }

        for(size_t n:barrier){

            auto supportersCopy = supporters;
            auto barrierCopy = barrier;
            auto electedCopy = elected;
            electedCopy.insert(delegations[n]);

            for(size_t i:barrier){
                if(delegations[n] == delegations[i]){
                    supportersCopy.insert(V_v[i].begin(), V_v[i].end());
                    supportersCopy.insert(i);
                    
                    barrierCopy.erase(i);
                    barrierCopy.insert(N_v[i].begin(), N_v[i].end());
                }
            }
            
            if (!PJRhelp(l, k, electedCopy, supportersCopy, barrierCopy, V_v, N_v, delegations)) {
                return false;
            }
        }
        return true;
    }

    bool isPJR(const std::set<size_t>& W, size_t k) {

        BitMatrix out(m.size());
        BitMatrix copy = m;
        copy.isLiquidProfileFast(out);
        std::vector<size_t> delegations(V, 0);
        for(size_t v = 0; v < V; v++) {
            for(size_t c = 0; c < C; c++) {
                if(out.at(v, c)) {
                    delegations[v] = c;
                }
            }
        }

        std::vector<std::set<size_t>> A_v_cap_W(V, std::set<size_t> {});
        std::vector<std::set<size_t>> V_v(V, std::set<size_t> {});
        std::vector<std::set<size_t>> N_v(V, std::set<size_t> {});

        for(size_t v = 0; v < V; v++) {
            for(size_t c = 0; c < C; c++) {
                if(m.at(v, c) && W.count(c)) {
                    A_v_cap_W[v].insert(c);
                }
            }
        }

        for(size_t v = 0; v < V; v++) {
            for(size_t r = 0; r < V; r++) {
                if((A_v_cap_W[r] == A_v_cap_W[v] && m.at(r, v))) {
                    V_v[v].insert(r);
                }
            }
        }

        for(size_t v = 0; v < V; v++) {
            for(size_t r = 0; r < V; r++) {
                size_t i = delegations[r];
                if(!V_v[v].count(r) && V_v[v].count(i)) {
                    N_v[v].insert(r);
                }
            }
        }

        std::vector<std::set<size_t>> elected(V, std::set<size_t> {});
        std::vector<size_t> candidateCount(V, 0);
        std::vector<size_t> rootedHere(V, 0);

        for(size_t v = 0; v < V; v++) {
            for(size_t c = 0; c < C; c++) {
                candidateCount[v] += m.at(v, c);
                rootedHere[c] += m.at(v, c);
            }
            for(auto w:W){
                if(m.at(v, w)) {
                    elected[v].insert(w);
                }
            }
        }

        for(size_t l = 1; l <= k; l++) {
            for(size_t v = 0; v < V; v++) {
                if(((m.at(v, v) && candidateCount[v] >= l) || (!m.at(v, v) && candidateCount[v] >= l - 1))
                && rootedHere[v] * k >= l * V) {

                    // voters in rootedHere are l-cohesive

                    std::set<size_t> alreadyElected = elected[v];
                    if(W.count(v)) {alreadyElected.insert(v);}

                    if(!PJRhelp(l, k, alreadyElected, V_v[v], N_v[v], V_v, N_v, delegations)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool isEJR(const std::set<size_t>& W, size_t k) {

        std::vector<size_t> electedCount(V, 0);
        std::vector<size_t> candidateCount(V, 0);

        for(size_t v = 0; v < V; v++) {
            for(size_t c = 0; c < C; c++) {
                candidateCount[v] += m.at(v, c);
            }
            for(auto w:W){
                electedCount[v] += m.at(v, w);
            }
        }

        for(size_t l = 1; l <= k; l++) {
            std::vector<size_t> unsatisfiedRootedHere(V, 0);
            for(size_t v = 0; v < V; v++) {
                if(electedCount[v] < l) {
                    for(size_t c = 0; c < C; c++) {
                        unsatisfiedRootedHere[c] += m.at(v, c);
                    }
                }
            }

            for(size_t v = 0; v < V; v++) {
                if((m.at(v, v) && candidateCount[v] >= l) || (!m.at(v, v) && candidateCount[v] >= l - 1)) {
                    if (unsatisfiedRootedHere[v] * k >= l * V) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool isEJRplus(const std::set<size_t>& W, size_t k) {
        // count of candidates in committee approved by each voter
        std::vector<size_t> electedCount(V, 0);
        for(size_t v = 0; v < V; v++) {
            for(size_t w:W){
                if (m.at(v, w)) {
                    electedCount[v]++;
                }
            }
        }

        for(size_t c = 0; c < C; c++){
            // filtering elected candidates
            bool elected = false;
            for(size_t w:W){
                if (w == c){
                    elected = true;
                    break;
                }
            }
            if (elected){
                continue;
            }

            for(size_t l = 1; l <= k; l++) {
                size_t underrepresentedSupportersCount = 0;
                for(size_t v = 0; v < V; v++) {
                    // counting all the voters approving c, which are also underrepresented.
                    if(m.at(v, c) && electedCount[v] < l){
                        underrepresentedSupportersCount++;
                    }
                }
                if(underrepresentedSupportersCount*k >= l*V){
                    return false;
                }
            }
        }

        return true;
    }

    bool isPR(const std::set<size_t>& W, size_t k) {
        if(V%k != 0) {
            return false;
        }
        max_flow<size_t> network(1+V+V+1);
        for(size_t v = 0; v < V; v++){
            network.add(0, v+1, 1);
            for(size_t w:W){
                if(m.at(v, w)){
                    network.add(1+v, 1+V+w, 1);
                }
            }
        }
        for(size_t w:W){
            network.add(1+V+w, 1+V+V, V/k);
        }

        return (V == network.calc(0, 1+V+V));
    }

    bool CShelper(size_t k, size_t t, size_t maxCandidate, const std::vector<size_t> & approvedCount, const std::vector<size_t> & electedCount) {
        // if number of candidates in T is larger than k, we can stop
        if (k < t){
            return true;
        }

        // try to add all possible candidates to T
        for(size_t c = maxCandidate; c < C; c++){

            // add candidate c to counter of approved candidates in alternative committee
            std::vector<size_t> newApprovedCount = approvedCount;
            for(size_t v = 0; v < V; v++) {
                if(m.at(v, c)) {
                    newApprovedCount[v]++;
                }
            }

            // count how many voters would be happier with T than W
            // also count potential maximum number of voters if they would approve all the candidates added to T later
            size_t countMoreHappyWithT = 0;
            size_t countPotentialMoreHappyWithT = 0;
            for(size_t v = 0; v < V; v++) {
                if(newApprovedCount[v] > electedCount[v]){
                    countMoreHappyWithT++;
                }
                if(newApprovedCount[v]+k-(t+1) > electedCount[v]){
                    countPotentialMoreHappyWithT++;
                }
            }

            // check the rule
            if (countMoreHappyWithT*k >= (t+1)*V){
                return false;
            }
            // if not enough people would be happy with best potential T, we can stop this branch anyway
            if (countPotentialMoreHappyWithT*k < (t+1)*V){
                continue;
            }
            if(!CShelper(k, t+1, c+1, newApprovedCount, electedCount)){
                return false;
            }
        }
        return true;
    }

    bool isCS(const std::set<size_t>& W, size_t k) {
        std::vector<size_t> electedCount(V, 0);
        std::vector<size_t> approvedCount(V, 0);

        // count approved count for each voter in W and set this count to 0 in T
        for(size_t v = 0; v < V; v++) {
            for(size_t w:W){
                if (m.at(v, w)) {
                    electedCount[v]++;
                }
            }
        }
        return CShelper(k, 0, 0, approvedCount, electedCount);
    }

    bool issJR(const std::set<size_t>& W, size_t k) {

        for(size_t c = 0; c < C; c++) {

            // count group of voters supporting c
            // check if everyone from this group is represented
            size_t countGroup = 0;
            bool everyoneRepresented = true;
            for(size_t v = 0; v < V; v++) {
                if(m.at(v, c)) {
                    countGroup++;

                    bool voterRepresented = false;
                    for(size_t w:W) {
                        if (m.at(v, w)) {
                            voterRepresented = true;
                            break;
                        }
                    }
                    everyoneRepresented &= voterRepresented;
                }
            }
            if(countGroup*k > V && !everyoneRepresented) {
                return false;
            }

        }
        return true;
    }

    bool isLR(const std::set<size_t>& W, size_t k) {

        std::vector<size_t> supportersCount(C, 0); // \sum_{v \in V} \mathbf{1}_{c \in A_v}
        std::vector<size_t> CAndA_CCount(C, 0); // \lvert c \cup A_c \rvert
        std::vector<size_t> electedCAndA_CCount(C, 0); // \lvert (c \cup A_c) \cap W \rvert

        for(size_t v = 0; v < V; v++) {
            for(size_t c = 0; c < C; c++) {
                supportersCount[c] += m.at(v, c);
                CAndA_CCount[c] += ((v == c) || m.at(c, v));
            }
        }

        for(size_t w:W) {
            for(size_t c = 0; c < C; c++) {
                electedCAndA_CCount[c] += ((w == c) || m.at(c, w));
            }
        }

        for(size_t l = 1; l <= k; l++) {
            for(size_t c = 0; c < C; c++) {

                if(supportersCount[c] * k < l * V) {
                    continue;
                }

                if(CAndA_CCount[c] >= l) {
                    if(electedCAndA_CCount[c] < l){
                        return false;
                    }
                }

                else {
                    if(electedCAndA_CCount[c] < CAndA_CCount[c]) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool issEJR(const std::set<size_t>& W, size_t k) {
        
        std::vector<size_t> electedCount(V, 0);
        std::vector<size_t> candidateCount(V, 0);
        std::vector<std::vector<size_t>> rootedHere(V);

        for(size_t v = 0; v < V; v++) {
            for(size_t c = 0; c < C; c++) {
                candidateCount[v] += m.at(v, c);
                if (m.at(v, c)) {
                    rootedHere[c].push_back(v);
                }
            }
            for(auto w:W){
                electedCount[v] += m.at(v, w);
            }
        }

        for(size_t l = 1; l <= k; l++) {
            for(size_t v = 0; v < V; v++) {
                if(((m.at(v, v) && candidateCount[v] >= l) || (!m.at(v, v) && candidateCount[v] >= l - 1))
                && rootedHere[v].size() * k >= l * V) {
                    for(auto r:rootedHere[v]) {
                        if (electedCount[r] < l) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    bool issEJRplus(const std::set<size_t>& W, size_t k) {

        // count how many candidates approved by voter were elected
        std::vector<size_t> electedCount(V, 0);
        std::vector<std::vector<size_t>> supporters(V);

        for(size_t v = 0; v < V; v++) {
            for(size_t w:W){
                if (m.at(v, w)) {
                    electedCount[v]++;
                }
            }
            for(size_t c = 0; c < C; c++){
                if (m.at(v, c)) {
                    supporters[c].push_back(v);
                }
            }
        }
        
        for(size_t l = 1; l <= k; l++) {
            for(size_t c = 0; c < C; c++){
                if(W.count(c) == 1) {
                    continue;
                }
                if (supporters[c].size() * k < l * V) { 
                    continue;
                }

                for(auto s:supporters[c]) {
                    if(electedCount[s] < l) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

};