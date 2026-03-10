#pragma once

#include <iostream>
#include <vector>
#include "bitmatrix.h"
#include "other/maxflow2.h"

template<typename T>
void printVector(const std::vector<T> & v) { 
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

    bool isJRFast(const std::vector<size_t>& W, size_t k) {
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

    bool PJRhelper(const std::vector<size_t>& committee, size_t k, size_t maxCandidate, size_t l, std::vector<bool> voters) {
        
        if (l > k) {
            // we already checked all (k and less)-cohesive groups 
            return true;
        }
        //count voters
        size_t voterCount = 0;
        for(bool b : voters) {
            if(b) {
                voterCount++;
            } 
        }
        // pruning the non l-cohesive groups
        if(voterCount * k < l * V){
            return true;
        }

        if (l > 0){
            return false;
            //TODO check the largest group of voters supporting just l-1 winners is big enough
            // probably do it by recursion through all supsets of W up to l-1 size
        }

        // recursive check for smaller groups of voters
        for(size_t c = maxCandidate; c < C; c++) {
            std::vector<bool> votersNew = voters;
            for (size_t v = 0; v < V; v++) {
                votersNew[v] = votersNew[v] && m.at(v, c);
            }
            if (!PJRhelper(committee, k, c+1, l+1, votersNew)) {
                return false;
            }
        }
        return true;
    }

    bool isPJRFast(const std::vector<size_t>& W, size_t k) {
        std::vector<bool> voters(V, true);
        return PJRhelper(W, k, 0, 1, voters);
    }

    bool EJRhelper(const std::vector<size_t>& W, size_t k, size_t maxCandidate, size_t l, std::vector<bool> voters) {
        if (l > k) {
            // we already checked all (k and less)-cohesive groups 
            return true;
        }
        if (l > 0){
            //count voters
            size_t voterCount = 0;
            for(bool b : voters) {
                if(b) {
                    voterCount++;
                } 
            }
            // pruning the non l-cohesive groups
            if(voterCount * k < l * V){
                return true;
            }

            // count how many voters are not represented by at least l elected 
            size_t countNotRepresented = 0;
            for(size_t v = 0; v < V; v++){
                if(voters[v]){
                    size_t countRepresentors = 0;
                    for (size_t w:W) {
                        if(m.at(v,w)) {
                            countRepresentors++;
                            if(countRepresentors == l) {
                                break;
                            }
                        }
                    }
                    if(countRepresentors < l){
                        countNotRepresented++;
                    }
                }
            }

            //check if we have just found l cohesive and not enough represented group.
            if(countNotRepresented * k > l * V){
                return false;
            }
        }

        // recursive check for l+1 cohesive groups of voters
        for(size_t c = maxCandidate; c < C; c++) {

            std::vector<bool> votersNew = voters;

            for (size_t v = 0; v < V; v++) {
                votersNew[v] = votersNew[v] && m.at(v, c);
            }

            if (!EJRhelper(W, k, c+1, l+1, votersNew)) {
                return false;
            }
        }
        return true;
    }

    bool isEJRFast(const std::vector<size_t>& W, size_t k) {
        std::vector<bool> voters(V, true);
        return EJRhelper(W, k, 0, 0, voters);
    }

    bool isEJRplusFast(const std::vector<size_t>& W, size_t k) {
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

    size_t IRhelper(size_t s, size_t k, size_t maxCandidate, const std::vector<bool> & NS, const std::vector<bool>& Ai) {
        
        size_t res = s;

        for(size_t c = maxCandidate; c < C; c++) {

            if(!Ai[c]){
                continue;
            }

            std::vector<bool> newNS = NS;
            size_t countNewNS = 0;
            for(size_t v = 0; v < V; v++){
                if (!m.at(v,c)){
                    newNS[v] = false;
                }
                countNewNS += newNS[v];
            }
            if(countNewNS * k < (s+1) * V) {
                continue;
            }
            res = std::max(res, IRhelper(s+1, k, c+1, newNS, Ai));
        }
        return res;
    }

    bool isIRFast(const std::vector<size_t>& W, size_t k) {
        for(size_t v = 0; v < V; v++) {

            std::vector<bool> Ai(C, false);
            for(size_t c = 0; c < C; c++){
                if(m.at(v, c)){
                    Ai[c] = true;
                }
            }

            size_t representCount = 0;
            for(size_t w:W) {
                if(m.at(v, w)) {
                    representCount++;
                }
            }
            std::vector<bool> NS(V, true);
            if(representCount < IRhelper(0, k, 0, NS, Ai)){
                return false;
            }
        }
        return true;
    }

    bool isPRFast(const std::vector<size_t>& W, size_t k) {
        if(V%k != 0) {
            std::cout << "-";
            return false;
        }
        max_flow<1000, size_t> network;

        for(size_t v = 0; v < V; v++){
            network.add(0, v+1, 1);
            for(size_t w:W){
                if(m.at(v, w)){
                    network.add(1+v, 1+V+w, 1);
                }
            }
        }
        for(size_t w:W){
            network.add(1+V+w, 1+V+k, V/k);
        }

        return (V == network.calc(0, 1+V+k));
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

    bool isCS(const std::vector<size_t>& W, size_t k) {
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

    bool isSJRFast(const std::vector<size_t>& W, size_t k) {

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

    bool isLRFast(const std::vector<size_t>& W, size_t k) {

        for(size_t c = 0; c < C; c++) {

            // count supporters of candidate c
            size_t supportersCount = 0;
            for(size_t v = 0; v < V; v++){
                supportersCount += m.at(v, c);
            }

            // this should be correct (+-1) TODO discuss this against the definition.
            size_t l = (supportersCount * k) / V;
            
            std::vector<size_t> cAndUpC;
            cAndUpC.push_back(c);
            for(size_t s = 0; s < C; s++) {
                if(m.at(c, s)){
                    cAndUpC.push_back(s);
                }
            }

            if (cAndUpC.size() < l) {
                for(size_t e:cAndUpC){
                    bool found = false;
                    for(size_t w:W){
                        if(e == w){
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        return false;
                    }
                }
            }
            else {
                size_t foundCount = 0;
                for(size_t e:cAndUpC){
                    bool found = false;
                    for(size_t w:W){
                        if(e == w){
                            found = true;
                            break;
                        }
                    }
                    if(found){
                        foundCount++;
                    }
                    if(foundCount >= l){
                        break;
                    }
                }
                if (foundCount < l) {
                    return false;
                }
            }
        }
        return true;
    }

    bool SEJRhelper(const std::vector<size_t>& W, size_t k, size_t maxCandidate, size_t l, std::vector<bool> voters) {
        if (l > k) {
            // we already checked all (k and less)-cohesive groups 
            return true;
        }
        if (l > 0){
            //count voters
            size_t voterCount = 0;
            for(bool b : voters) {
                if(b) {
                    voterCount++;
                } 
            }
            // pruning the non l-cohesive groups
            if(voterCount * k < l * V) {
                return true;
            }

            // group is l cohesive, so it should hold, that everyone has at least l candidates in W
            else {
                for(size_t v = 0; v < V; v++) {
                    if(voters[v]) {
                        size_t count = 0;
                        for(size_t w:W){
                            if(m.at(v, w)){
                                count++;
                            }
                        }
                        if(count < l){
                            return false;
                        }
                    }
                }
            }
        }

        // recursive check for l+1 cohesive groups of voters
        for(size_t c = maxCandidate; c < C; c++) {

            std::vector<bool> votersNew = voters;

            for (size_t v = 0; v < V; v++) {
                votersNew[v] = votersNew[v] && m.at(v, c);
            }

            if (!SEJRhelper(W, k, c+1, l+1, votersNew)) {
                return false;
            }
        }
        return true;
    }

    bool isSEJRFast(const std::vector<size_t>& W, size_t k) {
        std::vector<bool> voters(V, true);
        return SEJRhelper(W, k, 0, 0, voters);
    }

    bool isSEJRPlusFast(const std::vector<size_t>& W, size_t k) {

        // count how many candidates approved by voter were elected
        std::vector<size_t> electedCount(V, 0);
        for(size_t v = 0; v < V; v++) {
            for(size_t w:W){
                if (m.at(v, w)) {
                    electedCount[v]++;
                }
            }
        }

        for(size_t c = 0; c < C; c++){
            // filtering candidates in W
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

            // find supporters for this candidate
            std::vector<size_t> supporters;
            for(size_t v = 0; v < V; v++){
                if(m.at(v, c)) {
                    supporters.push_back(v);
                }
            }
            
            // this should be correct (+-1) TODO discuss this against the definition.
            size_t l = (supporters.size() * k) / V;

            for(size_t s : supporters){
                if (electedCount[s] < l) {
                    return false;
                }
            }
        }
        return true;
    }

};