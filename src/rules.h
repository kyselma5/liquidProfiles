#pragma once

#include <iostream>
#include <vector>
#include <set>
#include "bitmatrix.h"

size_t hammingDistance(const std::set<size_t>& A, const std::set<size_t>& B) {
    size_t distance = 0;

    auto itA = A.begin();
    auto itB = B.begin();

    while (itA != A.end() && itB != B.end()) {
        if (*itA < *itB) {
            distance++;
            ++itA;
        }
        else if (*itB < *itA) {
            distance++;
            ++itB;
        }
        else {
            ++itA;
            ++itB;
        }
    }

    distance += std::distance(itA, A.end());
    distance += std::distance(itB, B.end());

    return distance;
}

class Rules
{
private:
    BitMatrix m_matrix;
    size_t n;
public:
    Rules(const BitMatrix & matrix): m_matrix(matrix), n(matrix.size()){}

    std::set<size_t> approvalVoting(size_t committeeSize, std::set<size_t> winners = {}) {
        std::vector<size_t> candidateScores(n, 0);
        for(size_t i = 0; i < n; ++i){
            for(size_t j = 0; j < n; ++j){
                if(m_matrix.at(j, i)) {
                    candidateScores[i]++;
                }
            }
        }
        std::vector<size_t> res(n);

        for (size_t i = 0; i < n; ++i)
            res[i] = i;

        std::sort(res.begin(), res.end(), [&](size_t a, size_t b) {
            return candidateScores[a] > candidateScores[b];
        });

        for (auto w:res) {
            if(winners.size() >= committeeSize){
                break;
            }
            winners.insert(w);
        }
        return winners;
    }

    std::set<size_t> sequentialPhragmen(size_t committeeSize) {
        std::vector<double> budget(n, 0.0);
        std::vector<bool> selected(n, false);
        std::set<size_t> winners;

        while (winners.size() < committeeSize) {
            double bestTime = std::numeric_limits<double>::max();
            size_t bestCandidate = n;

            for (size_t c = 0; c < n; ++c) {
                if (selected[c]) continue;

                double sumBudget = 0.0;
                size_t supporters = 0;

                for (size_t v = 0; v < n; ++v) {
                    if (m_matrix.at(v, c)) {
                        sumBudget += budget[v];
                        supporters++;
                    }
                }

                if (supporters == 0) continue;

                double t = (1.0 - sumBudget) / supporters;

                if (t < bestTime) {
                    bestTime = t;
                    bestCandidate = c;
                }
            }

            for (size_t v = 0; v < n; ++v) {
                budget[v] += bestTime;
            }

            winners.insert(bestCandidate);
            selected[bestCandidate] = true;

            for (size_t v = 0; v < n; ++v) {
                if (m_matrix.at(v, bestCandidate)) {
                    budget[v] = 0.0;
                }
            }
        }

        return winners;
    }

    double computeScore(const std::set<size_t>& committee, const std::vector<double>& w) {
        double totalScore = 0.0;

        for (size_t v = 0; v < n; ++v) {
            size_t approved = 0;

            for (size_t c : committee) {
                if (m_matrix.at(v, c)) {
                    approved++;
                }
            }

            for (size_t i = 0; i < approved; ++i) {
                totalScore += w[i];
            }
        }

        return totalScore;
    }

    void generate(size_t start, size_t depth, size_t committeeSize, std::set<size_t>& current,
                    std::set<size_t>& best, double& bestScore, const std::vector<double>& w) {

        if (depth == committeeSize) {
            double score = computeScore(current, w);
            if (score > bestScore) {
                bestScore = score;
                best = current;
            }
            return;
        }

        for (size_t i = start; i < n; ++i) {
            current.insert(i);
            generate(i + 1, depth + 1, committeeSize, current, best, bestScore, w);
            current.erase(i);
        }
    }

    std::set<size_t> thieleRule(size_t committeeSize, const std::vector<double>& w) {

        std::set<size_t> current, best;
        double bestScore = -1;

        generate(0, 0, committeeSize, current, best, bestScore, w);

        return best;
    }

    std::set<size_t> thieleGreedy(size_t committeeSize, const std::vector<double>& w, std::set<size_t> winners = {}) {

        std::vector<size_t> approvedCount(n, 0);
        std::vector<bool> used(n, false);

        for (size_t step = 0; step < committeeSize; ++step) {
            double bestGain = -1;
            size_t bestCandidate = n;

            for (size_t c = 0; c < n; ++c) {
                if (used[c]) continue;

                double gain = 0.0;

                for (size_t v = 0; v < n; ++v) {
                    if (m_matrix.at(v, c)) {
                        gain += w[approvedCount[v]];
                    }
                }

                if (gain > bestGain) {
                    bestGain = gain;
                    bestCandidate = c;
                }
            }

            winners.insert(bestCandidate);
            used[bestCandidate] = true;

            for (size_t v = 0; v < n; ++v) {
                if (m_matrix.at(v, bestCandidate)) {
                    approvedCount[v]++;
                }
            }
        }

        return winners;
    }

    double computeRho(size_t c, const std::vector<double>& budget) {
        std::vector<double> supporters;

        for (size_t v = 0; v < n; ++v) {
            if (m_matrix.at(v, c)) {
                supporters.push_back(budget[v]);
            }
        }

        if (supporters.empty()) {
            return std::numeric_limits<double>::infinity();
        }

        std::sort(supporters.begin(), supporters.end());

        double prefixSum = 0.0;
        size_t m = supporters.size();

        for (size_t i = 0; i < m; ++i) {
            double rho = (1.0 - prefixSum) / (m - i);

            if (rho <= supporters[i]) {
                return rho;
            }

            prefixSum += supporters[i];
        }

        return std::numeric_limits<double>::infinity();
    }

    std::set<size_t> MES(size_t committeeSize) {

        std::vector<double> budget(n, double(committeeSize) / n);
        std::vector<bool> selected(n, false);
        std::set<size_t> winners;

        for (size_t step = 0; step < committeeSize; ++step) {
            double bestRho = std::numeric_limits<double>::infinity();
            size_t bestCandidate = n;

            for (size_t c = 0; c < n; ++c) {
                if (selected[c]) continue;

                double rho = computeRho(c, budget);

                if (rho < bestRho) {
                    bestRho = rho;
                    bestCandidate = c;
                }
            }
            if (bestCandidate == n) {
                break;
            }

            winners.insert(bestCandidate);
            selected[bestCandidate] = true;

            for (size_t v = 0; v < n; ++v) {
                if (m_matrix.at(v, bestCandidate)) {
                    budget[v] -= std::min(budget[v], bestRho);
                }
            }
        }

        return winners;
    }

    std::set<size_t> GJCR(size_t committeeSize) {

        std::set<size_t> winners;
        std::vector<bool> selected(n, false);
        std::vector<size_t> approvedCount(n, 0);

        for (size_t l = committeeSize; l >= 1; --l) {

            bool added = true;

            while (added) {
                added = false;

                for (size_t c = 0; c < n; ++c) {
                    if (selected[c]) continue;

                    size_t count = 0;

                    for (size_t v = 0; v < n; ++v) {
                        if (m_matrix.at(v, c) && approvedCount[v] < l) {
                            count++;
                        }
                    }

                    if (count * committeeSize >= l * n) {
                        winners.insert(c);
                        selected[c] = true;

                        for (size_t v = 0; v < n; ++v) {
                            if (m_matrix.at(v, c)) {
                                approvedCount[v]++;
                            }
                        }

                        added = true;
                        break;
                    }
                }
            }
        }
        return winners;
    }
};
