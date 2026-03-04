#pragma once

#include <iostream>
#include <vector>
#include "bitmatrix.h"

class AxiomChecker {
private:
    const BitMatrix& m;
    size_t n;

public:
    AxiomChecker(const BitMatrix& matrix) : m(matrix), n(matrix.size()) {}

    size_t approvedCount(size_t voter, const std::vector<size_t>& W) {
        size_t cnt = 0;
        for (auto c : W)
            if (m.at(voter, c)) cnt++;
        return cnt;
    }

    size_t groupApprovedUnion(const std::vector<size_t>& group, const std::vector<size_t>& W) {
        std::vector<bool> seen(n, false);

        for (auto v : group) {
            for (auto c : W) {
                if (m.at(v, c))
                    seen[c] = true;
            }
        }

        return std::count(seen.begin(), seen.end(), true);
    }

    size_t countVotersApprovingAll(const std::vector<size_t>& S) {
        size_t count = 0;

        for (size_t v = 0; v < n; ++v) {
            bool ok = true;
            for (auto c : S) {
                if (!m.at(v, c)) {
                    ok = false;
                    break;
                }
            }
            if (ok) count++;
        }

        return count;
    }

    size_t maxClaimForVoter(size_t voter, size_t k) {
        std::vector<size_t> Ai;

        for (size_t c = 0; c < n; ++c)
            if (m.at(voter, c))
                Ai.push_back(c);

        size_t best = 0;
        size_t mAi = Ai.size();

        for (size_t mask = 1; mask < (1ULL << mAi); ++mask) {
            std::vector<size_t> S;

            for (size_t i = 0; i < mAi; ++i)
                if (mask & (1ULL << i))
                    S.push_back(Ai[i]);

            size_t sSize = S.size();
            size_t supporters = countVotersApprovingAll(S);

            if (supporters * k >= sSize * n) {
                best = std::max(best, sSize);
            }
        }

        return best;
    }

    bool isJR(const std::vector<size_t>& W, size_t k) {
        size_t quota = n / k;

        for (size_t mask = 1; mask < (1ULL << n); ++mask) {
            std::vector<size_t> group;

            for (size_t i = 0; i < n; ++i)
                if (mask & (1ULL << i))
                    group.push_back(i);

            if (group.size() < quota) continue;

            for (size_t c = 0; c < n; ++c) {
                bool allApprove = true;

                for (auto v : group)
                    if (!m.at(v, c))
                        allApprove = false;

                if (!allApprove) continue;

                bool satisfied = false;

                for (auto v : group)
                    if (approvedCount(v, W) >= 1)
                        satisfied = true;

                if (!satisfied)
                    return false;
            }
        }

        return true;
    }

    bool isPJR(const std::vector<size_t>& W, size_t k) {
        for (size_t ell = 1; ell <= k; ++ell) {
            size_t quota = ell * n / k;

            for (size_t mask = 1; mask < (1ULL << n); ++mask) {
                std::vector<size_t> group;

                for (size_t i = 0; i < n; ++i)
                    if (mask & (1ULL << i))
                        group.push_back(i);

                if (group.size() < quota) continue;

                size_t common = 0;

                for (size_t c = 0; c < n; ++c) {
                    bool ok = true;
                    for (auto v : group)
                        if (!m.at(v, c)) ok = false;
                    if (ok) common++;
                }

                if (common < ell) continue;

                size_t represented = groupApprovedUnion(group, W);

                if (represented < ell)
                    return false;
            }
        }

        return true;
    }

    bool isEJR(const std::vector<size_t>& W, size_t k) {
        for (size_t ell = 1; ell <= k; ++ell) {
            size_t quota = ell * n / k;

            for (size_t mask = 1; mask < (1ULL << n); ++mask) {
                std::vector<size_t> group;

                for (size_t i = 0; i < n; ++i)
                    if (mask & (1ULL << i))
                        group.push_back(i);

                if (group.size() < quota) continue;

                size_t common = 0;

                for (size_t c = 0; c < n; ++c) {
                    bool ok = true;
                    for (auto v : group)
                        if (!m.at(v, c)) ok = false;
                    if (ok) common++;
                }

                if (common < ell) continue;

                bool exists = false;

                for (auto v : group)
                    if (approvedCount(v, W) >= ell)
                        exists = true;

                if (!exists)
                    return false;
            }
        }

        return true;
    }
    bool isEJRplus(const std::vector<size_t>& W, size_t k) {
        std::vector<bool> inW(n, false);
        for (auto c : W) inW[c] = true;

        for (size_t c = 0; c < n; ++c) {
            if (inW[c]) continue;

            for (size_t ell = 1; ell <= k; ++ell) {
                size_t quota = ell * n / k;

                for (size_t mask = 1; mask < (1ULL << n); ++mask) {
                    std::vector<size_t> group;

                    for (size_t i = 0; i < n; ++i)
                        if (mask & (1ULL << i))
                            group.push_back(i);

                    if (group.size() < quota) continue;

                    bool allApproveC = true;
                    for (auto v : group)
                        if (!m.at(v, c))
                            allApproveC = false;

                    if (!allApproveC) continue;

                    bool bad = true;

                    for (auto v : group)
                        if (approvedCount(v, W) >= ell)
                            bad = false;

                    if (bad)
                        return false;
                }
            }
        }

        return true;
    }

    bool isIR(const std::vector<size_t>& W, size_t k) {
        for (size_t v = 0; v < n; ++v) {

            size_t claim = maxClaimForVoter(v, k);

            size_t approvedInW = 0;
            for (auto c : W)
                if (m.at(v, c))
                    approvedInW++;

            if (approvedInW < claim) {
                return false;
            }
        }

        return true;
    }
};