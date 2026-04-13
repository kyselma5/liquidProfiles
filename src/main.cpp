#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <mutex>
#include <thread>
#include <cassert>

#include "bitmatrix.h"
#include "rules.h"
#include "axioms.h"
#include "configLoader.h"

struct Stats{
    double JR_t = 0;
    double PJR_t = 0;
    double EJR_t = 0;
    double EJRplus_t = 0;
    double IR_t = 0;
    double LR_t = 0;
    double PR_t = 0;
    double CS_t = 0;
    double sJR_t = 0;
    double sEJR_t = 0;
    double sEJRplus_t = 0;
};

void printProgressBar(size_t current, size_t total) {
    const int barWidth = 40;

    float progress = (float)current / total;
    int pos = barWidth * progress;

    std::cout << "\r[";
    for (int i = 0; i < barWidth; i++) {
        if (i < pos)
            std::cout << "█";
        else
            std::cout << " ";
    }

    std::cout << "] " << int(progress * 100.0) << "%";
    std::cout.flush();
}

BitMatrix generateDelegatinMatrix(size_t size, std::random_device& rd) {
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, size - 1);
    std::vector<size_t> v(size);
    for (size_t j = 0; j < size; j++) {
        v[j] = dist(gen);
    }

    return BitMatrix(v);
}

void generateRandomDelegationMatrixes(const Config & cfg, std::random_device& rd) {

    std::cout << "GENERATING\n";
    std::ofstream out(cfg.matrixFile);
    out << cfg.numberOfMatrixes << " " << cfg.matrixSize << "\n\n";

    size_t prevPercentage = 0;

    for (size_t i = 0; i < cfg.numberOfMatrixes; i++) {
        if(100*(i+1)/cfg.numberOfMatrixes != prevPercentage) {
            prevPercentage = 100*(i+1)/cfg.numberOfMatrixes;
            printProgressBar(i+1, cfg.numberOfMatrixes);
        }
        generateDelegatinMatrix(cfg.matrixSize, rd).print(out);
        out << "\n";
    }
    std::cout << std::endl;
}

void workWithMatrix(BitMatrix & m, const Config & cfg, std::ostream & os,
                    const std::vector<double> & CC, const std::vector<double> & PAV, Stats & s) {
    auto tc = m.transitiveClosure();
    BitMatrix out(m.size());
    tc.isLiquidProfileFast(out);

    Rules r(tc);
    AxiomChecker a(tc);

    // sinkcount vs rules distances
    // different committee sizes 
    for(const auto committeeSize : cfg.committeeSizes) {
        auto resAV = r.approvalVoting(committeeSize);
        auto resSP = r.sequentialPhragmen(committeeSize);
        auto resCC = r.thieleGreedy(committeeSize, CC);
        auto resPAV = r.thieleGreedy(committeeSize, PAV);
        auto resMES = r.MES(committeeSize);
        auto resGJCR = r.GJCR(committeeSize);

        auto resMES_AV = r.approvalVoting(committeeSize, resMES);
        auto resGJCR_AV = r.approvalVoting(committeeSize, resGJCR);

        os << committeeSize << "," 
        << tc.m_condensed.maxTreeSize << "," 
        << tc.m_condensed.maxCycleSize << "," 
        << tc.m_condensed.maxDistanceToCycle << "," 
        << tc.m_condensed.sinkCount << "," 
        << tc.m_condensed.sourceCount << "," 
        << tc.m_condensed.isolatedCount << "," 
        << hammingDistance(resAV, resSP) << ","
        << hammingDistance(resAV, resCC) << ","
        << hammingDistance(resAV, resPAV) << ","
        << hammingDistance(resAV, resMES_AV) << ","
        << hammingDistance(resAV, resGJCR_AV) << ","

        << hammingDistance(resSP, resCC) << ","
        << hammingDistance(resSP, resPAV) << ","
        << hammingDistance(resSP, resMES_AV) << ","
        << hammingDistance(resSP, resGJCR_AV) << ","

        << hammingDistance(resCC, resPAV) << ","
        << hammingDistance(resCC, resMES_AV) << ","
        << hammingDistance(resCC, resGJCR_AV) << ","

        << hammingDistance(resPAV, resMES_AV) << ","
        << hammingDistance(resPAV, resGJCR_AV) << ","

        << hammingDistance(resMES_AV, resGJCR_AV);

        std::chrono::steady_clock::time_point start;
        bool JR_AV = false, PJR_AV = false, EJR_AV = false, EJRplus_AV = false, IR_AV = false, LR_AV = false, CS_AV = false, PR_AV = false, sJR_AV = false, sEJR_AV = false, sEJRplus_AV = false;
        bool JR_SP = false, PJR_SP = false, EJR_SP = false, EJRplus_SP= false, IR_SP = false, LR_SP = false, CS_SP = false, PR_SP = false, sJR_SP = false, sEJR_SP = false, sEJRplus_SP = false;
        bool JR_CC = false, PJR_CC = false, EJR_CC = false, EJRplus_CC= false, IR_CC = false, LR_CC = false, CS_CC = false, PR_CC = false, sJR_CC = false, sEJR_CC = false, sEJRplus_CC = false;
        bool JR_PAV = false, PJR_PAV = false, EJR_PAV = false, EJRplus_PAV= false, IR_PAV = false, LR_PAV = false, CS_PAV = false, PR_PAV = false, sJR_PAV = false, sEJR_PAV = false, sEJRplus_PAV = false;
        bool JR_MES = false, PJR_MES = false, EJR_MES = false, EJRplus_MES= false, IR_MES = false, LR_MES = false, CS_MES = false, PR_MES = false, sJR_MES = false, sEJR_MES = false, sEJRplus_MES = false;
        bool JR_GJCR = false, PJR_GJCR = false, EJR_GJCR = false, EJRplus_GJCR= false, IR_GJCR = false, LR_GJCR = false, CS_GJCR = false, PR_GJCR = false, sJR_GJCR = false, sEJR_GJCR = false, sEJRplus_GJCR = false;
        
        if(cfg.JR) {
            start = std::chrono::steady_clock::now();
            JR_AV = a.isJR(resAV, committeeSize);
            JR_SP = a.isJR(resSP, committeeSize);
            JR_CC = a.isJR(resCC, committeeSize);
            JR_PAV = a.isJR(resPAV, committeeSize);
            JR_MES = a.isJR(resMES, committeeSize);
            JR_GJCR = a.isJR(resGJCR, committeeSize);
            s.JR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << JR_AV << ","<< JR_SP << ","<< JR_CC << ","<< JR_PAV << ","<< JR_MES << ","<< JR_GJCR;
        }
        if(cfg.PJR) {
            start = std::chrono::steady_clock::now();
            PJR_AV = a.isPJR(resAV, committeeSize);
            PJR_SP = a.isPJR(resSP, committeeSize);
            PJR_CC = a.isPJR(resCC, committeeSize);
            PJR_PAV = a.isPJR(resPAV, committeeSize);
            PJR_MES = a.isPJR(resMES, committeeSize);
            PJR_GJCR = a.isPJR(resGJCR, committeeSize);
            s.PJR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << PJR_AV << ","<< PJR_SP << ","<< PJR_CC << ","<< PJR_PAV << ","<< PJR_MES << ","<< PJR_GJCR;
 
        }
        if(cfg.EJR) {
            start = std::chrono::steady_clock::now();
            EJR_AV = a.isEJR(resAV, committeeSize);
            EJR_SP = a.isEJR(resSP, committeeSize);
            EJR_CC = a.isEJR(resCC, committeeSize);
            EJR_PAV = a.isEJR(resPAV, committeeSize);
            EJR_MES = a.isEJR(resMES, committeeSize);
            EJR_GJCR = a.isEJR(resGJCR, committeeSize);
            s.EJR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << EJR_AV << ","<< EJR_SP << ","<< EJR_CC << ","<< EJR_PAV << ","<< EJR_MES << ","<< EJR_GJCR;
        }
        if(cfg.EJRplus) {
            start = std::chrono::steady_clock::now();
            EJRplus_AV = a.isEJRplus(resAV, committeeSize);
            EJRplus_SP = a.isEJRplus(resSP, committeeSize);
            EJRplus_CC = a.isEJRplus(resCC, committeeSize);
            EJRplus_PAV = a.isEJRplus(resPAV, committeeSize);
            EJRplus_MES = a.isEJRplus(resMES, committeeSize);
            EJRplus_GJCR = a.isEJRplus(resGJCR, committeeSize);
            s.EJRplus_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << EJRplus_AV << ","<< EJRplus_SP << ","<< EJRplus_CC << ","<< EJRplus_PAV << ","<< EJRplus_MES << ","<< EJRplus_GJCR;
        }
        if(cfg.IR) {
            start = std::chrono::steady_clock::now();
            IR_AV = a.isIR(resAV, committeeSize);
            IR_SP = a.isIR(resSP, committeeSize);
            IR_CC = a.isIR(resCC, committeeSize);
            IR_PAV = a.isIR(resPAV, committeeSize);
            IR_MES = a.isIR(resMES, committeeSize);
            IR_GJCR = a.isIR(resGJCR, committeeSize);
            s.IR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << IR_AV << ","<< IR_SP << ","<< IR_CC << ","<< IR_PAV << ","<< IR_MES << ","<< IR_GJCR;
        }
        if(cfg.LR) {
            start = std::chrono::steady_clock::now();
            LR_AV = a.isLR(resAV, committeeSize);
            LR_SP = a.isLR(resSP, committeeSize);
            LR_CC = a.isLR(resCC, committeeSize);
            LR_PAV = a.isLR(resPAV, committeeSize);
            LR_MES = a.isLR(resMES, committeeSize);
            LR_GJCR = a.isLR(resGJCR, committeeSize);
            s.LR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << LR_AV << ","<< LR_SP << ","<< LR_CC << ","<< LR_PAV << ","<< LR_MES << ","<< LR_GJCR;
        }
        if(cfg.CS) {
            start = std::chrono::steady_clock::now();
            CS_AV = a.isCS(resAV, committeeSize);
            CS_SP = a.isCS(resSP, committeeSize);
            CS_CC = a.isCS(resCC, committeeSize);
            CS_PAV = a.isCS(resPAV, committeeSize);
            CS_MES = a.isCS(resMES, committeeSize);
            CS_GJCR = a.isCS(resGJCR, committeeSize);
            s.CS_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << CS_AV << ","<< CS_SP << ","<< CS_CC << ","<< CS_PAV << ","<< CS_MES << ","<< CS_GJCR;
        }
        if(cfg.PR) {
            start = std::chrono::steady_clock::now();
            PR_AV = a.isPR(resAV, committeeSize);
            PR_SP = a.isPR(resSP, committeeSize);
            PR_CC = a.isPR(resCC, committeeSize);
            PR_PAV = a.isPR(resPAV, committeeSize);
            PR_MES = a.isPR(resMES, committeeSize);
            PR_GJCR = a.isPR(resGJCR, committeeSize);
            s.PR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << PR_AV << ","<< PR_SP << ","<< PR_CC << ","<< PR_PAV << ","<< PR_MES << ","<< PR_GJCR;
        }
        if(cfg.sJR) {
            start = std::chrono::steady_clock::now();
            sJR_AV = a.issJR(resAV, committeeSize);
            sJR_SP = a.issJR(resSP, committeeSize);
            sJR_CC = a.issJR(resCC, committeeSize);
            sJR_PAV = a.issJR(resPAV, committeeSize);
            sJR_MES = a.issJR(resMES, committeeSize);
            sJR_GJCR = a.issJR(resGJCR, committeeSize);
            s.sJR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << sJR_AV << ","<< sJR_SP << ","<< sJR_CC << ","<< sJR_PAV << ","<< sJR_MES << ","<< sJR_GJCR;
        }
        if(cfg.sEJR) {
            start = std::chrono::steady_clock::now();
            sEJR_AV = a.issEJR(resAV, committeeSize);
            sEJR_SP = a.issEJR(resSP, committeeSize);
            sEJR_CC = a.issEJR(resCC, committeeSize);
            sEJR_PAV = a.issEJR(resPAV, committeeSize);
            sEJR_MES = a.issEJR(resMES, committeeSize);
            sEJR_GJCR = a.issEJR(resGJCR, committeeSize);
            s.sEJR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << sEJR_AV << ","<< sEJR_SP << ","<< sEJR_CC << ","<< sEJR_PAV << ","<< sEJR_MES << ","<< sEJR_GJCR;
        }
        if(cfg.sEJRplus) {
            start = std::chrono::steady_clock::now();
            sEJRplus_AV = a.issEJRplus(resAV, committeeSize);
            sEJRplus_SP = a.issEJRplus(resSP, committeeSize);
            sEJRplus_CC = a.issEJRplus(resCC, committeeSize);
            sEJRplus_PAV = a.issEJRplus(resPAV, committeeSize);
            sEJRplus_MES = a.issEJRplus(resMES, committeeSize);
            sEJRplus_GJCR = a.issEJRplus(resGJCR, committeeSize);
            s.sEJRplus_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
            os << "," << sEJRplus_AV << ","<< sEJRplus_SP << ","<< sEJRplus_CC << ","<< sEJRplus_PAV << ","<< sEJRplus_MES << ","<< sEJRplus_GJCR;
        }

        os << "\n";

        if(cfg.doSanityCheck) {
            if(cfg.sEJRplus && cfg.sEJR) {
                assert(sEJRplus_AV == sEJR_AV);
                assert(sEJRplus_SP == sEJR_SP);
                assert(sEJRplus_CC == sEJR_CC);
                assert(sEJRplus_PAV == sEJR_PAV);
                assert(sEJRplus_MES == sEJR_MES);
                assert(sEJRplus_GJCR == sEJR_GJCR);
            }
            if(cfg.sEJRplus && cfg.IR) {
                assert(sEJRplus_AV == IR_AV);
                assert(sEJRplus_SP == IR_SP);
                assert(sEJRplus_CC == IR_CC);
                assert(sEJRplus_PAV == IR_PAV);
                assert(sEJRplus_MES == IR_MES);
                assert(sEJRplus_GJCR == IR_GJCR);
            }
            if(cfg.sEJRplus && cfg.LR) {
                assert(sEJRplus_AV == LR_AV);
                assert(sEJRplus_SP == LR_SP);
                assert(sEJRplus_CC == LR_CC);
                assert(sEJRplus_PAV == LR_PAV);
                assert(sEJRplus_MES == LR_MES);
                assert(sEJRplus_GJCR == LR_GJCR);
            }
            if(cfg.sEJR && cfg.sJR) {
                assert(!sEJRplus_AV || sJR_AV);
                assert(!sEJRplus_SP || sJR_SP);
                assert(!sEJRplus_CC || sJR_CC);
                assert(!sEJRplus_PAV || sJR_PAV);
                assert(!sEJRplus_MES || sJR_MES);
                assert(!sEJRplus_GJCR || sJR_GJCR);
            }
            if(cfg.PR && cfg.sJR) {
                assert(!PR_AV || sJR_AV);
                assert(!PR_SP || sJR_SP);
                assert(!PR_CC || sJR_CC);
                assert(!PR_PAV || sJR_PAV);
                assert(!PR_MES || sJR_MES);
                assert(!PR_GJCR || sJR_GJCR);
            }
            if(cfg.CS && cfg.EJR) {
                assert(!CS_AV || EJR_AV);
                assert(!CS_SP || EJR_SP);
                assert(!CS_CC || EJR_CC);
                assert(!CS_PAV || EJR_PAV);
                assert(!CS_MES || EJR_MES);
                assert(!CS_GJCR || EJR_GJCR);
            }
            if(cfg.sEJR && cfg.EJR) {
                assert(!sEJR_AV || EJR_AV);
                assert(!sEJR_SP || EJR_SP);
                assert(!sEJR_CC || EJR_CC);
                assert(!sEJR_PAV || EJR_PAV);
                assert(!sEJR_MES || EJR_MES);
                assert(!sEJR_GJCR || EJR_GJCR);
            }
            if(cfg.EJR && cfg.JR) {
                assert(!EJR_AV || JR_AV);
                assert(!EJR_SP || JR_SP);
                assert(!EJR_CC || JR_CC);
                assert(!EJR_PAV || JR_PAV);
                assert(!EJR_MES || JR_MES);
                assert(!EJR_GJCR || JR_GJCR);
            }
            if(cfg.sJR && cfg.JR) {
                assert(!sJR_AV || JR_AV);
                assert(!sJR_SP || JR_SP);
                assert(!sJR_CC || JR_CC);
                assert(!sJR_PAV || JR_PAV);
                assert(!sJR_MES || JR_MES);
                assert(!sJR_GJCR || JR_GJCR);
            }
            if(cfg.JR && cfg.PJR) {
                assert(!PJR_AV || JR_AV);
                assert(!PJR_SP || JR_SP);
                assert(!PJR_CC || JR_CC);
                assert(!PJR_PAV || JR_PAV);
                assert(!PJR_MES || JR_MES);
                assert(!PJR_GJCR || JR_GJCR);
            }     
            if(cfg.PJR && cfg.EJR) {
                assert(!EJR_AV || PJR_AV);
                assert(!EJR_SP || PJR_SP);
                assert(!EJR_CC || PJR_CC);
                assert(!EJR_PAV || PJR_PAV);
                assert(!EJR_MES || PJR_MES);
                assert(!EJR_GJCR || PJR_GJCR);
            }
        }
    }
}

void processMatricesFromFile(const Config & cfg){
    std::ifstream in(cfg.matrixFile);
    std::ofstream os(cfg.resultsFile);

    std::cout << "PROCESS" << std::endl;

    // --- CSV HEADER ---
    os << "committeeSize,maxTreeSize,maxCycleSize,maxDistanceToCycle,"
       << "sinkCount,sourceCount,isolatedCount,"
       << "hammingDist_AV_to_SP,hammingDist_AV_to_CC,hammingDist_AV_to_PAV,hammingDist_AV_to_MES_AV,hammingDist_AV_to_GJCR_AV,"
       << "hammingDist_SP_to_CC,hammingDist_SP_to_PAV,hammingDist_SP_to_MES_AV,hammingDist_SP_to_GJCR_AV,"
       << "hammingDist_CC_to_PAV,hammingDist_CC_to_MES_AV,hammingDist_CC_to_GJCR_AV,"
       << "hammingDist_PAV_to_MES_AV,hammingDist_PAV_to_GJCR_AV,"
       << "hammingDist_MES_AV_to_GJCR_AV";
       
    if(cfg.JR) {os << ",JR_AV,JR_SP,JR_CC,JR_PAV,JR_MES,JR_GJCR";}
    if(cfg.PJR) {os << ",PJR_AV,PJR_SP,PJR_CC,PJR_PAV,PJR_MES,PJR_GJCR";}
    if(cfg.EJR) {os << ",EJR_AV,EJR_SP,EJR_CC,EJR_PAV,EJR_MES,EJR_GJCR";}
    if(cfg.EJRplus) {os << ",EJR+_AV,EJR+_SP,EJR+_CC,EJR+_PAV,EJR+_MES,EJR+_GJCR";}
    if(cfg.IR) {os << ",IR_AV,IR_SP,IR_CC,IR_PAV,IR_MES,IR_GJCR";}
    if(cfg.LR) {os << ",LR_AV,LR_SP,LR_CC,LR_PAV,LR_MES,LR_GJCR";}
    if(cfg.CS) {os << ",CS_AV,CS_SP,CS_CC,CS_PAV,CS_MES,CS_GJCR";}
    if(cfg.PR) {os << ",PR_AV,PR_SP,PR_CC,PR_PAV,PR_MES,PR_GJCR";}
    if(cfg.sJR) {os << ",sJR_AV,sJR_SP,sJR_CC,sJR_PAV,sJR_MES,sJR_GJCR";}
    if(cfg.sEJR) {os << ",sEJR_AV,sEJR_SP,sEJR_CC,sEJR_PAV,sEJR_MES,sEJR_GJCR";}
    if(cfg.sEJRplus) {os << ",sEJR+_AV,sEJR+_SP,sEJR+_CC,sEJR+_PAV,sEJR+_MES,sEJR+_GJCR";}
    os << "\n";

    // --- INPUT HEADER ---
    size_t count, size;
    in >> count >> size;
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // --- PRECOMPUTE ---
    std::vector<double> CC(size, 0), PAV(size, 0);
    CC[0] = 1;
    PAV[0] = 1;
    for (size_t i = 1; i < size; i++) {
        PAV[i] = 1.0 / (i+1);
    }

    std::vector<std::string> outputs(count);

    // --- QUEUE + SYNC ---
    std::queue<std::pair<size_t, std::string>> queue;
    std::mutex mtx;
    std::condition_variable cv_not_empty, cv_not_full;

    const size_t MAX_QUEUE_SIZE = 50;
    bool doneReading = false;

    // --- TIMING ---
    std::mutex timeMutex;
    Stats s;

    // --- THREADS ---
    size_t numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;

    std::mutex coutMutex;

    // =========================
    // WORKERS
    // =========================
    for (size_t t = 0; t < numThreads; t++) {
        workers.emplace_back([&]() {
            Stats ls;

            while (true) {
                std::pair<size_t, std::string> task;

                {
                    std::unique_lock lock(mtx);
                    cv_not_empty.wait(lock, [&]() {
                        return !queue.empty() || doneReading;
                    });

                    if (queue.empty() && doneReading)
                        break;

                    task = std::move(queue.front());
                    queue.pop();

                    lock.unlock();
                    cv_not_full.notify_one();
                }

                BitMatrix m(task.second);
                std::stringstream ss;

                workWithMatrix(m, cfg, ss, CC, PAV, ls);

                outputs[task.first] = ss.str();
            }

            std::lock_guard lock(timeMutex);
            s.JR_t += ls.JR_t;
            s.PJR_t += ls.PJR_t;
            s.EJR_t += ls.EJR_t;
            s.EJRplus_t += ls.EJRplus_t;
            s.IR_t += ls.IR_t;
            s.PR_t += ls.PR_t;
            s.CS_t += ls.CS_t;
            s.sJR_t += ls.sJR_t;
            s.sEJR_t += ls.sEJR_t;
            s.sEJRplus_t += ls.sEJRplus_t;
        });
    }

    // =========================
    // READER THREAD
    // =========================
    std::thread reader([&]() {
        std::string line;
        size_t prevPercentage = 0;

        for (size_t i = 0; i < count; i++) {
            if(100*(i+1)/count != prevPercentage) {
                prevPercentage = 100*(i+1)/count;
                printProgressBar(i+1, count);
            }
            std::stringstream matrixStream;

            for (size_t r = 0; r < size; r++) {
                std::getline(in, line);
                matrixStream << line << "\n";
            }
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            std::unique_lock lock(mtx);
            cv_not_full.wait(lock, [&]() {
                return queue.size() < MAX_QUEUE_SIZE;
            });

            queue.emplace(i, matrixStream.str());

            lock.unlock();
            cv_not_empty.notify_one();
        }

        {
            std::lock_guard lock(mtx);
            doneReading = true;
        }
        cv_not_empty.notify_all();
    });

    // =========================
    // JOIN
    // =========================
    reader.join();
    for (auto &t : workers)
        t.join();

    // =========================
    // OUTPUT
    // =========================
    for (const auto& line : outputs) {
        os << line;
    }

    std::cout << "\nDONE\n";

    auto totalTime = (s.JR_t+s.PJR_t+s.EJR_t+s.EJRplus_t+s.IR_t+s.LR_t+s.PR_t+s.CS_t+s.sJR_t+s.sEJR_t+s.sEJRplus_t);

    std::cout << "\n\nTotal axiom time - " << totalTime/1000 << std::endl;
    totalTime /= 100;

    std::cout << s.JR_t/totalTime << " % JR_time\n";
    std::cout << s.PJR_t/totalTime << " % PJR_time\n";
    std::cout << s.EJR_t/totalTime << " % EJR_time\n";
    std::cout << s.EJRplus_t/totalTime << " % EJRplus_time\n";
    std::cout << s.IR_t/totalTime << " % IR_time\n";
    std::cout << s.LR_t/totalTime << " % LR_time\n";
    std::cout << s.PR_t/totalTime << " % PR_time\n";
    std::cout << s.CS_t/totalTime << " % CS_time\n";
    std::cout << s.sJR_t/totalTime << " % SJR_time\n";
    std::cout << s.sEJR_t/totalTime << " % SEJR_time\n";
    std::cout << s.sEJRplus_t/totalTime << " % SEJRplus_time\n";
}

bool endsWith(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

int main(int argc, char* argv[]) {

    if(argc != 2) {
        std::cout << "specify path to configuration file\n";
        return 1;
    }

    Config cfg = loadConfig(argv[1]);

    printConfig(cfg);

    if (cfg.matrixFile.empty()) {
        std::cout << "no file to read matrixes from or generate matrixes to\n";
        return 1;
    }

    if(cfg.generate) {
        if(cfg.numberOfMatrixes == 0) {
            std::cout << "number of matrixes to generate not set [-n]\n";
            return 1;
        }
        if(cfg.matrixSize == 0) {
            std::cout << "size of matrixes to generate not set [-s]\n";
            return 1;
        }
        std::random_device rd;
        generateRandomDelegationMatrixes(cfg, rd);
    }
    if(cfg.process) {
        if(cfg.resultsFile == "") {
            std::cout << "no file to save results was given (needs to be .csv file)\n";
            return 1;
        }
        if(cfg.committeeSizes.size() == 0) {
            std::cout << "no elected committee not set [-c]\n";
            return 1;
        }
        processMatricesFromFile(cfg);
    }
    
    return 0;
}
