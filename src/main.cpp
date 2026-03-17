#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <mutex>
#include <thread>
#include "bitmatrix.h"
#include "rules.h"
#include "axioms.h"


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

void generateRandomDelegationMatrixes(size_t count, size_t size, std::random_device& rd, const std::string& fileName) {

    std::cout << "GENERATING\n";
    std::ofstream out(fileName);
    out << count << " " << size << "\n\n";

    size_t prevPercentage = 0;

    for (size_t i = 0; i < count; i++) {
        if(100*i/count != prevPercentage) {
            prevPercentage = 100*i/count;
            printProgressBar(i+1, count);
        }
        generateDelegatinMatrix(size, rd).print(out);
        out << "\n";
    }
    std::cout << std::endl;
}

void workWithMatrix(BitMatrix & m, size_t committeeSize, std::ostream & os,
                    const std::vector<double> & CC, const std::vector<double> & PAV,
                    double &JR_t, double &EJR_t, double &EJRplus_t,
                    double &PR_t, double &CS_t, double &SJR_t, double &SEJRplus_t) {
    auto tc = m.transitiveClosure();
    BitMatrix out(m.size());
    tc.isLiquidProfileFast(out);

    Rules r(tc);
    AxiomChecker a(tc);

    auto resAV = r.approvalVoting(committeeSize);
    auto resSP = r.sequentialPhragmen(committeeSize);
    auto resCC = r.thieleGreedy(committeeSize, CC);
    auto resPAV = r.thieleGreedy(committeeSize, PAV);
    auto resMES = r.MES(committeeSize);
    auto resGJCR = r.GJCR(committeeSize);

    auto resSP_AV = r.approvalVoting(committeeSize, resSP);
    auto resCC_AV = r.approvalVoting(committeeSize, resCC);
    auto resPAV_AV = r.approvalVoting(committeeSize, resPAV);
    auto resMES_AV = r.approvalVoting(committeeSize, resMES);
    auto resGJCR_AV = r.approvalVoting(committeeSize, resGJCR);

    os << tc.m_condensed.maxTreeSize << "," 
       << tc.m_condensed.maxCycleSize << "," 
       << tc.m_condensed.maxDistanceToCycle << "," 
       << tc.m_condensed.sinkCount << "," 
       << tc.m_condensed.sourceCount << "," 
       << tc.m_condensed.isolatedCount << "," 
       << hammingDistance(resAV, resSP_AV) << ","
       << hammingDistance(resAV, resCC_AV) << ","
       << hammingDistance(resAV, resPAV_AV) << ","
       << hammingDistance(resAV, resMES_AV) << ","
       << hammingDistance(resAV, resGJCR_AV) << ",";
    auto start = std::chrono::steady_clock::now();
    os << a.isJR(resAV, committeeSize) << ","
       << a.isJR(resSP, committeeSize) << ","
       << a.isJR(resCC, committeeSize) << ","
       << a.isJR(resPAV, committeeSize) << ","
       << a.isJR(resMES, committeeSize) << ","
       << a.isJR(resGJCR, committeeSize) << ",";
    JR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
    os << a.isEJR(resAV, committeeSize) << ","
       << a.isEJR(resSP, committeeSize) << ","
       << a.isEJR(resCC, committeeSize) << ","
       << a.isEJR(resPAV, committeeSize) << ","
       << a.isEJR(resMES, committeeSize) << ","
       << a.isEJR(resGJCR, committeeSize) << ",";
    EJR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
    os << a.isEJRplus(resAV, committeeSize) << ","
       << a.isEJRplus(resSP, committeeSize) << ","
       << a.isEJRplus(resCC, committeeSize) << ","
       << a.isEJRplus(resPAV, committeeSize) << ","
       << a.isEJRplus(resMES, committeeSize) << ","
       << a.isEJRplus(resGJCR, committeeSize) << ",";
    EJRplus_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
    os << a.isPR(resAV, committeeSize) << ","
       << a.isPR(resSP, committeeSize) << ","
       << a.isPR(resCC, committeeSize) << ","
       << a.isPR(resPAV, committeeSize) << ","
       << a.isPR(resMES, committeeSize) << ","
       << a.isPR(resGJCR, committeeSize) << ",";
    PR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
    os << a.isCS(resAV, committeeSize) << ","
       << a.isCS(resSP, committeeSize) << ","
       << a.isCS(resCC, committeeSize) << ","
       << a.isCS(resPAV, committeeSize) << ","
       << a.isCS(resMES, committeeSize) << ","
       << a.isCS(resGJCR, committeeSize) << ",";
    CS_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
    os << a.isSJR(resAV, committeeSize) << ","
       << a.isSJR(resSP, committeeSize) << ","
       << a.isSJR(resCC, committeeSize) << ","
       << a.isSJR(resPAV, committeeSize) << ","
       << a.isSJR(resMES, committeeSize) << ","
       << a.isSJR(resGJCR, committeeSize) << ",";
    SJR_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
    os << a.isSEJRPlus(resAV, committeeSize) << ","
       << a.isSEJRPlus(resSP, committeeSize) << ","
       << a.isSEJRPlus(resCC, committeeSize) << ","
       << a.isSEJRPlus(resPAV, committeeSize) << ","
       << a.isSEJRPlus(resMES, committeeSize) << ","
       << a.isSEJRPlus(resGJCR, committeeSize) << std::endl;
    SEJRplus_t += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
}

void processMatricesFromFile(size_t committeeSize,
                            const std::string& infilename,
                            const std::string& outfilename){
    std::ifstream in(infilename);
    std::ofstream os(outfilename);

    std::cout << "PROCESS" << std::endl;

    // --- CSV HEADER ---
    os << "maxTreeSize,maxCycleSize,maxDistanceToCycle,"
       << "sinkCount,sourceCount,isolatedCount,"
       << "hammingDist_AV_to_SP_AV,"
       << "hammingDist_AV_to_CC_AV,"
       << "hammingDist_AV_to_PAV_AV,"
       << "hammingDist_AV_to_MES_AV,"
       << "hammingDist_AV_to_GJCR_AV,"
       << "JR_AV,JR_SP,JR_CC,JR_PAV,JR_MES,JR_GJCR,"
       << "EJR_AV,EJR_SP,EJR_CC,EJR_PAV,EJR_MES,EJR_GJCR,"
       << "EJR+_AV,EJR+_SP,EJR+_CC,EJR+_PAV,EJR+_MES,EJR+_GJCR,"
       << "PR_AV,PR_SP,PR_CC,PR_PAV,PR_MES,PR_GJCR,"
       << "CS_AV,CS_SP,CS_CC,CS_PAV,CS_MES,CS_GJCR,"
       << "SJR_AV,SJR_SP,SJR_CC,SJR_PAV,SJR_MES,SJR_GJCR,"
       << "SEJR+_AV,SEJR+_SP,SEJR+_CC,SEJR+_PAV,SEJR+_MES,SEJR+_GJCR\n";

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
        PAV[i] = 1.0 / i;
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
    double JR_time = 0, EJR_time = 0, EJRplus_time = 0;
    double PR_time = 0, CS_time = 0, SJR_time = 0, SEJRplus_time = 0;

    // --- THREADS ---
    size_t numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;

    std::mutex coutMutex;

    // =========================
    // WORKERS
    // =========================
    for (size_t t = 0; t < numThreads; t++) {
        workers.emplace_back([&]() {
            double lJR=0, lEJR=0, lEJRp=0, lPR=0, lCS=0, lSJR=0, lSEJRp=0;

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

                workWithMatrix(m, committeeSize, ss, CC, PAV,
                               lJR, lEJR, lEJRp, lPR, lCS, lSJR, lSEJRp);

                outputs[task.first] = ss.str();
            }

            std::lock_guard lock(timeMutex);
            JR_time += lJR;
            EJR_time += lEJR;
            EJRplus_time += lEJRp;
            PR_time += lPR;
            CS_time += lCS;
            SJR_time += lSJR;
            SEJRplus_time += lSEJRp;
        });
    }

    // =========================
    // READER THREAD
    // =========================
    std::thread reader([&]() {
        std::string line;

        for (size_t i = 0; i < count; i++) {

            printProgressBar(i+1, count);

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

    auto totalTime = (JR_time + EJR_time + EJRplus_time + PR_time + CS_time + SJR_time + SEJRplus_time);

    std::cout << "\n\nTotal axiom time - " << totalTime/1000 << std::endl;
    totalTime /= 100;
    std::cout <<"JR_time - " << JR_time/totalTime << " % EJR_time - " 
    << EJR_time/totalTime << " % EJRplus_time - " << EJRplus_time/totalTime << " % PR_time - " 
    << PR_time/totalTime << " % CS_time - " << CS_time/totalTime << " % SJR_time - " 
    << SJR_time/totalTime << " % SEJRplus_time - " << SEJRplus_time/totalTime << " % " << std::endl;
}

bool endsWith(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

int main(int argc, char* argv[]) {

    bool generate = false;
    bool process = false;

    std::string matrixFile = "";
    std::string resultsFile = "";

    size_t committeeSize = 0;
    size_t numberOfMatrixes = 0;
    size_t matrixSize = 0;

    for (int i = 1; i < argc; i++) {

        std::string arg = argv[i];

        if (arg == "-g" || arg == "--generate") {
            generate = true;
        }
        else if (arg == "-p" || arg == "--process") {
            process = true;
        }
        else if (endsWith(arg, ".txt")) {
            matrixFile = arg;
        }
        else if (endsWith(arg, ".csv")) {
            resultsFile = arg;
        }
        else if (arg == "-s") {
            matrixSize = stoul(argv[++i]);
        }        
        else if (arg == "-n") {
            numberOfMatrixes = stoul(argv[++i]);
        }        
        else if (arg == "-c") {
            committeeSize = stoul(argv[++i]);
        }
        else {
            std::cerr << "Unknown argument: " << arg << std::endl;
        }
    }

    if(matrixFile == "") {
        std::cout << "no file to read matrixes from or generate matrixes to (needs to be .txt file)\n";
        return 1;
    }
    if(generate) {
        if(numberOfMatrixes == 0) {
            std::cout << "number of matrixes to generate not set [-n]\n";
            return 1;
        }
        if(matrixSize == 0) {
            std::cout << "size of matrixes to generate not set [-s]\n";
            return 1;
        }
        std::random_device rd;
        generateRandomDelegationMatrixes(numberOfMatrixes, matrixSize, rd, matrixFile);
    }
    if(process) {
        if(resultsFile == "") {
            std::cout << "no file to save results was given (needs to be .csv file)\n";
            return 1;
        }
        if(committeeSize == 0) {
            std::cout << "elected committee not set [-c]\n";
            return 1;
        }
        processMatricesFromFile(committeeSize, matrixFile, resultsFile);
    }
    
    return 0;
}
