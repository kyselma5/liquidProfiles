#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include "bitmatrix.h"
#include "rules.h"
#include "axioms.h"

double JR_time = 0;
double EJR_time = 0;
double EJRplus_time = 0;
double PR_time = 0;
double CS_time = 0;
double SJR_time = 0;
double SEJRplus_time = 0;


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
            printProgressBar(i, count);
        }
        generateDelegatinMatrix(size, rd).print(out);
        out << "\n";
    }
    printProgressBar(count, count);
    std::cout << std::endl;
}

void workWithMatrix(BitMatrix & m, size_t committeeSize, std::ofstream & os, const std::vector<double> & CC, const std::vector<double> & PAV) {
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
    auto start = std::chrono::high_resolution_clock::now();
    os << a.isJR(resAV, committeeSize) << ","
       << a.isJR(resSP, committeeSize) << ","
       << a.isJR(resCC, committeeSize) << ","
       << a.isJR(resPAV, committeeSize) << ","
       << a.isJR(resMES, committeeSize) << ","
       << a.isJR(resGJCR, committeeSize) << ",";
    JR_time += (std::chrono::high_resolution_clock::now() - start).count();
    start = std::chrono::high_resolution_clock::now();
    os << a.isEJR(resAV, committeeSize) << ","
       << a.isEJR(resSP, committeeSize) << ","
       << a.isEJR(resCC, committeeSize) << ","
       << a.isEJR(resPAV, committeeSize) << ","
       << a.isEJR(resMES, committeeSize) << ","
       << a.isEJR(resGJCR, committeeSize) << ",";
    EJR_time += (std::chrono::high_resolution_clock::now() - start).count();
    start = std::chrono::high_resolution_clock::now();
    os << a.isEJRplus(resAV, committeeSize) << ","
       << a.isEJRplus(resSP, committeeSize) << ","
       << a.isEJRplus(resCC, committeeSize) << ","
       << a.isEJRplus(resPAV, committeeSize) << ","
       << a.isEJRplus(resMES, committeeSize) << ","
       << a.isEJRplus(resGJCR, committeeSize) << ",";
    EJRplus_time += (std::chrono::high_resolution_clock::now() - start).count();
    start = std::chrono::high_resolution_clock::now();
    os << a.isPR(resAV, committeeSize) << ","
       << a.isPR(resSP, committeeSize) << ","
       << a.isPR(resCC, committeeSize) << ","
       << a.isPR(resPAV, committeeSize) << ","
       << a.isPR(resMES, committeeSize) << ","
       << a.isPR(resGJCR, committeeSize) << ",";
    PR_time += (std::chrono::high_resolution_clock::now() - start).count();
    start = std::chrono::high_resolution_clock::now();
    os << a.isCS(resAV, committeeSize) << ","
       << a.isCS(resSP, committeeSize) << ","
       << a.isCS(resCC, committeeSize) << ","
       << a.isCS(resPAV, committeeSize) << ","
       << a.isCS(resMES, committeeSize) << ","
       << a.isCS(resGJCR, committeeSize) << ",";
    CS_time += (std::chrono::high_resolution_clock::now() - start).count();
    start = std::chrono::high_resolution_clock::now();
    os << a.isSJR(resAV, committeeSize) << ","
       << a.isSJR(resSP, committeeSize) << ","
       << a.isSJR(resCC, committeeSize) << ","
       << a.isSJR(resPAV, committeeSize) << ","
       << a.isSJR(resMES, committeeSize) << ","
       << a.isSJR(resGJCR, committeeSize) << ",";
    SJR_time += (std::chrono::high_resolution_clock::now() - start).count();
    start = std::chrono::high_resolution_clock::now();
    os << a.isSEJRPlus(resAV, committeeSize) << ","
       << a.isSEJRPlus(resSP, committeeSize) << ","
       << a.isSEJRPlus(resCC, committeeSize) << ","
       << a.isSEJRPlus(resPAV, committeeSize) << ","
       << a.isSEJRPlus(resMES, committeeSize) << ","
       << a.isSEJRPlus(resGJCR, committeeSize) << std::endl;
    SEJRplus_time += (std::chrono::high_resolution_clock::now() - start).count();
}

void processMatricesFromFile(size_t committeeSize, const std::string& infilename, const std::string& outfilename) {

    std::cout << "PROCESSING\n";
    std::ifstream in(infilename);
    std::ofstream os(outfilename);

    os << "maxTreeSize" << "," 
       << "maxCycleSize" << "," 
       << "maxDistanceToCycle" << "," 
       << "sinkCount" << "," 
       << "sourceCount" << "," 
       << "isolatedCount" << "," 
       << "hammingDist_AV_to_SP_AV" << ","
       << "hammingDist_AV_to_CC_AV" << ","
       << "hammingDist_AV_to_PAV_AV" << ","
       << "hammingDist_AV_to_MES_AV" << ","
       << "hammingDist_AV_to_GJCR_AV" << ","
       << "JR_AV" << ","
       << "JR_SP" << ","
       << "JR_CC" << ","
       << "JR_PAV" << ","
       << "JR_MES" << ","
       << "JR_GJCR" << ","
       << "EJR_AV" << ","
       << "EJR_SP" << ","
       << "EJR_CC" << ","
       << "EJR_PAV" << ","
       << "EJR_MES" << ","
       << "EJR_GJCR" << ","
       << "EJR+_AV" << ","
       << "EJR+_SP" << ","
       << "EJR+_CC" << ","
       << "EJR+_PAV" << ","
       << "EJR+_MES" << ","
       << "EJR+_GJCR" << ","
       << "PR_AV" << ","
       << "PR_SP" << ","
       << "PR_CC" << ","
       << "PR_PAV" << ","
       << "PR_MES" << ","
       << "PR_GJCR" << ","
       << "CS_AV" << ","
       << "CS_SP" << ","
       << "CS_CC" << ","
       << "CS_PAV" << ","
       << "CS_MES" << ","
       << "CS_GJCR" << ","
       << "SJR_AV" << ","
       << "SJR_SP" << ","
       << "SJR_CC" << ","
       << "SJR_PAV" << ","
       << "SJR_MES" << ","
       << "SJR_GJCR" << ","
       << "SEJR+_AV" << ","
       << "SEJR+_SP" << ","
       << "SEJR+_CC" << ","
       << "SEJR+_PAV" << ","
       << "SEJR+_MES" << ","
       << "SEJR+_GJCR" << std::endl;

    size_t count, size;
    in >> count >> size;
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<double> CC(size, 0);
    CC[0] = 1;
    std::vector<double> PAV(size, 0);
    PAV[0] = 1;
    for(size_t i = 1; i < size; i++){
        PAV[i] = 1.0/i;
    }

    size_t prevPercentage = 0;

    for (size_t i = 0; i < count; i++) {
        if(100*i/count != prevPercentage) {
            prevPercentage = 100*i/count;
            printProgressBar(i, count);
        }
        std::stringstream matrixStream;
        std::string line;

        for (size_t r = 0; r < size; r++) {
            std::getline(in, line);
            matrixStream << line << "\n";
        }
        std::string matrixText = matrixStream.str();

        BitMatrix m(matrixText);

        workWithMatrix(m, committeeSize, os, CC, PAV);

        std::getline(in, line);
    }
    printProgressBar(count, count);

    std::cout << std::endl;
    auto totalTime = (JR_time + EJR_time + EJRplus_time + PR_time + CS_time + SJR_time + SEJRplus_time)/100;

    std::cout << "JR_time - " << JR_time/totalTime << " % EJR_time - " << EJR_time/totalTime 
                << " % EJRplus_time - " << EJRplus_time/totalTime << " % PR_time - " << PR_time/totalTime 
                << " % CS_time - " << CS_time/totalTime << " % SJR_time - " << SJR_time/totalTime 
                << " % SEJRplus_time - " << SEJRplus_time/totalTime << " % " << std::endl;
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
