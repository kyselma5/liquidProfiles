#pragma once

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <iostream>

struct Config {
    bool generate = false;
    bool process = false;
    bool JR = false;
    bool PJR = false;
    bool EJR = false;
    bool EJRplus = false;
    bool LR = false;
    bool CS = false;
    bool PR = false;
    bool sJR = false;
    bool sEJR = false;
    bool sEJRplus = false;

    bool doSanityCheck = false;

    std::string matrixFile;
    std::string resultsFile;

    size_t matrixSize = 0;
    size_t numberOfMatrixes = 0;

    std::vector<size_t> committeeSizes;
};

void printConfig(const Config& cfg) {
    std::cout << "===== CONFIG =====\n";

    std::cout << std::boolalpha; // true/false místo 1/0

    std::cout << "generate: " << cfg.generate << "\n";
    std::cout << "process: " << cfg.process << "\n";

    std::cout << "doSanityCheck: " << cfg.doSanityCheck << "\n";

    std::cout << "\n--- Rules ---\n";
    std::cout << "JR: " << cfg.JR << "\n";
    std::cout << "PJR: " << cfg.PJR << "\n";
    std::cout << "EJR: " << cfg.EJR << "\n";
    std::cout << "EJRplus: " << cfg.EJRplus << "\n";
    std::cout << "LR: " << cfg.LR << "\n";
    std::cout << "CS: " << cfg.CS << "\n";
    std::cout << "PR: " << cfg.PR << "\n";
    std::cout << "sJR: " << cfg.sJR << "\n";
    std::cout << "sEJR: " << cfg.sEJR << "\n";
    std::cout << "sEJRplus: " << cfg.sEJRplus << "\n";

    std::cout << "\n--- Files ---\n";
    std::cout << "matrixFile: " << cfg.matrixFile << "\n";
    std::cout << "resultsFile: " << cfg.resultsFile << "\n";

    std::cout << "\n--- Parameters ---\n";
    std::cout << "matrixSize: " << cfg.matrixSize << "\n";
    std::cout << "numberOfMatrixes: " << cfg.numberOfMatrixes << "\n";

    std::cout << "committeeSizes: ";
    for (size_t x : cfg.committeeSizes) {
        std::cout << x << " ";
    }
    std::cout << "\n";

    std::cout << "==================\n";
}

std::string trim(std::string s) {
    s.erase(0, s.find_first_not_of(" \t\r\n"));
    s.erase(s.find_last_not_of(" \t\r\n") + 1);
    return s;
}

bool toBool(std::string s) {
    s = trim(s);
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s == "true" || s == "1";
}

std::vector<size_t> parseList(const std::string& s) {
    std::vector<size_t> result;
    std::stringstream ss(s);
    std::string item;

    while (ss >> item) {
        std::stringstream ss2(item);
        std::string sub;

        while (std::getline(ss2, sub, ',')) {
            sub = trim(sub);
            if (!sub.empty())
                result.push_back(std::stoul(sub));
        }
    }
    return result;
}

Config loadConfig(const std::string& filename) {
    std::ifstream file(filename);
    Config cfg;

    if (!file) {
        throw std::runtime_error("Cannot open config file");
    }

    std::string line;
    while (std::getline(file, line)) {
        line = trim(line);

        if (line.empty() || line[0] == '#') continue;

        auto pos = line.find('=');
        if (pos == std::string::npos) continue;

        std::string key = trim(line.substr(0, pos));
        std::string value = trim(line.substr(pos + 1));

        if (key == "generate") cfg.generate = toBool(value);
        else if (key == "process") cfg.process = toBool(value);

        else if (key == "JR") cfg.JR = toBool(value);
        else if (key == "PJR") cfg.PJR = toBool(value);
        else if (key == "EJR") cfg.EJR = toBool(value);
        else if (key == "EJRplus") cfg.EJRplus = toBool(value);
        else if (key == "LR") cfg.LR = toBool(value);
        else if (key == "CS") cfg.CS = toBool(value);
        else if (key == "PR") cfg.PR = toBool(value);
        else if (key == "sJR") cfg.sJR = toBool(value);
        else if (key == "sEJR") cfg.sEJR = toBool(value);
        else if (key == "sEJRplus") cfg.sEJRplus = toBool(value);

        else if (key == "doSanityCheck") cfg.doSanityCheck = toBool(value);

        else if (key == "matrixFile") cfg.matrixFile = value;
        else if (key == "resultsFile") cfg.resultsFile = value;
        else if (key == "matrixSize") cfg.matrixSize = std::stoul(value);
        else if (key == "numberOfMatrixes") cfg.numberOfMatrixes = std::stoul(value);
        else if (key == "committeeSizes") cfg.committeeSizes = parseList(value);
    }

    return cfg;
}