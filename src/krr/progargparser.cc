/*
 * progargparser.cpp
 *
 *  Created on: 31.10.2016
 *      Author: klaus
 */

#include "progargparser.h"

ProgramArgumentParser::ProgramArgumentParser(int &argc, char **argv) {
    for (int i = 1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

const std::string& ProgramArgumentParser::getCmdOption(
        const std::string &option) const {
    auto itr = std::find(tokens.begin(), tokens.end(), option);
    if (itr != tokens.end() && ++itr != tokens.end())
        return *itr;
    return *(new std::string(""));
}

const std::string &ProgramArgumentParser::getCmdArgument(
        const int reversePosition) const {
    if (reversePosition >= tokens.size()) {
        return *(new std::string(""));
    } else {
        const std::string option = tokens[tokens.size() - reversePosition - 1];
        if (option.empty() || option[0] == '-')
            return *(new std::string(""));
        return tokens[tokens.size() - reversePosition - 1];
    }
}

const std::string &ProgramArgumentParser::getStringValue(
        const std::string &flag) const {
    for (const std::string token : tokens) {
        size_t pos = token.find(flag);
        if (pos != std::string::npos) {
            pos = token.find_first_of("=");
            return *(new std::string(token.substr(pos + 1)));
        }
    }
    return NULL;
}

const std::vector<int> &ProgramArgumentParser::getIntegerArrayValue(
        const std::string &flag) const {
    std::vector<int> *res = new std::vector<int>;
    std::string valStr = getStringValue(flag);
    std::stringstream trimmer;
    trimmer << valStr;
    valStr.clear();
    trimmer >> valStr;
    if (!valStr.empty()) {
        if (valStr[0] != '[' && valStr[valStr.length() - 1] != ']') {
            std::cerr
                    << "Integeger array hast to start with '[' and end with ']!"
                    << std::endl;
            return *res;
        }
        char sep = ',';
        for (size_t p = 0, q = 0; p != valStr.npos; p = q) {
            std::string element = valStr.substr(p + (p != 0),
                    (q = valStr.find(sep, p + 1)) - p - (p != 0));
            if (p == 0)
                element = element.substr(1, element.length() - 1);
            if (element[element.length() - 1] == ']')
                element = element.substr(0, element.length() - 1);
            res->push_back(std::atoi(element.c_str()));
        }
    }
    return *res;
}

const int ProgramArgumentParser::getIntegerValue(
        const std::string &flag) const {
    const std::string valStr = getStringValue(flag);
    if (!valStr.empty()) {
        return std::atoi(valStr.c_str());
    }
    return 0;
}

const double ProgramArgumentParser::getDoubleValue(
        const std::string &flag) const {
    const std::string valStr = getStringValue(flag);
    if (!valStr.empty()) {
        return std::atof(valStr.c_str());
    }
    return 0.0;
}

bool ProgramArgumentParser::cmdOptionExists(const std::string &option) const {
    for (const std::string token : tokens) {
        if (token.find(option) != std::string::npos)
            return true;
    }
    return false;
}

