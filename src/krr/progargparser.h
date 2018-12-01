/*
 * progargparser.h
 *
 *  Created on: 31.10.2016
 *      Author: klaus
 */

#ifndef PROGARGPARSER_H_
#define PROGARGPARSER_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>

class ProgramArgumentParser {
public:

    ProgramArgumentParser(int &argc, char **argv);

    const std::string& getCmdOption(const std::string &option) const;

    const std::string &getCmdArgument(const int reversePosition) const;

    const std::string &getStringValue(const std::string &flag) const;

    const std::vector<int> &getIntegerArrayValue(const std::string &flag) const;

    const double getDoubleValue(const std::string &flag) const;

    const int getIntegerValue(const std::string &flag) const;

    bool cmdOptionExists(const std::string &option) const;

private:
    std::vector<std::string> tokens;
};

#endif /* PROGARGPARSER_H_ */
