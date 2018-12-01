/*
 * unit.cc
 *
 *  Created on: 13.11.2016
 *      Author: Klaus Steiner
 */

#include <iostream>

#include "quantity.h"

int main(int argc, char** argv) {
    double val = 5.0;
    unit::Quantity<unit::eV> energy =  val * unit::milli * "eV"_u;
    std::cout << energy() << std::endl;
    return 0;
}

