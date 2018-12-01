/*
 * utiles.cc
 *
 *  Created on: 21.09.2016
 *      Author: klaus
 */

#include <iostream>

const size_t getArrayPosRowMajor(const size_t *pos, const size_t *dims,
        const size_t rank) {
    size_t arrayPos = 0;
    for (size_t i = 1; i <= rank; i++) {
        size_t offset = 1;
        for (size_t j = i + 1; j <= rank; j++)
            offset *= dims[j - 1];
        arrayPos += offset * pos[i - 1];
    }
    return arrayPos;
}

const size_t getArrayPosColMajor(const size_t *pos, const size_t *dims,
        const size_t rank) {
    size_t arrayPos = 0;
    for (size_t i = 1; i <= rank; i++) {
        size_t offset = 1;
        for (size_t j = 1; j <= i - 1; j++)
            offset *= dims[j - 1];
        arrayPos += offset * pos[i - 1];
    }
    return arrayPos;
}

