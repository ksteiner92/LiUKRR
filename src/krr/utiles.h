/*
 * utiles.h
 *
 *  Created on: 21.09.2016
 *      Author: Klaus Steiner
 */

#ifndef UTILES_H_
#define UTILES_H_

const size_t getArrayPosColMajor(const size_t *pos, const size_t *dims,
        const size_t rank);

const size_t getArrayPosRowMajor(const size_t *pos, const size_t *dims,
        const size_t rank);

#endif /* UTILES_H_ */
