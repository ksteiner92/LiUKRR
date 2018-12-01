/*
 * dataset.h
 *
 *  Created on: 19.09.2016
 *      Author: Klaus Steiner
 */

#ifndef DATASET_H_
#define DATASET_H_

#include <bzlib.h>
#include <iostream>
#include <vector>
#include <matio.h>
#include <memory>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstring>
#include <random>
#include <fstream>
#include <map>

#include "machine.h"
#include "utiles.h"
#include "types.h"

template<class T> class Machine;

template<class T> class DataSet {
public:

    /**
     * Constructor of the DataSet class
     *
     * @param fileName The file or directory path where the
     * data are stored.
     */
    DataSet(const std::string fileName);

    /**
     * Virtual desctructor of the DataSet class.
     */
    virtual ~DataSet();

    /**
     * Creates an aligned sub set through a start and a end index.
     *
     * @param start The start index
     * @param end The end index
     */
    void setSubSet(const size_t start, const size_t end);

    /**
     * Creates sub set with random elements and size @size
     *
     * @param size The number of elements the sub set should have.
     */
    void setRandomSubSet(const size_t size);

    /**
     * Creates a sub set by using all elements
     */
    void setSubSetToAllExcludedElements();

    const std::string &getFileName() const;

    const size_t getMaxNumAtoms() const;

    const size_t getTotalNumCompounds() const;

    const size_t getNumCompoundsWithRestriction();

    void setRestrictOnNumAtoms(const size_t numAtoms);

    const size_t getNumAtomsWithRestriction() const;

    void getStructures(Structures<T> &structures) const;

    void getEnergies(T* &energies) const;

protected:
    std::string fileName;
    SubSet *subSet;
    bool *includeSetElements;
    Structures<T> restStructures;
    Structures<T> structures;

    inline void freeSubSet();

    inline void initIncludeSetElements();

};

class QMDataSet_7: public DataSet<float> {
public:
    QMDataSet_7(const std::string &fileName);

    void parseCoulombMatrices(std::vector<DescriptorMatrix<float>> &matrices);

protected:

    void parseStructures();

private:

    MatFileHandle *openMatFile(const mat_acc &acc = MAT_ACC_RDONLY) const;

};

class QMDataSet_7b: public QMDataSet_7 {
};

class VASPDataSet: public DataSet<double> {

    enum DirType {
        FINISHED, WAIT, RUN, SUBDIR
    };

    typedef std::function<bool(std::istream &buffer)> ParserType;

public:
    VASPDataSet(const std::string &workingDir, const bool felixVASPFileHeader =
            false);

    void preProcess(const std::string &vaspDir);

protected:

    void parseStructures();

private:
    bool compressed;
    std::string fileEnding;
    std::vector<std::string> sortedFileNames;
    bool felixVASPFileHeader;

    const bool preProcessDir(const char *dir,
            const VASPDataSet::DirType currentDirType);

    inline DIR *cdAndOpenDir(const char *dirName) const;

    const bool init();

    const bool parseDirectoryTree(const char* dirname,
            const std::string &relativePath,
            std::vector<std::string> &sortedFileNames);

    const bool parseVASPFile(const char *fileName, ParserType parser,
            const bool compressed = true) const;

};

#endif /* DATASET_H_ */
