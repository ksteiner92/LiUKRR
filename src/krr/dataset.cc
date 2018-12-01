/*
 * dataset.cc
 *
 *  Created on: 19.09.2016
 *      Author: Klaus Steiner
 */

#include "dataset.h"

template<class T> DataSet<T>::DataSet(const std::string fileName) :
        fileName(fileName), subSet(NULL), includeSetElements(NULL) {
}

template<class T> DataSet<T>::~DataSet() {
    freeSubSet();
    if (includeSetElements != NULL)
        delete[] includeSetElements;
}

template<class T> const size_t DataSet<T>::getMaxNumAtoms() const {
    return structures.numAtoms;
}

template<class T> const size_t DataSet<T>::getTotalNumCompounds() const {
    return structures.structures.size();
}

template<class T> inline void DataSet<T>::freeSubSet() {
    if (subSet != NULL) {
        if (subSet->elements != NULL)
            delete[] subSet->elements;
        delete subSet;
        subSet = NULL;
    }
}

template<class T> inline void DataSet<T>::initIncludeSetElements() {
    if (includeSetElements != NULL)
        delete[] includeSetElements;
    const size_t numCompounds = getNumCompoundsWithRestriction();
    includeSetElements = new bool[numCompounds];
    std::fill(includeSetElements, includeSetElements + numCompounds, false);
}

template<class T> void DataSet<T>::setSubSet(const size_t start,
        const size_t end) {
    if (start >= end || end <= start || end > getNumCompoundsWithRestriction())
        throw std::out_of_range(
                "Sub set size bigger than total set size (include the number of atom restriction)!");
    freeSubSet();
    initIncludeSetElements();
    const size_t setSize = end - start;
    subSet = new SubSet;
    subSet->size = setSize;
    subSet->aligned = true;
    subSet->elements = new size_t[subSet->size];
    size_t idx = 0;
    for (size_t i = start; i < end; i++) {
        subSet->elements[i - start] = i;
        includeSetElements[i] = true;
    }
}

template<class T> void DataSet<T>::setRandomSubSet(const size_t size) {
    const size_t numCompounds = getNumCompoundsWithRestriction();
    if (size > numCompounds)
        throw std::out_of_range(
                "Sub set size bigger than total set size (include the number of atom restriction)!");
    freeSubSet();
    initIncludeSetElements();
    std::uniform_int_distribution<size_t> dist(0, numCompounds - 1);
    std::random_device rd; // uses RDRND or /dev/urandom
    subSet = new SubSet;
    subSet->size = size;
    subSet->aligned = false;
    subSet->elements = new size_t[size];
    std::fill_n(subSet->elements, size, 0);
    for (size_t i = 0; i < size; i++) {
        while (true) {
            const size_t idx = dist(rd);
            if (!includeSetElements[idx]) {
                subSet->elements[i] = idx;
                includeSetElements[idx] = true;
                break;
            }
        }
    }
}

template<class T> void DataSet<T>::setSubSetToAllExcludedElements() {
    std::vector<size_t> excludedElements;
    for (size_t i = 0; i < getNumCompoundsWithRestriction(); i++)
        if (!includeSetElements[i])
            excludedElements.push_back(i);
    freeSubSet();
    subSet = new SubSet;
    subSet->size = excludedElements.size();
    subSet->aligned = false;
    subSet->elements = new size_t[subSet->size];
    for (size_t i = 0; i < excludedElements.size(); i++)
        subSet->elements[i] = excludedElements[i];
}

template<class T> const std::string &DataSet<T>::getFileName() const {
    return fileName;
}

template<class T> void DataSet<T>::setRestrictOnNumAtoms(
        const size_t numAtoms) {
    restStructures.structures.clear();
    restStructures.numAtoms = numAtoms;
    if (numAtoms == 0)
        return;
    for (const Structure<T> structure : structures.structures)
        if (structure.atoms.size() <= numAtoms)
            restStructures.structures.push_back(structure);
    std::cout << "[Info] Restricted data set on compounds with " << numAtoms
            << " atoms." << std::endl;
    std::cout << "       Total number are now: "
            << restStructures.structures.size() << std::endl;
}

template<class T> const size_t DataSet<T>::getNumAtomsWithRestriction() const {
    return restStructures.structures.empty() ?
            structures.numAtoms : restStructures.numAtoms;
}

template<class T> const size_t DataSet<T>::getNumCompoundsWithRestriction() {
    return restStructures.structures.empty() ?
            structures.structures.size() : restStructures.structures.size();
}

template<class T> void DataSet<T>::getStructures(
        Structures<T> &structures) const {
    structures.numAtoms = getNumAtomsWithRestriction();
    if (subSet == NULL) {
        for (const Structure<T> structure : (
                restStructures.structures.empty() ?
                        this->structures.structures : restStructures.structures)) {
            structures.structures.push_back(structure);
        }
    } else {
        if (restStructures.structures.empty())
            for (size_t i = 0; i < subSet->size; i++)
                structures.structures.push_back(
                        this->structures.structures[subSet->elements[i]]);
        else
            for (size_t i = 0; i < subSet->size; i++)
                structures.structures.push_back(
                        restStructures.structures[subSet->elements[i]]);
    }
}

template<class T> void DataSet<T>::getEnergies(T* &energies) const {
    if (subSet == NULL) {
        size_t i = 0;
        for (const Structure<T> structure : (
                restStructures.structures.empty() ?
                        this->structures.structures : restStructures.structures)) {
            energies[i] = structure.energy;
            i++;
        }
    } else {
        if (restStructures.structures.empty())
            for (size_t i = 0; i < subSet->size; i++)
                energies[i] =
                        this->structures.structures[subSet->elements[i]].energy;
        else
            for (size_t i = 0; i < subSet->size; i++)
                energies[i] =
                        restStructures.structures[subSet->elements[i]].energy;
    }
}

VASPDataSet::VASPDataSet(const std::string &workingDir,
        const bool felixVASPFileHeader) :
        DataSet<double>(workingDir), compressed(true), fileEnding("bz2"), felixVASPFileHeader(
                felixVASPFileHeader) {
    parseStructures();
}

const bool VASPDataSet::parseDirectoryTree(const char* dirname,
        const std::string &relativePath,
        std::vector<std::string> &sortedFileNames) {
    DIR *dir = cdAndOpenDir(dirname);
    if (dir == NULL)
        return false;
    dirent *dent;
    struct stat st;
    while ((dent = readdir(dir))) {
        lstat(dent->d_name, &st);
        std::string newRelativePath(relativePath);
        newRelativePath.append("/").append(dent->d_name);
        if (!S_ISDIR(st.st_mode) && !std::string(dent->d_name).empty()) {
            std::string name(dent->d_name);
            size_t idx = name.find_last_of(".");
            if (idx < (name.length() - 1) && name.substr(idx + 1) == fileEnding)
                sortedFileNames.push_back(newRelativePath);
        } else if (std::string(dent->d_name) != "."
                && std::string(dent->d_name) != "..") {
            std::stringstream newDirName;
            newDirName << dirname << "/" << dent->d_name;
            parseDirectoryTree(newDirName.str().c_str(), newRelativePath,
                    sortedFileNames);
        }
    }
    closedir(dir);
    return true;
}

void VASPDataSet::parseStructures() {
    structures.structures.clear();
    std::string dirStr(fileName);
    std::stringstream trimmer;
    trimmer << dirStr;
    dirStr.clear();
    trimmer >> dirStr;
    if (dirStr[dirStr.size() - 1] == '/')
        dirStr.pop_back();
    fileName = dirStr;
    std::cout << "[Info] Initializing VASP data set...";
    std::ifstream fs;
    std::stringstream ss;
    ss << dirStr << "/info";
    fs.open(ss.str().c_str(), std::ios::in);
    if (fs) {
        std::string line;
        if (!std::getline(fs, line))
            throw "Info file of VASP data set contains no entries!";
        //numCompounds = std::atoi(line.c_str());
        if (!std::getline(fs, line))
            throw "Info file of VASP data set contains no information about number of atoms!";
        structures.numAtoms = std::atoi(line.c_str());
        if (!std::getline(fs, line))
            throw "Info file of VASP data set contains no information about file ending!";
        fileEnding = line;
        if (std::getline(fs, line)) {
            if (line == "compressed")
                compressed = true;
            else if (line == "uncompressed")
                compressed = false;
            else
                throw "Info file of VASP data set has an invalid info about compressen (compressed or uncompressed)!";
        }
        fs.close();
    } else {
        throw "Could not open info file!";
    }
    std::string relativePath("");
    if (!parseDirectoryTree(fileName.c_str(), relativePath, sortedFileNames))
        throw "Could not parse data directory!";
    std::sort(sortedFileNames.begin(), sortedFileNames.end());
    int count = 0;
    double averageEnergy = 0.0;
    for (size_t idx = 0; idx < sortedFileNames.size(); idx++) {
        const std::string name = sortedFileNames[idx];
        std::stringstream path;
        path << fileName << name;
        if (parseVASPFile(path.str().c_str(),
                [count, name, this](std::istream &buffer) {
                    Structure<double> structure;
                    std::string line;
                    if (std::getline(buffer, line)) {
                        std::string E0_str;
                        if(!felixVASPFileHeader) {
                            const size_t n = line.find_last_of(" ");
                            E0_str = line.substr(n + 1);
                            const size_t nameEnd = line.find_first_of(" ");
                            structure.name = line.substr(0, nameEnd);
                        } else {
                            const size_t start = line.find("energy");
                            const size_t end = line.find("symetry");
                            const size_t endSymm = line.find(" ", end + 9);
                            E0_str = line.substr(start + 7, end - 8 - start);
                            structure.name = line.substr(end + 9, endSymm - end - 9);
                        }
                        structure.energy = std::atof(E0_str.c_str());
                    } else {
                        std::cerr << "Structure file " << name << " contains no energy!"
                        << std::endl;
                        return false;
                    }
                    double scale = 1.0;
                    if (std::getline(buffer, line)) {
                        scale = std::atof(line.c_str());
                    } else {
                        std::cerr << "Structure file " << name << " contains no scaling factor!"
                        << std::endl;
                        return false;
                    }
                    std::fill_n(structure.base, 9, 0.0);
                    for (int i = 0; i < 3; i++) {
                        if (std::getline(buffer, line)) {
                            std::stringstream ss(line);
                            std::string coordinate;
                            int j = 0;
                            while (ss >> coordinate) {
                                structure.base[i * 3 + j] = std::atof(
                                        coordinate.c_str()) * scale;
                                j++;
                            }
                        } else {
                            std::cerr << "Structure file " << name
                            << " contains not enough basis vectors!"
                            << std::endl;
                            return false;
                        }
                    }
                    std::vector<Atom<double>> atoms;
                    if (std::getline(buffer, line)) {
                        std::stringstream ss(line);
                        std::string symbole;
                        while (ss >> symbole) {
                            const auto it = PeriodicTable.find(symbole);
                            if (it != PeriodicTable.end()) {
                                Atom<double> atom;
                                atom.Z = it->second;
                                std::fill_n(atom.r, 3, 0.0);
                                atoms.push_back(atom);
                            } else {
                                std::cerr << "Element Symbol '" << symbole
                                << "' not recognized!" << std::endl;
                                return false;
                            }
                        }
                    } else {
                        std::cerr << "Structure file " << name
                        << " contains no information about atom types!"
                        << std::endl;
                        return false;
                    }
                    if (std::getline(buffer, line)) {
                        std::stringstream ss(line);
                        std::string numStr;
                        int i = 0;
                        while (ss >> numStr) {
                            int num = std::atoi(numStr.c_str());
                            for (int j = 0; j < num; j++)
                            structure.atoms.push_back(atoms[i]);
                            i++;
                        }
                    } else {
                        std::cerr << "Structure file " << name
                        << " contains no information about number of atoms!"
                        << std::endl;
                        return false;
                    }
                    std::string latticeType;
                    bool directLattice = true;
                    if (!std::getline(buffer, latticeType)) {
                        std::cerr << "Structure file " << name
                        << " contains no information of latice type!"
                        << std::endl;
                        return false;
                    }
                    std::transform(latticeType.begin(),
                            latticeType.end(),
                            latticeType.begin(),
                            ::tolower);
                    if (latticeType != "direct" && latticeType != "d") {
                        directLattice = false;
                        std::cerr << "Structure file " << name
                        << " is not a direct lattice! (" << latticeType << ")"
                        << std::endl;
                    }
                    int i = 0;
                    while (std::getline(buffer, line)) {
                        std::stringstream ss(line);
                        std::string coordinate;
                        int j = 0;
                        double a[3] = {0.0, 0.0, 0.0};
                        while ((ss >> coordinate) && (j < 3)) {
                            a[j] = std::atof(coordinate.c_str());
                            if (directLattice)
                            a[j] = a[j] + (a[j] < 0 ? +1 : -1) * (((int) a[j]) + (a[j] < 0 ? 1 : 0));
                            j++;
                        }
                        if (directLattice) {
                            for (j = 0; j < 3; j++) {
                                for (int k = 0; k < 3; k++) {
                                    structure.atoms[i].r[j] += a[k] * structure.base[k * 3 + j];
                                }
                            }
                        } else {
                            for (j = 0; j < 3; j++)
                            structure.atoms[i].r[j] = scale * a[j];

                        }
                        i++;
                    }
                    structures.structures.push_back(structure);
                    return true;
                }, compressed)) {
            averageEnergy += structures.structures[count].energy;
            count++;
        }
    }
    averageEnergy /= structures.structures.size();
    for (Structure<double> structure : structures.structures)
        structure.energy -= averageEnergy;
    std::cout << "[OK]" << std::endl;
    std::cout.flush();
}

void VASPDataSet::preProcess(const std::string &vaspDir) {
    std::string dir(vaspDir);
    std::stringstream trimmer;
    trimmer << dir;
    dir.clear();
    trimmer >> dir;
    if (dir[dir.size() - 1] == '/')
        dir.pop_back();
    if (!preProcessDir(dir.c_str(), SUBDIR)) {
        throw "Pre-process failed for directory " + dir;
    }
}

inline DIR *VASPDataSet::cdAndOpenDir(const char *dirName) const {
    DIR *dir = NULL;
    if (chdir(dirName)) {
        std::cerr << "No permission to access directory " << dirName
                << std::endl;
        return NULL;
    }
    dir = opendir(".");
    if (dir == NULL) {
        std::cerr << "Could not open directory " << dirName << std::endl;
        return NULL;
    }
    return dir;
}

const bool VASPDataSet::preProcessDir(const char *dirName,
        VASPDataSet::DirType currentDirType) {
    DIR *dir = cdAndOpenDir(dirName);
    if (dir == NULL)
        return false;
    dirent *dent;
    struct stat st;
    while ((dent = readdir(dir))) {
        lstat(dent->d_name, &st);
        bool isDir = S_ISDIR(st.st_mode);
        if (isDir && currentDirType != RUN && std::string(dent->d_name) != "."
                && std::string(dent->d_name) != "..") {
            std::string dirNameStr(dent->d_name);
            int lastDotPos = dirNameStr.find_last_of(".");
            std::string typeName = dirNameStr.substr(lastDotPos + 1);
            if (typeName == "finished")
                currentDirType = FINISHED;
            else if (typeName == "waitstart")
                currentDirType = WAIT;
            else if (dirNameStr.find(".run") != std::string::npos)
                currentDirType = RUN;
            else
                currentDirType = SUBDIR;
            if (currentDirType != WAIT && currentDirType != RUN) {
                std::stringstream newDirName;
                newDirName << dirName << "/" << dent->d_name;
                preProcessDir(newDirName.str().c_str(), currentDirType);
            } else if (currentDirType == RUN) {
                std::stringstream newDirName;
                newDirName << dirName << "/" << dent->d_name;
                double E0 = 0.0;
                if (parseVASPFile(
                        newDirName.str().append("/OSZICAR.relax-final.bz2").c_str(),
                        [&E0](std::istream &buffer) {
                            std::string line;
                            while (std::getline(buffer, line)) {
                                size_t pos = 0;
                                while ((pos = line.find("E0=", pos)) != std::string::npos) {
                                    if (pos + 4 < line.length()) {
                                        pos += 4;
                                        std::string sub = line.substr(pos);
                                        std::stringstream subStream(sub);
                                        std::string valueStr;
                                        if (std::getline(subStream, valueStr, ' ')) {
                                            const double E0_new = std::atof(valueStr.c_str());
                                            if (E0_new < E0)
                                            E0 = E0_new;
                                        }
                                    } else
                                    break;
                                }
                            }
                            return true;
                        })) {
                    std::cout << "E0=" << E0 << std::endl;
                    //struct stat fileStat;
                    /*char result[ PATH_MAX];
                     ssize_t count = readlink("/proc/self/exe", result,
                     PATH_MAX);
                     std::string progPath(result, (count > 0) ? count : 0);
                     progPath = progPath.substr(0, progPath.find_last_of("/"));
                     std::string cmd(
                     "bzcat "
                     + newDirName.str().append(
                     "/XDATCAR.relax-final.bz2") + " > "
                     + progPath + "/.XDATCAR && " + progPath
                     + "/run_isotropy " + progPath
                     + "/.XDATCAR > " + progPath
                     + "/.XDATCAR_sym");*/
                    /*std::string cmd(
                     "bzcat "
                     + newDirName.str().append(
                     "/XDATCAR.relax-final.bz2")
                     + " > .XDATCAR && ./run_isotropy > .XDATCAR_sym");*/
                    if (!std::system(
                            ("bzcat "
                                    + newDirName.str().append(
                                            "/XDATCAR.relax-final.bz2")
                                    + " > .XDATCAR").c_str())) {
                    }
                    //FILE* in = popen(command.c_str(), "r");
                    /*std::cout << cmd << std::endl;
                     if (!std::system(cmd.c_str())) {
                     std::cerr << "Run Isotropy faild on file '" << fileName
                     << "'" << std::endl;
                     return false;
                     }*/
                    if (parseVASPFile(".XDATCAR_sym",
                            [E0](std::istream &buffer) {
                                std::string line;
                                while (std::getline(buffer, line)) {
                                    std::cout << line << std::endl;
                                }
                                return true;
                            })) {
                        return true;
                    } else {
                        std::cerr << "Could not parse 'XDATCAR.relax-final.bz2'"
                                << std::endl;
                        return false;
                    }
                } else {
                    std::cerr << "Could not parse 'OSZICAR.relax-final.bz2'"
                            << std::endl;
                    return false;
                }
            }
        }
    }
    closedir(dir);
    return true;
}

const bool VASPDataSet::parseVASPFile(const char *fileName, ParserType parser,
        const bool compressed) const {
    if (compressed) {
        FILE *vasp = fopen(fileName, "r");
        if (vasp == NULL) {
            std::cerr << "Could not open VASP file '" << fileName << "'"
                    << std::endl;
            return false;
        }
        int bzerror;
        BZFILE* vaspBZ = BZ2_bzReadOpen(&bzerror, vasp, 0, 0, NULL, 0);
        if (bzerror != BZ_OK) {
            BZ2_bzReadClose(&bzerror, vaspBZ);
            fclose(vasp);
            std::cerr << "Could not open compressed VASP file '" << fileName
                    << "'" << std::endl;
            return false;
        }
        std::stringstream buffer;
        const int chunkSize = 2048;
        char chunk[chunkSize];
        int readedBytes = 0;
        while (bzerror == BZ_OK) {
            readedBytes = BZ2_bzRead(&bzerror, vaspBZ, chunk, chunkSize);
            for (int i = 0; i < readedBytes; i++)
                buffer << chunk[i];
        }
        if (!parser(buffer)) {
            BZ2_bzReadClose(&bzerror, vaspBZ);
            fclose(vasp);
            std::cerr << "Could not parse compressed VASP file '" << fileName
                    << "'" << std::endl;
            return false;
        }
        if (bzerror != BZ_STREAM_END) {
            BZ2_bzReadClose(&bzerror, vaspBZ);
            fclose(vasp);
            std::cerr << "Could not read compressed VASP file '" << fileName
                    << "'" << std::endl;
            return false;
        }
        BZ2_bzReadClose(&bzerror, vaspBZ);
        fclose(vasp);
    } else {
        std::ifstream buffer;
        buffer.open(fileName, std::ios::in);
        if (buffer) {
            if (!parser(buffer)) {
                std::cerr << "Could not parse VASP file '" << fileName << "'"
                        << std::endl;
                buffer.close();
                return false;
            }
            buffer.close();
        }
    }
    return true;
}

QMDataSet_7::QMDataSet_7(const std::string &fileName) :
        DataSet<float>(fileName) {
    parseStructures();
}

MatFileHandle *QMDataSet_7::openMatFile(const mat_acc &acc) const {
    MatFileHandle *matfp = NULL;
    matfp = Mat_Open(fileName.c_str(), acc);
    if (matfp == NULL) {
        throw "Could not open MAT file";
    }
    return matfp;
}

void QMDataSet_7::parseStructures() {
    MatFileHandle* matfp = openMatFile();
    MatVariable *matvarNuc = Mat_VarRead(matfp, "Z");
    if (matvarNuc == NULL) {
        throw "Variable 'Z' not found, or error "
                "reading MAT file!";
    }
    float *nuclearCharges = (float *) matvarNuc->data;
    MatVariable* matvarR = Mat_VarRead(matfp, "R");
    if (matvarR == NULL) {
        throw "Variable 'R' not found, or error "
                "reading MAT file!";
    }
    float* r = (float *) matvarR->data;
    MatVariable* matvarE = Mat_VarRead(matfp, "T");
    if (matvarE == NULL) {
        throw "Variable 'T' not found, or error "
                "reading MAT file!";
    }
    std::cout << "[Info] Copying formation engeries... ";
    float *energies = (float *) matvarE->data;
    size_t pos[3];
    size_t posZ[2];
    structures.numAtoms = matvarR->dims[1];
    for (size_t i = 0; i < matvarR->dims[0]; i++) {
        pos[0] = i;
        posZ[0] = pos[0];
        Structure<float> structure;
        structure.energy = energies[i];
        for (pos[1] = 0; pos[1] < matvarR->dims[1]; pos[1]++) {
            posZ[1] = pos[1];
            int Z =
                    nuclearCharges[getArrayPosColMajor(posZ, matvarNuc->dims, 2)];
            if (Z == 0)
                break;
            Atom<float> atom;
            pos[2] = 0;
            atom.r[0] = r[getArrayPosColMajor(pos, matvarR->dims, 3)];
            pos[2] = 1;
            atom.r[1] = r[getArrayPosColMajor(pos, matvarR->dims, 3)];
            pos[2] = 2;
            atom.r[2] = r[getArrayPosColMajor(pos, matvarR->dims, 3)];
            atom.Z = Z;
            structure.atoms.push_back(atom);
        }
        structures.structures.push_back(structure);
    }
    Mat_VarFree(matvarR);
    Mat_VarFree(matvarNuc);
    Mat_VarFree(matvarE);
    Mat_Close(matfp);
}

void QMDataSet_7::parseCoulombMatrices(
        std::vector<DescriptorMatrix<float>> &matrices) {
    /* Open input mat file */
    MatFileHandle *matfp = openMatFile();
    MatVariable *matvar = Mat_VarRead(matfp, "X");
    if (matvar == NULL) {
        throw "Variable 'X' not found, or error reading MAT file!";
    }
    std::cout << "[Info] Checking coulomb matrices... ";
    if (matvar->rank != 3) {
        throw std::invalid_argument(
                "Given matlab variable is not a coulomb matrix or just one molecule!");
    }
    if (matvar->dims[1] != matvar->dims[2]) {
        throw std::invalid_argument(
                "Given coulomb matrix is not a square matrix!");
    }
    std::cout << "[OK]" << std::endl << "[Info] Loading Coulomb Matrices... ";
    const size_t numAtoms = matvar->dims[1];
    const size_t descSize = numAtoms * numAtoms;
    float *c = (float *) matvar->data;
    size_t pos[matvar->rank];
    for (size_t i = 0; i < (subSet != NULL ? subSet->size : matvar->dims[0]);
            i++) {
        pos[0] = subSet != NULL ? subSet->elements[i] : i;
        float* matrix = new float[descSize];
        std::fill_n(matrix, descSize, 0.0);
        DescriptorMatrix<float> descriptor;
        descriptor.numAtoms = -1;
        size_t offset = 0;
        for (pos[1] = 0; pos[1] < numAtoms; pos[1]++) {
            for (pos[2] = 0; pos[2] < numAtoms; pos[2]++) {
                const float val = c[getArrayPosColMajor(pos, matvar->dims,
                        matvar->rank)];
                if (descriptor.numAtoms < 0 && val == 0)
                    descriptor.numAtoms = pos[2];
                matrix[pos[1] * numAtoms + pos[2]] = val;
            }
        }
        descriptor.m = std::shared_ptr<float>(matrix);
        matrices.push_back(descriptor);
    }
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    std::cout << "[OK]" << std::endl;
    std::cout.flush();
}

template class DataSet<float> ;
template class DataSet<double> ;
