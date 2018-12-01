/*
 * analyzer.cc
 *
 *  Created on: 06.12.2016
 *      Author: Klaus Steiner
 */

#include <iostream>
#include <chrono>
#include <fstream>

#include "machine.h"
#include "dataset.h"
#include "descriptor.h"
#include "progargparser.h"

#define CUT_OUTLIERS 0

struct RunConfiguration {
    std::string destDir;
    size_t *biosPoints;
    size_t numBiosPoints;
    size_t numRunsPerBiosPoint;
};

enum SetType {
    VASP, VASP_FELIX, QM7
};

void printHelp(const std::string &programName) {
    std::cout << "Usage: " << programName
            << " <Data Type> [Options] <Input File/Dir> <Output Dir>"
            << std::endl;
    std::cout << "Data Type:" << std::endl;
    std::cout << "  --data=QM7:\tError analyze for QM7 data set." << std::endl;
    std::cout << "  --data=VASP:\tError analyze for VASP data set."
            << std::endl;
    std::cout
            << "  --data=VASP_FELIX:\tError analyze for VASP data set with header from Felix."
            << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --RunsPerBiosPoint=<Integer>:\t" << std::endl;
    std::cout
            << "  --BiosPoints=<Array>:\ti.e. must be given in the format [1, 10, 100]"
            << std::endl;
    std::cout << "  --sorting=<Name>:\t" << std::endl;
    std::cout << "                   \tROW_L1_DESC: default" << std::endl;
    std::cout << "                   \tROW_L2_DESC: " << std::endl;
    std::cout
            << "  --maxAtoms=<Integer>:\tRestricts the maximum numbers per crystal."
            << std::endl;
    std::cout << "<Input File/Dir>:\t" << std::endl;
    std::cout << "<Output Dir>:\t" << std::endl;
}

template<class T> void analyze(Machine<T> &qm, DataSet<T> *data,
        const RunConfiguration &conf) {
    const auto begin = std::chrono::steady_clock::now();
    Structures<T> structures;
    data->getStructures(structures);
    std::ofstream fs_structures_info;
    fs_structures_info.open((conf.destDir + "/structures_info.dat").c_str(),
            std::ios::out);
    fs_structures_info << "# name" << std::setw(30) << "num atoms"
            << std::setw(30) << "real value" << std::endl;
    for (const Structure<T> structure : structures.structures)
        fs_structures_info << structure.name << std::setw(30)
                << structure.atoms.size() << std::setw(30) << structure.energy
                << std::endl;
    fs_structures_info.close();
    std::ofstream fs_mae;
    fs_mae.open((conf.destDir + "/mae.dat").c_str(), std::ios::out);
    fs_mae << "# size" << std::setw(30) << "MAE" << std::setw(30)
            << "MAE per atom" << std::endl;
    for (size_t i = 0; i < conf.numBiosPoints; i++) {
        const size_t n = conf.biosPoints[i];
        std::cout << "---- Running " << conf.numRunsPerBiosPoint
                << " random set(s) with size " << n << " ----" << std::endl;
        PredictionResults average;
        average.mae = 0.0;
        average.maePerAtom = 0.0;
        for (size_t j = 0; j < conf.numRunsPerBiosPoint; j++) {

            std::ofstream fs_error;
            std::stringstream ss;
            ss << conf.destDir << "/error_" << n << ".dat." << j;
            fs_error.open(ss.str().c_str(), std::ios::out);
            fs_error << "# name" << std::setw(30) << "num atoms"
                    << std::setw(30) << "prediction" << std::setw(30)
                    << "real value" << std::setw(30) << "error" << std::endl;
            data->setRandomSubSet(n);
            qm.train(data);
            data->setSubSetToAllExcludedElements();
            PredictionResults results;
            qm.predict(data, results);

#if CUT_OUTLIERS
            double av_error = 0.0;
            double av_error_per_atom = 0.0;
            for (int k = 0; k < results.results.size(); k++) {
                const double abs_error = fabs(results.results[i].error);
                const double abs_error_per_atom = fabs(
                        results.results[i].errorPerAtom);
                if (abs_error >= 10.0) {
                    results.results.erase(results.results.begin() + i);
                } else {
                    av_error += abs_error;
                    av_error_per_atom += abs_error_per_atom;
                }
            }
            average.mae = av_error / results.results.size();
            average.maePerAtom = av_error_per_atom / results.results.size();
#endif

            std::cout << "\tMAE: " << results.mae << std::endl;
            std::cout << "\tMAE per atom: " << results.maePerAtom << std::endl;
            const auto end = std::chrono::steady_clock::now();
            std::cout << "\tTime: "
                    << std::chrono::duration_cast<std::chrono::seconds>(
                            end - begin).count() << "s" << std::endl;
            average.mae += results.mae;
            average.maePerAtom += results.maePerAtom;
            for (const Prediction pred : results.results) {
                fs_error << pred.name << std::setw(30) << pred.numAtoms
                        << std::setw(30) << pred.value << std::setw(30)
                        << (pred.value + pred.error) << std::setw(30)
                        << pred.error << std::endl;
            }
            fs_error.close();
        }
        average.mae /= conf.numRunsPerBiosPoint;
        average.maePerAtom /= conf.numRunsPerBiosPoint;
        std::cout << std::endl << "-------" << std::endl << "\tAverage MAE: "
                << average.mae << std::endl;
        std::cout << "\tAverage MAE per atom: " << average.maePerAtom
                << std::endl;
        const auto end = std::chrono::steady_clock::now();
        std::cout << "\tTime: "
                << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()
                << "s" << std::endl << std::endl;
        fs_mae << n << std::setw(30) << average.mae << std::setw(30)
                << average.maePerAtom << std::endl;
    }
    fs_mae.close();
}

int main(int argc, char** argv) {
    std::cout << "Start analayzing..." << std::endl;
    SetType type;
    RunConfiguration conf;
    conf.numRunsPerBiosPoint = 3;
    ProgramArgumentParser programArgs(argc, argv);
    if (!programArgs.cmdOptionExists("--data=")) {
        std::cerr << "Could not find data set type flag!" << std::endl;
        printHelp(argv[0]);
        return 1;
    }
    std::string setType = programArgs.getStringValue("--data");
    if (setType.empty()) {
        std::cerr << "Data set type has no value!" << std::endl;
        printHelp(argv[0]);
        return 1;
    }
    if (setType == "VASP") {
        type = VASP;
    } else if (setType == "VASP_FELIX") {
        type = VASP_FELIX;
    } else if (setType == "QM7") {
        type = QM7;
    } else {
        std::cerr << "Data set type has invalid value!" << std::endl;
        printHelp(argv[0]);
        return 1;
    }
    if (programArgs.cmdOptionExists("--RunsPerBiosPoint")) {
        int runsPerBiosPoint = programArgs.getIntegerValue(
                "--RunsPerBiosPoint");
        if (runsPerBiosPoint > 0) {
            conf.numRunsPerBiosPoint = runsPerBiosPoint;
        } else {
            std::cerr << "--RunsPerBiosPoint has invalid value!" << std::endl;
            printHelp(argv[0]);
            return 1;
        }
    }
    size_t maxAtoms = 0;
    if (programArgs.cmdOptionExists("--maxAtoms")) {
        maxAtoms = programArgs.getIntegerValue("--maxAtoms");
        if (maxAtoms < 0) {
            std::cerr << "--maxAtoms has invalid value!" << std::endl;
            printHelp(argv[0]);
            return 1;
        }
    }
    const bool defaultConfForBiosPoints = !programArgs.cmdOptionExists(
            "--BiosPoints");
    if (!defaultConfForBiosPoints) {
        std::vector<int> biosPoints = programArgs.getIntegerArrayValue(
                "--BiosPoints");
        if (!biosPoints.empty()) {
            conf.numBiosPoints = biosPoints.size();
            conf.biosPoints = new size_t[conf.numBiosPoints];
            for (int i = 0; i < conf.numBiosPoints; i++)
                conf.biosPoints[i] = biosPoints[i];
        } else {
            std::cerr << "--BiosPoints has invalid value!" << std::endl;
            printHelp(argv[0]);
            return 1;
        }
    }
    const std::string inputFile = programArgs.getCmdArgument(1);
    if (inputFile.empty()) {
        std::cerr << "[ERROR] Must specify a input file or directory."
                << std::endl;
        printHelp(argv[0]);
        return 1;
    }
    if (inputFile.find("--data") != std::string::npos
            || inputFile.find("--de") != std::string::npos
            || inputFile.find("--BiosPoints") != std::string::npos
            || inputFile.find("--RunsPerBiosPoint") != std::string::npos
            || inputFile.find("--sorting") != std::string::npos) {
        std::cerr << "[ERROR] Must specify a input file or directory."
                << std::endl;
        printHelp(argv[0]);
        return 1;
    }
    conf.destDir = programArgs.getCmdArgument(0);
    if (conf.destDir.empty()) {
        std::cerr << "[ERROR] Must specify a output directory." << std::endl;
        printHelp(argv[0]);
        return 1;
    }
    if (conf.destDir.find("--data") != std::string::npos
            || conf.destDir.find("--de") != std::string::npos
            || conf.destDir.find("--BiosPoints") != std::string::npos
            || conf.destDir.find("--RunsPerBiosPoint") != std::string::npos
            || conf.destDir.find("--sorting") != std::string::npos) {
        std::cerr << "[ERROR] Must specify a output directory." << std::endl;
        printHelp(argv[0]);
        return 1;
    }
    try {
        if (type == QM7) {
            DataSet<float> *data = new QMDataSet_7(inputFile);
            data->setRestrictOnNumAtoms(maxAtoms);
            if (defaultConfForBiosPoints) {
                conf.numBiosPoints = 6;
                conf.biosPoints = new size_t[6];
                const size_t maxTrainSetSize =
                        data->getNumCompoundsWithRestriction() * 0.75;
                const double k = pow((double) maxTrainSetSize, 1.0 / 6.0);
                conf.biosPoints[0] = 1;
                conf.biosPoints[1] = pow(k, 2);
                conf.biosPoints[2] = pow(k, 3);
                conf.biosPoints[3] = pow(k, 4);
                conf.biosPoints[4] = pow(k, 5);
                conf.biosPoints[5] = maxTrainSetSize;
            }
            Descriptor<float> *descriptor = new CoulombDescriptor<float>;
            descriptor->setSorting(SortingAlgorithm<float>::ROW_L1_DESC);
            Machine<float> qm(descriptor, Kernel<float>::LAPLACE);
            qm.setSigma(2.5e3);
            qm.setLambda(1.0e-6);
            analyze<float>(qm, data, conf);
            delete data;
            delete descriptor;
        } else if (type == VASP || type == VASP_FELIX) {
            VASPDataSet *data = new VASPDataSet(inputFile,
                    type == VASP_FELIX ? true : false);
            data->setRestrictOnNumAtoms(maxAtoms);
            if (defaultConfForBiosPoints) {
                conf.numBiosPoints = 6;
                conf.biosPoints = new size_t[6];
                const size_t maxTrainSetSize =
                        data->getNumCompoundsWithRestriction() * 0.75;
                const double k = pow((double) maxTrainSetSize, 1.0 / 6.0);
                conf.biosPoints[0] = 1;
                conf.biosPoints[1] = pow(k, 2);
                conf.biosPoints[2] = pow(k, 3);
                conf.biosPoints[3] = pow(k, 4);
                conf.biosPoints[4] = pow(k, 5);
                conf.biosPoints[5] = maxTrainSetSize;
            }
            Descriptor<double> *descriptor = new SineDescriptor<double>;
            descriptor->setSorting(SortingAlgorithm<double>::ROW_L1_DESC);
            Machine<double> qm2(descriptor, Kernel<double>::LAPLACE);
            qm2.setSigma(4.0e4);
            qm2.setLambda(1.0e-4);
            analyze<double>(qm2, data, conf);
            delete data;
            delete descriptor;
        }
    } catch (const char *e) {
        std::cerr << e << std::endl;
        std::cerr.flush();
        delete e;
    }
    return 0;
}

