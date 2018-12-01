/*
 * error.cc
 *
 *  Created on: 14.10.2016
 *      Author: klaus
 */

#include <iostream>
#include <chrono>
#include <fstream>

#include "machine.h"
#include "dataset.h"
#include "descriptor.h"
#include "progargparser.h"

struct RunConfiguration {
    std::string destDir;
    size_t *biosPoints;
    size_t numBiosPoints;
    size_t numRunsPerBiosPoint;
    double de;
};

enum SetType {
    VASP, VASP_FELIX, QM7
};

template<class T> void error(Machine<T> &qm, DataSet<T> *data,
        const RunConfiguration &conf) {
    const auto begin = std::chrono::steady_clock::now();
    std::ofstream fs_mae;
    fs_mae.open((conf.destDir + "/mae.dat").c_str(), std::ios::out);
    fs_mae << "# size" << std::setw(30) << "MAE" << std::setw(30)
            << "MAE per atom" << std::endl;
    for (size_t i = 0; i < conf.numBiosPoints; i++) {
        const size_t n = conf.biosPoints[i];
        std::cout << "---- Running " << conf.numRunsPerBiosPoint
                << " random set(s) with size " << n << " ----" << std::endl;

        std::ofstream fs_error;
        std::stringstream ss;
        ss << conf.destDir << "/error_" << n << ".dat";
        fs_error.open(ss.str().c_str(), std::ios::out);
        fs_error << "# start" << std::setw(30) << "end" << std::setw(30)
                << "count" << std::endl;

        std::ofstream fs_errorPerAtom;
        std::stringstream ssPerAtom;
        ssPerAtom << conf.destDir << "/errorPerAtom_" << n << ".dat";
        fs_errorPerAtom.open(ssPerAtom.str().c_str(), std::ios::out);
        fs_errorPerAtom << "# start" << std::setw(30) << "end" << std::setw(30)
                << "count" << std::endl;

        PredictionResults average;
        average.mae = 0.0;
        average.maePerAtom = 0.0;
        for (Prediction pred : average.results) {
            pred.error = 0.0;
            pred.errorPerAtom = 0.0;
            pred.value = 0.0;
            pred.valuePerAtom = 0.0;
        }
        for (size_t j = 0; j < conf.numRunsPerBiosPoint; j++) {
            data->setRandomSubSet(n);
            qm.train(data);
            data->setSubSetToAllExcludedElements();
            PredictionResults results;
            qm.predict(data, results);
            std::cout << "\tMAE: " << results.mae << " kcal/mol" << std::endl;
            std::cout << "\tMAE: " << results.mae * 0.043e3 << " meV"
                    << std::endl;
            std::cout << "\tMAE per atom: " << results.maePerAtom << " kcal/mol"
                    << std::endl;
            std::cout << "\tMAE per atom: " << results.maePerAtom * 0.043e3
                    << " meV" << std::endl;
            const auto end = std::chrono::steady_clock::now();
            std::cout << "\tTime: "
                    << std::chrono::duration_cast<std::chrono::seconds>(
                            end - begin).count() << "s" << std::endl;
            average.mae += results.mae;
            average.maePerAtom += results.maePerAtom;
            for (size_t i = 0; i < n; i++) {
                if (average.results.size() <= i)
                    average.results = results.results;
                else {
                    average.results[i].error += results.results[i].error;
                    average.results[i].errorPerAtom +=
                            results.results[i].errorPerAtom;
                    average.results[i].valuePerAtom +=
                            results.results[i].valuePerAtom;
                    average.results[i].value += results.results[i].value;
                }
            }
        }
        average.mae /= conf.numRunsPerBiosPoint;
        average.maePerAtom /= conf.numRunsPerBiosPoint;
        std::vector<double> sortedError;
        std::vector<double> sortedErrorPerAtom;
        for (Prediction pred : average.results) {
            pred.error /= conf.numRunsPerBiosPoint;
            pred.errorPerAtom /= conf.numRunsPerBiosPoint;
            sortedError.push_back(pred.error);
            sortedErrorPerAtom.push_back(pred.errorPerAtom);
            pred.value /= conf.numRunsPerBiosPoint;
            pred.valuePerAtom /= conf.numRunsPerBiosPoint;
        }
        std::cout << std::endl << "-------" << std::endl << "\tAverage MAE: "
                << average.mae << " kcal/mol" << std::endl;
        std::cout << "\tAverage MAE: " << average.mae * 0.043e3 << " meV"
                << std::endl;
        std::cout << "\tAverage MAE per atom: " << average.maePerAtom
                << " kcal/mol" << std::endl;
        std::cout << "\tAverage MAE per atom: " << average.maePerAtom * 0.043e3
                << " meV" << std::endl;
        const auto end = std::chrono::steady_clock::now();
        std::cout << "\tTime: "
                << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()
                << "s" << std::endl << std::endl;
        fs_mae << n << std::setw(30) << average.mae << std::setw(30)
                << average.maePerAtom << std::endl;
        std::sort(sortedError.begin(), sortedError.end(),
                [](double i,double j) {return ( i < j);});
        std::sort(sortedErrorPerAtom.begin(), sortedErrorPerAtom.end(),
                [](double i,double j) {return ( i < j);});
        std::map<int, int> histogram;
        for (const double error : sortedError) {
            int step = error / conf.de;
            std::map<int, int>::iterator it = histogram.find(step);
            if (it == histogram.end()) {
                histogram.insert(std::pair<int, int>(step, 1));
            } else {
                it->second++;
            }
        }
        for (const auto &ent1 : histogram) {
            fs_error << ent1.first * conf.de + conf.de / 2.0 << std::setw(30)
                    << ent1.first * conf.de - conf.de / 2.0 << std::setw(30)
                    << ent1.second << std::endl;
        }
        fs_error.close();

        std::map<int, int> histogramPerAtom;
        for (const double error : sortedErrorPerAtom) {
            int step = error / conf.de;
            std::map<int, int>::iterator it = histogramPerAtom.find(step);
            if (it == histogramPerAtom.end()) {
                histogramPerAtom.insert(std::pair<int, int>(step, 1));
            } else {
                it->second++;
            }
        }
        for (const auto &ent1 : histogramPerAtom) {
            fs_errorPerAtom << ent1.first * conf.de + conf.de / 2.0
                    << std::setw(30) << ent1.first * conf.de - conf.de / 2.0
                    << std::setw(30) << ent1.second << std::endl;
        }
        fs_errorPerAtom.close();
    }
    fs_mae.close();
}

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
    std::cout << "  --BiosPoints=<Array>:\ti.e. must be given in the format [1, 10, 100]" << std::endl;
    std::cout << "  --de=<Float>:\t" << std::endl;
    std::cout << "  --sorting=<Name>:\tComming soon" << std::endl;
    std::cout << "                   \tROW_L1_DESC: default" << std::endl;
    std::cout << "                   \tROW_L2_DESC:"<< std::endl;
    std::cout << "<Input File/Dir>:\t" << std::endl;
    std::cout << "<Output Dir>:\t" << std::endl;
}

int main(int argc, char** argv) {
    SetType type;
    RunConfiguration conf;
    conf.de = 0.5;
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
            std::cout << "numRunsPerBiosPoint: " << runsPerBiosPoint
                    << std::endl;
        } else {
            std::cerr << "--RunsPerBiosPoint has invalid value!" << std::endl;
            printHelp(argv[0]);
            return 1;
        }
    }
    if (programArgs.cmdOptionExists("--de")) {
        double de = programArgs.getDoubleValue("--de");
        if (de > 0) {
            conf.de = de;
        } else {
            std::cerr << "--de has invalid value!" << std::endl;
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
    if (defaultConfForBiosPoints) {
        conf.biosPoints = new size_t[7];
        conf.biosPoints[0] = 1;
        conf.biosPoints[1] = 10;
        conf.biosPoints[2] = 100;
        conf.biosPoints[3] = 500;
        conf.biosPoints[4] = 1000;
        conf.biosPoints[5] = 3000;
        //conf.biosPoints[6] = 4000;
        conf.numBiosPoints = 6;
    }
    try {
        if (type == QM7) {
            DataSet<float> *data = new QMDataSet_7(inputFile);
            Descriptor<float> *descriptor = new CoulombDescriptor<float>;
            descriptor->setSorting(SortingAlgorithm<float>::ROW_L1_DESC);
            Machine<float> qm(descriptor, Kernel<float>::LAPLACE);
            qm.setSigma(2.5e3);
            qm.setLambda(1.0e-6);
            error<float>(qm, data, conf);
            delete data;
            delete descriptor;
        } else if (type == VASP || type == VASP_FELIX) {
            VASPDataSet *data = new VASPDataSet(inputFile,
                    type == VASP_FELIX ? true : false);
            Descriptor<double> *descriptor = new SineDescriptor<double>;
            descriptor->setSorting(SortingAlgorithm<double>::ROW_L1_DESC);
            Machine<double> qm2(descriptor, Kernel<double>::LAPLACE);
            qm2.setSigma(4.0e4);
            qm2.setLambda(1.0e-4);
            error<double>(qm2, data, conf);
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

