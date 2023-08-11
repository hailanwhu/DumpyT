#include <iostream>
#include <cstdio>
#include <vector>
#include <thread>
#include <chrono>
#include <set>
#include <random>
#include "../include/DataStructures/DumpyNode.h"
#include "../include/DataStructures/GraphConstruction.h"
#include "../include/Searchers/DumpySearcher.h"
#include "../include/Searchers/DumpyInMemorySearch.h"
#include "../include/DataStructures/FullAryTreeNode.h"
#include "../include/Utils/FileUtil.h"
#include "../include/Utils/MathUtil.h"
#include "../include/Utils/TimeSeriesUtil.h"
#include "../include/Utils/SaxUtil.h"
using namespace std;

vector<vector<int>>* loadGraphSkeleton(){
    int vd = 0;
    for(int i=1; i<=Const::bitsReserve;++i)
        vd += MathUtil::nChooseK(Const::segmentNum, i);
    auto nnList = new vector<vector<int>>(Const::vertexNum, vector<int>(vd, -1));

    if(!FileUtil::checkFileExists(Const::graphfn.c_str())){
        cout << "File not exists!" << Const::graphfn << endl;
        exit(-1);
    }
    FILE *f = fopen(Const::graphfn.c_str(), "rb");

    for(int i=1;i<Const::vertexNum;++i)
        fread(&((*nnList)[i][0]), sizeof(int), vd, f);

    return nnList;

}

void constructGraph(){
    GraphConstruction::buildAndSave2Disk();
}

void buildInMemoryIndexDumpy(){
    FullAryTreeNode* root = FullAryTreeNode::BuildFullAryTree();
    cout << "build finish" << endl;
    root->getIndexStats();
    root->save2Disk(Const::memoryidxfn);
}

void searchInMemoryThreshold() {
    FullAryTreeNode* root = FullAryTreeNode::loadFromDisk(Const::memoryidxfn);
    cout << "load finish." << endl;
    float *queries = FileUtil::readQueries();
    for(int i=0;i<Const::query_num;++i){
        Const::logPrint("Query " + to_string(i) +":");
        cout << "Threshold: " << Const::threshold << endl;
        long ts_count = 0;
        vector<PqItemSeries*> *approxKnn = DumpyInMemorySearch::thresholdSearch(root, queries + i * Const::tsLength, Const::threshold, &ts_count);
        Const::logPrint("Results:");
        cout << "Result Size: " << approxKnn->size() << endl;
        cout << "Calculate TS: " << ts_count << endl;
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << (*approxKnn)[j]->dist << "------" << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
//        for (int j = 0; j < approxKnn->size(); ++j) {
//            cout << j + 1 << ": " << (*approxKnn)[j]->dist << endl;
//        }
    }
}

void searchInMemory(){
    FullAryTreeNode* root = FullAryTreeNode::loadFromDisk(Const::memoryidxfn);
    cout << "load finish." << endl;
    float *queries = FileUtil::readQueries();
    cout << fixed <<setprecision(6);
    int search_num  = Const::series_num * Const::search_ratio;
    for(int i=0;i<Const::query_num;++i){
        Const::logPrint("Query " + to_string(i) +":");
        vector<PqItemSeries*> *approxKnn = DumpyInMemorySearch::approxSearch(root, queries + i * Const::tsLength, Const::k, search_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << (*approxKnn)[j]->dist << "------" << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}

void tsProfile() {
    float *queries = FileUtil::readQueries();
    cout << fixed <<setprecision(6);
    float *q1 = queries;
    float *q2 = queries + Const::tsLength;
    auto* t1 = new TimeSeries(q1);
    auto* t2 = new TimeSeries(q2);
    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(t2->sax))[i];
    cout << SaxUtil::LowerBound_Paa_iSax(t1->paa,sax, Const::bitsCardinality) << endl;
}

void statMemoryDumpy(){
    FullAryTreeNode* root = FullAryTreeNode::loadFromDisk(Const::memoryidxfn);
    cout << "load finish" << endl;
    root->getIndexStats();
}

void buildDumpy(){
    auto g = loadGraphSkeleton();
    DumpyNode* root = DumpyNode::BuildIndex(Const::datafn, Const::saxfn);
    root->save2Disk(Const::idxfn + "root.idx");
}

void approxSearchOneNode() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxSearch(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::idxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}

void approxSearchMoreNode() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxIncSearch(root, queries + i * Const::tsLength,
                                                                           Const::k, Const::idxfn,
                                                                           Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void approxSearchOneNodeDTW() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxSearchDTW(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::idxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}

void approxSearchMoreNodeDTW() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxIncSearchDTW(root, queries + i * Const::tsLength,
                                                                           Const::k, Const::idxfn,
                                                                           Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void buildDumpyFuzzy(){
    auto g = loadGraphSkeleton();
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode* root = DumpyNode::BuildIndexFuzzy(Const::datafn, Const::saxfn, Const::paafn, g);
    root->save2Disk(Const::fuzzyidxfn + "root.idx");
}

void approxSearchOneNodeFuzzy() {
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxSearch(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::fuzzyidxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void approxSearchMoreNodeFuzzy() {
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        auto start = chrono::system_clock::now();
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxIncSearchFuzzy(root, queries + i * Const::tsLength,
                                                                                Const::k, Const::fuzzyidxfn,
                                                                                Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void exactExprDumpy() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        auto start = chrono::system_clock::now();
        vector<PqItemSeries *> *exactKnn = DumpySearcher::exactSearch(root, queries + i * Const::tsLength, Const::k, g);
        Const::logPrint("Results:");
        for (int j = 0; j < exactKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*exactKnn)[j]->ts) << endl;
    }
}

void statIndexDumpy(){
    DumpyNode* root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    root->getIndexStats();
}

void storeTS() {
    FullAryTreeNode* root = FullAryTreeNode::loadFromDisk(Const::memoryidxfn);
    cout << "load finish." << endl;
    root->storeTSInLeafNode();
}

void statIndexDumpyFuzzy(){
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode* root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    root->getIndexStats();
}

std::set<int> generateDistinctIntegers(int count, int minRange, int maxRange) {
    std::set<int> distinctIntegers;

    // Check if the desired count is within the range
    if (count > (maxRange - minRange + 1)) {
        std::cout << "Cannot generate " << count << " distinct integers within the given range." << std::endl;
        return distinctIntegers;
    }

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> distribution(minRange, maxRange);

    while (distinctIntegers.size() < count) {
        int randomNumber = distribution(generator);
        distinctIntegers.insert(randomNumber);
    }

    return distinctIntegers;
}

void genTrainingTS(){
    FullAryTreeNode* root = FullAryTreeNode::loadFromDisk(Const::memoryidxfn);
    long f_size = FileUtil::getFileSize(Const::saxfn.c_str()), num = f_size / (sizeof(unsigned short) * Const::segmentNum);
    cout << "total ts count: " << num << endl;

    std::set<int> distinctIntegers = generateDistinctIntegers(2000000, 0, num);
    float *sampled = new float [2000000 * Const::tsLength];
    int i = 0;
    for (int offset : distinctIntegers) {
        copy(root->dataset + (long)offset * Const::tsLength,
             root->dataset + (long)offset * Const::tsLength + Const::tsLength,
             sampled + (long)i * Const::tsLength);
        i += 1;
    }
    cout << "i: " << i << endl;
    FILE *outf = fopen("sampled_2m.bin", "a");
    fwrite(sampled, sizeof(float), 2000000 * Const::tsLength, outf);
}

int main() {
    Const::readConfig();

    switch (Const::index) {
        case 0:
            constructGraph();
            break;
        case 1:
            switch (Const::ops) {
                case 0:
                    buildDumpy();
                    break;
                case 1:
                    approxSearchOneNode();
                    break;
                case 2:
                    exactExprDumpy();
                    break;
                case 3:
                    statIndexDumpy();
                    break;
                case 4:
                    approxSearchMoreNode();
                    break;
                case 5:
                    approxSearchOneNodeDTW();
                    break;
                case 6:
                    approxSearchMoreNodeDTW();
                    break;
                default:
                    break;
            }
            break;
        case 2:
            if(Const::ops == 0){
                buildDumpyFuzzy();
            }
            else if(Const::ops == 1){
                approxSearchOneNodeFuzzy();
            }else if(Const::ops == 3){
                statIndexDumpyFuzzy();
            }else if(Const::ops == 4){
                approxSearchMoreNodeFuzzy();
            }
            break;
        case 3:
            if(Const::ops == 0) buildInMemoryIndexDumpy();
            else if(Const::ops == 1) searchInMemory();
            else if(Const::ops == 4)    statMemoryDumpy();
            else if(Const::ops == 7) searchInMemoryThreshold();
            else if(Const::ops == 8) storeTS();
            break;
        case 4:
            SaxUtil::generateSaxFile(Const::datafn, Const::saxfn);
            SaxUtil::generatePaaFile(Const::datafn, Const::paafn);
            break;
        case 5:
            tsProfile();
            break;
        case 6:
            genTrainingTS();
            break;
        default:    break;
    }
}
