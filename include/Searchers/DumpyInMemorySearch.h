//
// Created by Zeyu Wang on 2021/9/24.
//

#ifndef DUMPY_DUMPYINMEMORYSEARCH_H
#define DUMPY_DUMPYINMEMORYSEARCH_H


#include <vector>
#include "../DataStructures/PqItemSeries.h"
#include "../DataStructures/FullAryTreeNode.h"

class DumpyInMemorySearch {
public:
    static vector<PqItemSeries *> *approxSearch(FullAryTreeNode *root, float *query, int k, int threshold);

    static vector<PqItemSeries *> *thresholdSearch(FullAryTreeNode *root, float *query, double threshold, long *ts_count);

    static void approxSearchSub(FullAryTreeNode *root, unsigned short *sax, TimeSeries *queryTs, int k, int threshold,
                                vector<PqItemSeries *> *heap);

    static void thresholdSearchSub(FullAryTreeNode* root, unsigned short *sax, TimeSeries* queryTs, double threshold, vector<PqItemSeries*>*heap, long *ts_count);

};


#endif //DUMPY_DUMPYINMEMORYSEARCH_H
