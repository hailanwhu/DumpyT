//
// Created by Zeyu Wang on 2021/8/7.
//

#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Utils/FileUtil.h"
#include <bitset>
#define distPeng(x,y) ((x-y)*(x-y))

using namespace std;

void TimeSeriesUtil::heap_data_copy(vector<PqItemSeries *> &heap){
    for(auto *pis:heap)
        pis->copyData();
}

bool TimeSeriesUtil::isSame(const PqItemSeries *ts1, const float *ts2)
{

    float *t1 = (ts1->ts);
    for(int i = 0; i < Const::tsLength; ++i) {
        if(abs(t1[i] - ts2[i])>1e-5)
            return false;
    }
    return true;;
}

int TimeSeriesUtil::intersectionTsSets(const vector<PqItemSeries *> *tsSet1, vector<float *> *tsSet2){
    int intersectionNum = 0, size= tsSet1->size();
    for(const PqItemSeries* currentTs : *tsSet1)
        for(int i=0;i<size;++i){
            if(isSame(currentTs, (*tsSet2)[i])) {
                intersectionNum++;
                break;
            }
        }

    return intersectionNum;
}

int TimeSeriesUtil::intersectionTsSetsCardinality(const vector<PqItemSeries *> &tsSet1, const vector<PqItemSeries *> &tsSet2){
    int intersectionNum = 0;
    for(PqItemSeries* currentTs : tsSet1)
        for (PqItemSeries* targetTs : tsSet2)
            if (isSame(currentTs, targetTs->ts)) {
                intersectionNum += 1;
                break;
            }
    return intersectionNum;
}


double TimeSeriesUtil::euclideanDist(const float* ts_1, const float* ts_2, int len) {
    double sum = 0, dp;
    for (int i = 0; i < len; i++) {
        dp = ts_1[i] - ts_2[i];
        sum += dp * dp;
    }
    return sum;
}

double TimeSeriesUtil::dtw(const float* A, const float* B, int len, int r, double bsf)
{

    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=numeric_limits<double>::max();

    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=numeric_limits<double>::max();

    for (i=0; i<len; i++)
    {
        k = max(0,r-i);
        min_cost = numeric_limits<double>::max();

        for(j=max(0,i-r); j<=min(len-1,i+r); j++, k++)
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[k]=distPeng(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = numeric_limits<double>::max();
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = numeric_limits<double>::max();
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = numeric_limits<double>::max();
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + distPeng(A[i],B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (i+r < len-1 && min_cost  >= bsf)
        {   free(cost);
            free(cost_prev);
            return min_cost ;
        }
        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}


double TimeSeriesUtil::euclideanDist(const float* ts_1, const float* ts_2, int len, double bound) {
    double sum = 0, dp;
    for (int i = 0; i < len && sum < bound; i++) {
        dp = ts_1[i] - ts_2[i];
        sum += dp * dp;
    }
    return sum;
}

template<typename ... Args>
string TimeSeriesUtil::str_format(const string &format, Args ... args)
{
    auto size_buf = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1;
    std::unique_ptr<char[]> buf(new(std::nothrow) char[size_buf]);

    if (!buf)
        return string{};

    std::snprintf(buf.get(), size_buf, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size_buf - 1);
}