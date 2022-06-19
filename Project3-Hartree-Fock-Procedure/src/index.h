#ifndef INDEX_H
#define IINDEX_H


#include<vector>

static int BigNum = 10000;
static std::vector<int> ioff = [](int BigNum_) -> std::vector<int>{
    std::vector<int> ioff_(BigNum_);
    ioff_[0] = 0;
    for (int i = 1; i < BigNum_; i++)
    {
        ioff_[i] = i + ioff_[i-1];
    }
    return ioff_;
}(BigNum);//Thanks to Ding Li from SEU / NUS / NTU










#endif
#define INDEX(i,j) (i>j) ? (ioff[i]+j) : (ioff[j]+i)