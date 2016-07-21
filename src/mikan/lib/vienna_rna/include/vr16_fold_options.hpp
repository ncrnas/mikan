#ifndef VR16_FOLD_OPTIONS_HPP_
#define VR16_FOLD_OPTIONS_HPP_

#include <vector>

namespace vr16 {

//
// Store base pairs
//
class VR16Bond {
public:
    // Declare variables
    int i;
    int j;

public:
    // Declare variables
    VR16Bond(): i(0), j(0) {}
};

//
// Options for ViennaRNA fold, duplexfold, and part_func
//
class VR16FoldOptions
{
public:
    // Declare variables
    bool mNoGU;                                    /* GU not allowed at all */
    bool mNoClosingGU;                             /* GU allowed only inside stacks */
    bool mTetraLoop;                               /* Fold with specially stable 4-loops */
    bool mTriLoop;
    int mEnergySet;                                /* 0 = BP; 1=any with GC; 2=any with AU parameters */
    int mDangles;                                  /* use dangling end energies */
    std::string mNonStandards;                     /* contains allowed non standard bases */
    double mTemperature;
    bool mJamesRule;                               /* interior loops of size 2 get energy 0.8Kcal and
                                                      no mismatches (no longer used) */
    int mLogML;                                    /* use logarithmic multiloop energy function */
    int mCutPoint;                                 /* first position of 2nd strand for co-folding */

    std::vector<VR16Bond > mBasePair;              /* list of base pairs */

    std::vector<double> mProb;                     /* base pairing prob. matrix */
    std::vector<int>  mIIndx;                      /* mProb[i,j] -> mProb[mIIndx[i]-j] */

    double mPfScale;                               /* scaling factor to avoid float overflows*/
    bool mFoldConstrained;                         /* fold with constraints */
    bool mDoBacktrack;                             /* calculate pair prob matrix in part_func() */
    bool mNoLonelyPairs;                           /* avoid helices of length 1 */
    char mBacktrackType;                           /* usually 'F'; 'C' require (1,N) to be bonded;
                                                      'M' seq is part of a multi loop */

public:
    // Define methods
    VR16FoldOptions() :
        mNoGU(false), mNoClosingGU(false), mTetraLoop(true), mTriLoop(false), mEnergySet(0),
        mDangles(1), mNonStandards(""), mTemperature(37.0), mJamesRule(true), mLogML(0), mCutPoint(0),
        mPfScale(-1), mFoldConstrained(false), mDoBacktrack(false), mNoLonelyPairs(false), mBacktrackType('F')
    {}
};

} // namespace vr16

#endif /* VR16_FOLD_OPTIONS_HPP_ */
