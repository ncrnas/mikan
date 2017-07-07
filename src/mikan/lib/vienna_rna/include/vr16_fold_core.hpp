#ifndef VR16_FOLD_CORE_HPP_
#define VR16_FOLD_CORE_HPP_

#include "vr16_energy.hpp"                // VR16EnergyParams
#include "vr16_fold.hpp"                  // VR16Fold
#include "vr16_fold_options.hpp"          // VR16FoldOptions
#include "vr16_pair_mat.hpp"              // VR16PairMat
#include "vr16_params.hpp"                // VR16Params, VR16ParamIL, VR16PFParams
#include <vector>
#include <string>

namespace vr16 {

//
// Return values for ViennaRNA fold
//
class VR16FoldRet {
public:
    // Declare variables
    std::vector<std::string> mStructure;
    std::vector<double> mFoldEnergy;

public:
    // Define methods
    VR16FoldRet() : mVecSize(0) {}

    // Method prototypes
    void init_ret_vals(int pSize);

private:
    int mVecSize;

};

//
// Main fold calculation
//
class VR16FoldWorkSpace {
public:
    // Declare variables
    double mTemperature;

public:
    VR16FoldWorkSpace() :
            mTemperature(37.0), mEn(), mParams(mEn, mTemperature), mPairMat(), mOpts(),
            mPpIL(mEn, mParams, mPairMat, mOpts), mFold(mEn, mParams, mPpIL, mPairMat, mOpts), mFoldRet() {
        init_workspace();
    }

    VR16FoldWorkSpace(double pTemp) :
            mTemperature(pTemp), mEn(), mParams(mEn, mTemperature), mPairMat(), mOpts(),
            mPpIL(mEn, mParams, mPairMat, mOpts), mFold(mEn, mParams, mPpIL, mPairMat, mOpts), mFoldRet() {
        init_workspace();
    }

    void set_backtrack(bool pBT) { mOpts.mDoBacktrack = pBT; }

    void set_fold_constraint(bool pBT) { mOpts.mFoldConstrained = pBT; }

    void init_arrays(int pLength) { mFold.init_arrays(pLength); }

    double get_fold_energy(int pRetIdx) { return mFoldRet.mFoldEnergy[pRetIdx]; }

    const std::string &get_structure(int pRetIdx) { return mFoldRet.mStructure[pRetIdx]; }

    void preppare_fold(int pSize);

    int calc_fold_energy(int pRetIdx, std::string &pS1);

    void print_fold_ret_vals(int pRetIdx);


private:
    void init_workspace();

private:
    VR16EnergyParams mEn;
    VR16Params mParams;
    VR16PairMat mPairMat;
    VR16FoldOptions mOpts;
    VR16ParamIL mPpIL;

    VR16Fold mFold;
    VR16FoldRet mFoldRet;

};

} // namespace vr16

#endif /* VR16_FOLD_CORE_HPP_ */
