#include "vr16_fold_core.hpp"
#include <iostream>

namespace vr16 {

//
// VR16FoldRet methods
//
void VR16FoldRet::init_ret_vals(int pSize) {
    if (pSize == mVecSize) {
        return;
    }

    mStructure.resize((unsigned) pSize);
    mFoldEnergy.resize((unsigned) pSize);

    mVecSize = pSize;
}

//
// VR16FoldWorkSpace methods
//
void VR16FoldWorkSpace::init_workspace() {
    mPairMat.make_pair_matrix(mOpts.mEnergySet, mOpts.mNonStandards, mOpts.mNoGU);
}

void VR16FoldWorkSpace::preppare_fold(int pSize) {
    mFoldRet.init_ret_vals(pSize);
}

void VR16FoldWorkSpace::print_fold_ret_vals(int pRetIdx) {
    std::cout << "Structure:     " << mFoldRet.mStructure[pRetIdx] << std::endl;
    std::cout << "Fold energy:   " << mFoldRet.mFoldEnergy[pRetIdx] << std::endl;
}

int VR16FoldWorkSpace::calc_fold_energy(int pRetIdx, std::string &pString) {
    double en;
    mFoldRet.mStructure[pRetIdx].resize(pString.size(), 0);

    en = mFold.fold(pString, mFoldRet.mStructure[pRetIdx]);
    mFoldRet.mFoldEnergy[pRetIdx] = en;

    return 0;
}


} // namespace vr16
