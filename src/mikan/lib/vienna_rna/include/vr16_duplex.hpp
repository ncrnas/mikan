#ifndef VR16_DUPLEX_HPP_
#define VR16_DUPLEX_HPP_

#include <mikan/lib/vienna_rna/include/vr16_energy.hpp>                // VR16EnergyParams
#include <mikan/lib/vienna_rna/include/vr16_fold_options.hpp>          // VR16FoldOptions
#include <mikan/lib/vienna_rna/include/vr16_pair_mat.hpp>              // VR16PairMat
#include <mikan/lib/vienna_rna/include/vr16_params.hpp>                // VR16Params, VR16ParamIL
#include <string>

namespace vr16 {

//
// ViennaRNA duplexfold
//
class VR16Duplex
{
public:
    // Constant values
    static const int FORBIDDEN = 9999;
    static const int BONUS = 10000;
    static const int TURN = 3;

    // Declare variables
    int mI;
    int mJ;
    std::string mStructure;
    int mL1;
    float mEnergy;
    int mE5p;
    int mE3p;

public:
    // Define methods
    VR16Duplex(VR16EnergyParams& pEn, VR16Params& pParams, VR16ParamIL& pPpIL, VR16PairMat& pPairMat,
            VR16FoldOptions& pOpts) :
        mI(0), mJ(0), mL1(0), mEnergy(0), mE5p(0), mE3p(0),
        mN1(0), mN2(0), mBonusGiven(0),
        mFivePrimeLength(0), mAutoFivePrimeLength(1), mDebugMode(0), mDelayFree(0),
        mEn(pEn), mParams(pParams), mPpIL(pPpIL), mPairMat(pPairMat), mOpts(pOpts)
    {}

    // Method prototypes
    int duplexfold(std::string &pS1, std::string &pS2, std::vector<int> &pArrayI, std::vector<int> &pArrayJ);

private:
    static const int STACK_BULGE1 = 1;   /* stacking energies for bulges of size 1 */
    static const int NEW_NINIO = 1;      /* new asymetry penalty */
    static const int MAXSECTORS = 500;   /* dimension for a backtrack array */
    static const int LOCALITY = 0;       /* locality parameter for base-pairs */
    static const int BONUS_SIZE = 10000;

    std::vector<std::vector<int> > mC;    /* energy array, given that i-j pair */
    std::vector<int> mS1;
    std::vector<int> mSS1;
    std::vector<int> mS2;
    std::vector<int> mSS2;

    int mN1;                               /* sequence lengths */
    int mN2;                               /* sequence lengths */

    int mBonusGiven;

    int mFivePrimeLength;
    int mAutoFivePrimeLength;

    int mDebugMode;
    int mDelayFree;

    VR16EnergyParams& mEn;
    VR16Params& mParams;
    VR16ParamIL& mPpIL;
    VR16PairMat& mPairMat;
    VR16FoldOptions& mOpts;

private:
    void encode_seq(std::string &pS1, std::string& pS2);
    int backtrack(int pI, int pJ, std::vector<int> &pArrayI, std::vector<int> &pArrayJ);

};

} // namespace vr16

#endif /* VR16_DUPLEX_HPP_ */
