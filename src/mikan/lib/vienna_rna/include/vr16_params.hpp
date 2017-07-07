#ifndef VR16_PARAMS_HPP_
#define VR16_PARAMS_HPP_

#include "vr16_energy.hpp"                // VR16EnergyParams
#include "vr16_fold_options.hpp"          // VR16FoldOptions
#include "vr16_fold_options.hpp"          // VR16FoldOptions
#include "vr16_pair_mat.hpp"              // VR16PairMat
#include <cmath>
#include <map>

namespace vr16 {

//
// Scaled parameters for ViennaRNA fold & duplexfold
//
class VR16Params {
public:
    // Constant values
    static const int NBPAIRS = 7;
    static const int MAXLOOP = 30; // must be MAXLOOP >= 30
    const double K0;
    const double GASCONST;                      /* in [cal/K] */
    const int INF;
    const int MAX_NINIO;                        /* maximum correction */

    // Declare variables
    int mId;
    int **mStack;
    int *mHairpin;
    int *mBulge;
    int *mInternalLoop;
    int ***mMismatchI;
    int ***mMismatchH;
    int ***mMismatchM;
    int **mDangle5;
    int **mDangle3;
    int ****mInt11;
    int *****mInt21;
    int ******mInt22;
    int *mFNinio;
    double mlxc;
    int mMLBase;
    int *mMLintern;
    int mMLclosing;
    int mTerminalAU;
    int mDuplexInit;
    int *mTetraEnergy;
//    std::map<std::string, int> mTetraloops;
    int *mTriloopE;
//    std::map<std::string, int> mTriloops;
    double mTemperature;

public:
    // Define methods
    VR16Params(VR16EnergyParams &pEn, double pTemprature) :
            K0(pEn.K0), GASCONST(pEn.GASCONST), INF(pEn.INF), MAX_NINIO(pEn.MAX_NINIO), mId(-1), mEn(pEn) {
        init_heap();
        init_parameters(pTemprature);
    }

    ~VR16Params() {
        free_heap();
    }

private:
    VR16EnergyParams &mEn;

private:
    void init_heap();

    void free_heap();

    void init_parameters(double pTemprature);

};

//
// Cache minimum values of internal loops
//
class VR16ParamIL {

public:
    // Constant values
    static const int MAXLOOP = 30;
    static const int NBPAIRS = 7;

    // Define methods
    VR16ParamIL(VR16EnergyParams &pEn, VR16Params &pParams, VR16PairMat &pPairMat, VR16FoldOptions &pOpts) :
            mEn(pEn), mParams(pParams), mPairMat(pPairMat), mOpts(pOpts), INF(pParams.INF) {
        init_heap();
        init_parameters();
    }

    ~VR16ParamIL() {
        free_heap();
    }

    int loop_energy(int pN1, int pN2, int pType1, int pType2, int pSi1, int pSj1, int pSp1, int pSq1);

public:
    int **mMinIL;

private:
    VR16EnergyParams &mEn;
    VR16Params &mParams;
    VR16PairMat &mPairMat;
    VR16FoldOptions &mOpts;
    const int INF;
    static const int TURN = 3;

private:
    void init_heap();

    void free_heap();

    void init_parameters();
};

//
// Parameters of sub-optimal structure calculation
//
class VR16PFParams {
public:
    // Constant values
    static const int NBPAIRS = 7;
    static const int MAXLOOP = 30;
    static const int PF_SCALE_10 = 10;
    static const int INITLEN = 50;

    // Declare variables
    int mInitLength;                   /* length in last call to init_pf_fold() */

    std::vector<double> mScale;

    std::vector<double> mExpHairpin;
    std::vector<double> mExpMLbase;

    double *mExpBulge;
    double *mExpInternal;

    double **mExpNinio;

    double mExpMLclosing;
    double *mExpMLintern;

    double mExpTermAU;
    double **mExpDangle5;
    double **mExpDangle3;

    double *mExpTetra;
    double *mExpTriloop;
    double **mExpStack;

    double ***mExpMismatchI;
    double ***mExpMismatchH;
    double ***mExpMismatchM;

    double ****mExpInt11;
    double *****mExpInt21;
    double ******mExpInt22;

    double ***mExpMLContrib;

public:
    // Define methods
    VR16PFParams(VR16EnergyParams &pEn, double pTemperature, VR16Params &pParams, VR16FoldOptions &pOpts) :
            mInitLength(0), mExpMLclosing(0), mExpTermAU(0), mInitTemp(0), mEn(pEn), mParams(pParams),
            mOpts(pOpts) {
        init_heap();
        scale_pf_params(INITLEN, pTemperature);
    }

    ~VR16PFParams() {
        free_heap();
    }

    // Method prototypes
    void scale_pf_params(unsigned int pLen, double pTemperature);

    void reset_scale();

private:
    double mInitTemp;                  /* temperature in last call to scale_pf_params */

    VR16EnergyParams &mEn;
    VR16Params &mParams;
    VR16FoldOptions &mOpts;

private:
    void init_heap();

    void free_heap();

    double pf_smooth(double pX);
};

} // namespace vr16

#endif /* VR16_PARAMS_HPP_ */
