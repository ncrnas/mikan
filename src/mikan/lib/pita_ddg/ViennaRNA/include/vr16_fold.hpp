#ifndef VR16_FOLD_HPP_
#define VR16_FOLD_HPP_

#include <mikan/lib/pita_ddg/ViennaRNA/include/vr16_energy.hpp>                // VR16EnergyParams
#include <mikan/lib/pita_ddg/ViennaRNA/include/vr16_fold_options.hpp>          // VR16FoldOptions
#include <mikan/lib/pita_ddg/ViennaRNA/include/vr16_pair_mat.hpp>              // VR16PairMat
#include <mikan/lib/pita_ddg/ViennaRNA/include/vr16_params.hpp>                // VR16Params, VR16ParamIL
#include <vector>
#include <string>

namespace vr16 {

//
// Partial structure for backtracking
//
struct VR16Sect
{
    int i;
    int j;
    int ml;
};

//
// Store C Array for ViennaRNA fold
//
class VR16FoldArrayC
{
public:
    // Constant values
    static const int LOCALITY = 0;          /* locality parameter for base-pairs */
    static const int TURN = 3;
    static const int FORBIDDEN = 9999;
    static const int BONUS = 10000;

    std::vector<int> mC;                    /* energy array, given that i-j pair */

public:
    // Define methods
    VR16FoldArrayC(VR16EnergyParams& pEn, VR16Params& pParams, VR16ParamIL& pPpIL, VR16PairMat& pPairMat,
            VR16FoldOptions& pOpts, std::vector<int>& pIdx) :
        mInitLength(-1), mIdDMLi2(-2), mIdDMLi1(1), mIdDMLi(0), mIdCC1(0), mIdCC(0),
        mEn(pEn), mParams(pParams), mPpIL(pPpIL), mPairMat(pPairMat), mOpts(pOpts), mIdex(pIdx)
    {}
    void set_dmli_val(int pJ, int decomp){set_dmlx_val(mIdDMLi, pJ, decomp);}

    // Method prototypes
    int init_arrays(int pStrLen);
    void fill_arrays(int pI, int pJ, std::vector<int>& pBP, std::vector<char>& pPtype, std::vector<int>& pS1,
            std::vector<int>& pFML, std::string &pString);
    void rotate_aux_idexes(int pStrLen);
    int hairpin_e(int pSize, int pType, int pSi1, int pSj1, int pIdx, std::string &pString);

private:
    std::vector<int> mCC;                   /* linear array for calculating canonical structures */
    std::vector<int> mCC1;                  /*   "     "        */

    std::vector<int> mDMLi;                 /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
    std::vector<int> mDMLi1;                /*             MIN(fML[i+1,k]+fML[k+1,j])  */
    std::vector<int> mDMLi2;                /*             MIN(fML[i+2,k]+fML[k+1,j])  */

    int mInitLength;

    int mIdDMLi2;
    int mIdDMLi1;
    int mIdDMLi;

    int mIdCC1;
    int mIdCC;

    VR16EnergyParams& mEn;
    VR16Params& mParams;
    VR16ParamIL& mPpIL;
    VR16PairMat& mPairMat;
    VR16FoldOptions& mOpts;
    std::vector<int>& mIdex;

private:
    void resize_arrays(int pStrLen);
    void clear_arrays();

    int get_init_bonus(int pI, int pJ, std::vector<int>& pBP);
    int get_type_1(int pI, int pJ, int pStrLen, std::vector<char>& pPtype, std::vector<int>& pBP);
    int calc_ml_decomp(int pI, int pJ, int pType1, std::vector<int>& pS1);
    int calc_coaxial_stack(int pI, int pJ, int pType1, std::vector<char>& pPtype, std::vector<int>& pFML);
    void calc_elementary_struct(int pI, int pJ, int pType1, int pNoClose, std::vector<char>& pPtype,
            std::vector<int>& pS1, int& pNewC, int& pStackEnergy);
    void set_c(int pI, int pJ, int pBonus, int pNewC, int pStackEnergy);

    int get_dmlx_val(int pDMLId, int pIdx);
    int get_ccx_val(int pCCId, int pIdx);
    void set_dmlx_val(int pDMLId, int pIdx, int pNewVal);
    void set_ccx_val(int pCCId, int pIdx, int pNewVal);

};

//
// Store ML array for ViennaRNA fold
//
class VR16FoldArrayML
{
public:
    // Constant values
    static const int TURN = 3;

    // Declare variables
    bool mUniqML;                           /* do ML decomposition uniquely (for subopt) */
    std::vector<int> mFML;                  /* multi-loop auxiliary energy array */

public:
    // Define methods
    VR16FoldArrayML(VR16EnergyParams& pEn, VR16Params& pParams, VR16ParamIL& pPpIL, VR16PairMat& pPairMat,
            VR16FoldOptions& pOpts, std::vector<int>& pIdx) :
        mUniqML(false), mInitLength(-1),
        mEn(pEn), mParams(pParams), mPpIL(pPpIL), mPairMat(pPairMat), mOpts(pOpts), mIdex(pIdx)
    {}

    // Method prototypes
    int init_arrays(int pStrLen);
    void fill_arrays(int pI, int pJ, int pStrLen, std::vector<char>& pPtype, std::vector<int>& pS1,
            VR16FoldArrayC& pArrayC);
    void reset_fmi(int pStrLen);

private:
    std::vector<int> mFM1;                  /* second ML array, only for subopt */
    std::vector<int> mFmi;                  /* holds row i of fML (avoids jumps in memory) */

    int mInitLength;

    VR16EnergyParams& mEn;
    VR16Params& mParams;
    VR16ParamIL& mPpIL;
    VR16PairMat& mPairMat;
    VR16FoldOptions& mOpts;
    std::vector<int>& mIdex;

private:
    void resize_arrays(int pSize);
    void clear_arrays();
};

//
// Process backtrack for ViennaRNA fold
//
class VR16FoldBackTrack
{
public:
    // Constant values
    static const int MAXSECTORS = 500;      /* dimension for a backtrack array */
    static const int TURN = 3;
    static const int FORBIDDEN = 9999;
    static const int BONUS = 10000;

public:
    // Define methods
    VR16FoldBackTrack(VR16EnergyParams& pEn, VR16Params& pParams, VR16ParamIL& pPpIL, VR16PairMat& pPairMat,
            VR16FoldOptions& pOpts, std::vector<int>& pIdx) :
        mEn(pEn), mParams(pParams), mPpIL(pPpIL), mPairMat(pPairMat), mOpts(pOpts), mIdex(pIdx)
    {}

    // Method prototypes
    void backtrack(std::string &pString, int pS, std::vector<char>& pPtype, std::vector<int>& pS1,
            std::vector<int>& pBP, VR16FoldArrayC& pArrayC, VR16FoldArrayML& pArrayML,
            std::vector<int>& pF5);
    void parenthesis_structure(std::string &pStructure, int pStrLen);

private:
    VR16Sect mSector[MAXSECTORS];           /* stack of partial structures for backtracking */

    VR16EnergyParams& mEn;
    VR16Params& mParams;
    VR16ParamIL& mPpIL;
    VR16PairMat& mPairMat;
    VR16FoldOptions& mOpts;
    std::vector<int>& mIdex;
};

//
// ViennaRNA fold
//
class VR16Fold
{
public:
    // Constant values
    static const int MAXSECTORS = 500;      /* dimension for a backtrack array */
    static const int TURN = 3;
    static const int FORBIDDEN = 9999;
    static const int BONUS = 10000;

public:
    // Define methods
    VR16Fold(VR16EnergyParams& pEn, VR16Params& pParams, VR16ParamIL& pPpIL, VR16PairMat& pPairMat,
            VR16FoldOptions& pOpts) :
        mInitLength(-1), mAlpha("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"),
        mEn(pEn), mParams(pParams), mPpIL(pPpIL), mPairMat(pPairMat), mOpts(pOpts), mIdex(),
        mArrayC(pEn, pParams, pPpIL, pPairMat, pOpts, mIdex), mArrayML(pEn, pParams, pPpIL, pPairMat, pOpts, mIdex),
        mBT(pEn, pParams, pPpIL, pPairMat, pOpts, mIdex)
    {}

    // Method prototypes
    float fold(std::string &pString, std::string &pStructure);
    int init_arrays(int pLength);

protected:
    std::vector<int> mF5;                   /* energy of 5' end */
    std::vector<char> mPtype;               /* precomputed array of pair types */
    std::vector<int> mS;
    std::vector<int> mS1;
    std::vector<int> mPairTable;

    int mInitLength;
    std::string mAlpha;
    std::vector<int> mBP;                   /* contains the structure constrainsts: BP[i]
                                               -1: | = base must be paired
                                               -2: < = base must be paired with j<i
                                               -3: > = base must be paired with j>i
                                               -4: x = base must not pair
                                               positive int: base is paired with int      */

    VR16EnergyParams& mEn;
    VR16Params& mParams;
    VR16ParamIL& mPpIL;
    VR16PairMat& mPairMat;
    VR16FoldOptions& mOpts;
    std::vector<int> mIdex;
    VR16FoldArrayC mArrayC;
    VR16FoldArrayML mArrayML;
    VR16FoldBackTrack mBT;


protected:
    void resize_arrays(int pStrLen);
    void clear_arrays();

    void encode_seq(std::string &pString);
    void make_ptypes(std::string &pStructure, std::vector<char>& pPtype, std::vector<int>& pBP,
            VR16FoldOptions& pOpts);
    int fill_arrays(std::string &pString);

    void calc_energy_of_fragments(int pStrLen, std::vector<char>& pPtype,
            std::vector<int>& pF5, VR16FoldOptions& pOpts);

    double finalize_energy_calc(int pStrLen, std::string &pStructure, int pEnergy, VR16FoldOptions& pOpts,
            std::vector<int>& pBP);

};

//
// VR16FoldArrayC inline methods
//
inline int VR16FoldArrayC::get_dmlx_val(int pDMLId, int pIdx)
{
    int retVal = 0;

    if (pDMLId == 0)
    {
        retVal = mDMLi[pIdx];
    }
    else if (pDMLId == 1)
    {
        retVal = mDMLi1[pIdx];
    }
    else if (pDMLId == 2)
    {
        retVal = mDMLi2[pIdx];
    }

    return retVal;
}

inline void VR16FoldArrayC::set_dmlx_val(int pDMLId, int pIdx, int pNewVal)
{
    if (pDMLId == 0)
    {
        mDMLi[pIdx] = pNewVal;
    }
    else if (pDMLId == 1)
    {
        mDMLi1[pIdx] = pNewVal;
    }
    else if (pDMLId == 2)
    {
        mDMLi2[pIdx] = pNewVal;
    }
}

inline int VR16FoldArrayC::get_ccx_val(int pCCId, int pIdx)
{
    int retVal = 0;

    if (pCCId == 0)
    {
        retVal = mCC[pIdx];
    }
    else if (pCCId == 1)
    {
        retVal = mCC1[pIdx];
    }

    return retVal;
}

inline void VR16FoldArrayC::set_ccx_val(int pCCId, int pIdx, int pNewVal)
{
    if (pCCId == 0)
    {
        mCC[pIdx] = pNewVal;
    }
    else if (pCCId == 1)
    {
        mCC1[pIdx] = pNewVal;
    }
}

} // namespace vr16

#endif /* VR16_FOLD_HPP_ */
