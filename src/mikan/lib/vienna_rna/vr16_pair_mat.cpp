#include "vr16_pair_mat.hpp"
#include <iostream>
#include "mk_memory.hpp"

namespace vr16 {

//
// VR16PairMat methods
//
void VR16PairMat::init_dat() {
    mLawAndOrder = "_ACGUTXKI";

/*     _  A  C  G  U  X  K  I
    {{ 0, 0, 0, 0, 0, 0, 0, 0},
     { 0, 0, 0, 0, 5, 0, 0, 5},
     { 0, 0, 0, 1, 0, 0, 0, 0},
     { 0, 0, 2, 0, 3, 0, 0, 0},
     { 0, 6, 0, 4, 0, 0, 0, 6},
     { 0, 0, 0, 0, 0, 0, 2, 0},
     { 0, 0, 0, 0, 0, 1, 0, 0},
     { 0, 6, 0, 0, 5, 0, 0, 0}}; */

    mBPPair[0][0] = 0;
    mBPPair[0][1] = 0;
    mBPPair[0][2] = 0;
    mBPPair[0][3] = 0;
    mBPPair[0][4] = 0;
    mBPPair[0][5] = 0;
    mBPPair[0][6] = 0;
    mBPPair[0][7] = 0;

    mBPPair[1][0] = 0;
    mBPPair[1][1] = 0;
    mBPPair[1][2] = 0;
    mBPPair[1][3] = 0;
    mBPPair[1][4] = 5;
    mBPPair[1][5] = 0;
    mBPPair[1][6] = 0;
    mBPPair[1][7] = 5;

    mBPPair[2][0] = 0;
    mBPPair[2][1] = 0;
    mBPPair[2][2] = 0;
    mBPPair[2][3] = 1;
    mBPPair[2][4] = 0;
    mBPPair[2][5] = 0;
    mBPPair[2][6] = 0;
    mBPPair[2][7] = 0;

    mBPPair[3][0] = 0;
    mBPPair[3][1] = 0;
    mBPPair[3][2] = 2;
    mBPPair[3][3] = 0;
    mBPPair[3][4] = 3;
    mBPPair[3][5] = 0;
    mBPPair[3][6] = 0;
    mBPPair[3][7] = 0;

    mBPPair[4][0] = 0;
    mBPPair[4][1] = 6;
    mBPPair[4][2] = 0;
    mBPPair[4][3] = 4;
    mBPPair[4][4] = 0;
    mBPPair[4][5] = 0;
    mBPPair[4][6] = 0;
    mBPPair[4][7] = 6;

    mBPPair[5][0] = 0;
    mBPPair[5][1] = 0;
    mBPPair[5][2] = 0;
    mBPPair[5][3] = 0;
    mBPPair[5][4] = 0;
    mBPPair[5][5] = 0;
    mBPPair[5][6] = 2;
    mBPPair[5][7] = 0;

    mBPPair[6][0] = 0;
    mBPPair[6][1] = 0;
    mBPPair[6][2] = 0;
    mBPPair[6][3] = 0;
    mBPPair[6][4] = 0;
    mBPPair[6][5] = 1;
    mBPPair[6][6] = 0;
    mBPPair[6][7] = 0;

    mBPPair[7][0] = 0;
    mBPPair[7][1] = 6;
    mBPPair[7][2] = 0;
    mBPPair[7][3] = 0;
    mBPPair[7][4] = 5;
    mBPPair[7][5] = 0;
    mBPPair[7][6] = 0;
    mBPPair[7][7] = 0;

    mRtype[0] = 0;
    mRtype[1] = 2;
    mRtype[2] = 1;
    mRtype[3] = 4;
    mRtype[4] = 3;
    mRtype[5] = 6;
    mRtype[6] = 5;
    mRtype[7] = 7;

}

int VR16PairMat::encode_char(char pChr, int pEnergySet) {
    /* return numerical representation of base used e.g. in pair[][] */
    int code;

    if (pEnergySet > 0) {
        code = (int) (pChr - 'A') + 1;
    } else {
        code = 0;
        unsigned i = 0;
        for (std::string::iterator it = mLawAndOrder.begin(); it != mLawAndOrder.end(); ++it) {
            if (pChr == (char) (*it)) {
                code = i;
                break;
            }
            ++i;
        }

        if (code > 4) {
            --code; /* make T and U equivalent */
        }
    }

    return code;
}

int VR16PairMat::make_pair_matrix(int pEnergySet, std::string &pNonStandards, bool pNoGU) {
    if (pEnergySet == 0) {
        for (unsigned i = 0; i < 5; ++i) {
            mAlias[i] = (int) i;
        }
        mAlias[5] = 3;     /* X <-> G */
        mAlias[6] = 2;     /* K <-> C */
        mAlias[7] = 0;     /* I <-> default base '@' */
        for (int i = 0; i < NBASES; ++i) {
            for (int j = 0; j < NBASES; ++j) {
                mPair[i][j] = mBPPair[i][j];
            }
        }
        if (pNoGU) {
            mPair[3][4] = 0;
            mPair[4][3] = 0;
        }
        if (pNonStandards.size() != 0) { /* allow nonstandard bp's */
            for (unsigned i = 0; i < pNonStandards.size(); i += 2) {
                mPair[encode_char(pNonStandards[i], pEnergySet)][encode_char(pNonStandards[i + 1], pEnergySet)] = 7;
            }
        }
        for (int i = 0; i < NBASES; ++i) {
            for (int j = 0; j < NBASES; ++j) {
                mRtype[mPair[i][j]] = mPair[j][i];
            }
        }
    } else {
        for (int i = 0; i <= MAXALPHA; ++i) {
            for (int j = 0; j <= MAXALPHA; ++j) {
                mPair[i][j] = 0;
            }
        }
        if (pEnergySet == 1) {
            for (int i = 1; i < MAXALPHA; ++i) {
                mAlias[i] = 3; /* A <-> G */
                ++i;
                mAlias[i] = 2; /* B <-> C */
            }
            for (int i = 1; i < MAXALPHA; ++i) {
                mPair[i][i + 1] = 2; /* AB <-> GC */
                ++i;
                mPair[i][i - 1] = 1; /* BA <-> CG */
            }
        } else if (pEnergySet == 2) {
            for (int i = 1; i < MAXALPHA; ++i) {
                mAlias[i] = 1; /* A <-> A*/
                ++i;
                mAlias[i] = 4; /* B <-> U */
            }
            for (int i = 1; i < MAXALPHA; ++i) {
                mPair[i][i + 1] = 5; /* AB <-> AU */
                ++i;
                mPair[i][i - 1] = 6; /* BA <-> UA */
            }
        } else if (pEnergySet == 3) {
            for (unsigned i = 1; i < MAXALPHA - 2; ++i) {
                mAlias[i] = 3; /* A <-> G */
                ++i;
                mAlias[i] = 2; /* B <-> C */
                ++i;
                mAlias[i] = 1; /* C <-> A */
                ++i;
                mAlias[i] = 4; /* D <-> U */
            }
            for (unsigned i = 1; i < MAXALPHA - 2; ++i) {
                mPair[i][i + 1] = 2; /* AB <-> GC */
                ++i;
                mPair[i][i - 1] = 1; /* BA <-> CG */
                ++i;
                mPair[i][i + 1] = 5; /* CD <-> AU */
                ++i;
                mPair[i][i - 1] = 6; /* DC <-> UA */
            }
        } else {
            std::cerr << "Error: Wrong energy_set." << std::endl;
            return 1;
        }

        for (int i = 0; i <= MAXALPHA; ++i) {
            for (int j = 0; j <= MAXALPHA; ++j) {
                mRtype[mPair[i][j]] = mPair[j][i];
            }
        }
    }

    return 0;
}

void VR16PairMat::init_heap() {

    mBPPair = mikan::create_2d_array<int>(NBASES, NBASES);

    mAlias = mikan::create_1d_array<int>(MAXALPHA + 1);
    mPair = mikan::create_2d_array<int>(MAXALPHA + 1, MAXALPHA + 1);
    mRtype = mikan::create_1d_array<int>(8);
}

void VR16PairMat::free_heap() {
    mikan::delete_2d_array<int>(mBPPair, NBASES);

    mikan::delete_1d_array<int>(mAlias);
    mikan::delete_2d_array<int>(mPair, MAXALPHA + 1);
    mikan::delete_1d_array<int>(mRtype);
}

} // namespace vr16
