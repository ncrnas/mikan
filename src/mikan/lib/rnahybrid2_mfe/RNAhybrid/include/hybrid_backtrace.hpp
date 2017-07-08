#ifndef HYBRID_BACKTRACE_HPP_
#define HYBRID_BACKTRACE_HPP_

#include <vector>
#include <set>
#include "hybrid_energy.hpp"
#include "hybrid_table.hpp"

namespace rh2 {

const char CHAR_A = 0;
const char CHAR_C = 1;
const char CHAR_G = 2;
const char CHAR_U = 3;
const char CHAR_N = 4;
const char CHAR_X = 5;

enum RH2SigID {
    SIGID_NULL = 0,
    SIGID__NTID,
    SIGID_Ult,
    SIGID_Ulb,
    SIGID_Eds,
    SIGID_Edt,
    SIGID_Edb,
    SIGID_Sr,
    SIGID_Bt,
    SIGID_Bb,
    SIGID_Il,
    SIGID_El,
    SIGID_Nil
};

enum RH2FuncID {
    FUNC_NONE = 0,
    FUNC_UNPAIRED_LEFT,
    FUNC_CLOSED
};

class RH2Signature {
public:
    RH2Signature() :
            mOpType(SIGID_NULL), mFuncType(FUNC_NONE), mAlgMFE(0),
            mI1(0), mJ1(0), mI2(0), mJ2(0),
            mA1(0), mA2(0), mA3(0), mA4(0), mA5(0), mA6(0),
            mCalcedVal(0) {}

    RH2Signature(RH2SigID pOpType, RH2FuncID pFuncType,
                 float pAlgMFE, int pI1 = 0, int pJ1 = 0, int pI2 = 0, int pJ2 = 0,
                 int pA1 = 0, int pA2 = 0, int pA3 = 0, int pA4 = 0, int pA5 = 0, int pA6 = 0) :
            mOpType(pOpType), mFuncType(pFuncType), mAlgMFE(pAlgMFE),
            mI1(pI1), mJ1(pJ1), mI2(pI2), mJ2(pJ2),
            mA1(pA1), mA2(pA2), mA3(pA3), mA4(pA4), mA5(pA5), mA6(pA6),
            mCalcedVal(0) {}

    RH2Signature(const RH2Signature &pSig) {
        mOpType = pSig.mOpType;
        mFuncType = pSig.mFuncType;
        mAlgMFE = pSig.mAlgMFE;
        mI1 = pSig.mI1;
        mJ1 = pSig.mJ1;
        mI2 = pSig.mI2;
        mJ2 = pSig.mJ2;
        mA1 = pSig.mA1;
        mA2 = pSig.mA2;
        mA3 = pSig.mA3;
        mA4 = pSig.mA4;
        mA5 = pSig.mA5;
        mA6 = pSig.mA6;
        mCalcedVal = pSig.mCalcedVal;
    }

    void set_values(RH2SigID pOpType, RH2FuncID pFuncType,
                    float pAlgMFE, int pI1 = 0, int pJ1 = 0, int pI2 = 0, int pJ2 = 0,
                    int pA1 = 0, int pA2 = 0, int pA3 = 0, int pA4 = 0, int pA5 = 0, int pA6 = 0) {
        mOpType = pOpType;
        mFuncType = pFuncType;
        mAlgMFE = pAlgMFE;
        mI1 = pI1;
        mJ1 = pJ1;
        mI2 = pI2;
        mJ2 = pJ2;
        mA1 = pA1;
        mA2 = pA2;
        mA3 = pA3;
        mA4 = pA4;
        mA5 = pA5;
        mA6 = pA6;
    }

    float get_mfe() { return mAlgMFE; }

    void set_mfe(float pMFE) { mAlgMFE = pMFE; }

    RH2FuncID get_func_id() { return mFuncType; }

    RH2SigID get_op_id() { return mOpType; }

    int get_i1() { return mI1; }

    int get_j1() { return mJ1; }

    int get_i2() { return mI2; }

    int get_j2() { return mJ2; }

    int get_a1() { return mA1; }

    int get_a2() { return mA2; }

    int get_a3() { return mA3; }

    int get_a4() { return mA4; }

    int get_a5() { return mA5; }

    int get_a6() { return mA6; }

private:
    RH2SigID mOpType;
    RH2FuncID mFuncType;
    float mAlgMFE;
    int mI1, mJ1, mI2, mJ2;
    int mA1, mA2, mA3, mA4, mA5, mA6;
    int mCalcedVal;

};

class RH2BackTrace {
public:
    void alloc_string(int max_target_len, int max_query_len) {
        int size = 2 * std::max(max_target_len, max_query_len);
        mBTStrng.resize(size, RH2Signature());
    }

    RH2BackTrace(RH2EnergyFunc &pEn, RH2Table &pTbl, int max_target_len, int max_query_len) :
            mCurPos(0), mEndPos(0), gx(0), mMinMfe(-65000),
            en(pEn), tbl(pTbl),
            mIloopUpperLimit(15),
            mBloopUpperLimit(15),
            mQueryHelixStart(0),
            mQueryHelixEnd(0),
            mTargetHelixStart(0),
            mTargetHelixEnd(0) {
        alloc_string(max_target_len, max_query_len);
    }

    void set_iloop_upper_limit(int l) { mIloopUpperLimit = l; }

    void set_bloop_upper_limit(int l) { mBloopUpperLimit = l; }

    void set_query_helix_start(int pStart) { mQueryHelixStart = pStart; }

    void set_query_helix_end(int pEnd) { mQueryHelixEnd = pEnd; }

    void set_target_helix_start(int pStart) { mTargetHelixStart = pStart; }

    void set_target_helix_end(int pEnd) { mTargetHelixEnd = pEnd; }

    void back_hybrid(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen, RH2Signature &pSig);

    void back_unpaired_left_bot(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen, RH2Signature &pSig);

    void back_closed(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen, RH2Signature &pSig);

    void reset() {
        mCurPos = 0;
        mEndPos = 0;
        gx = 0;
    }

    void build_bt_string(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen);

    int move_next_bt_string();

    RH2SigID get_op_type(int pIdx) { return mBTStrng[pIdx].get_op_id(); }

    int get_a1(int pIdx) { return mBTStrng[pIdx].get_a1(); }

    int get_a2(int pIdx) { return mBTStrng[pIdx].get_a2(); }

    int get_a3(int pIdx) { return mBTStrng[pIdx].get_a3(); }

    int get_a4(int pIdx) { return mBTStrng[pIdx].get_a4(); }

    int get_a5(int pIdx) { return mBTStrng[pIdx].get_a5(); }

    int get_a6(int pIdx) { return mBTStrng[pIdx].get_a6(); }

    int get_end_pos() { return mEndPos; }

    int get_gx() { return gx; }

    void add_processed_x() { mProcessedX.insert(gx); };

    void clear_processed_x() {
        mProcessedX.clear();
        mMinMfe = -65000;
    }

    float get_min_mfe() { return mMinMfe; }

private:
    std::vector<RH2Signature> mBTStrng;
    int mCurPos, mEndPos;

    int gx;
    std::set<int> mProcessedX;
    float mMinMfe;

    RH2EnergyFunc &en;
    RH2Table &tbl;

    int mIloopUpperLimit;
    int mBloopUpperLimit;

    int mQueryHelixStart;
    int mQueryHelixEnd;
    int mTargetHelixStart;
    int mTargetHelixEnd;

    int inpx(int i) { return tbl.inpx(i); }

    int inpy(int i) { return tbl.inpy(i); }

};

} // namespace rh2


#endif /* HYBRID_BACKTRACE_HPP_ */
