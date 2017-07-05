#ifndef HYBRID_CORE_HPP_
#define HYBRID_CORE_HPP_

#include <hybrid_backtrace.hpp>
#include <hybrid_energy.hpp>
#include <hybrid_prettyprint.hpp>
#include <hybrid_table.hpp>
#include <string>
#include <vector>

namespace rh2 {

class RH2RetValues {
public:
    RH2RetValues() :
            mMfe(0), mTargetPos(0), mTargetPos0(0),
            mTargetSub(""), mTarget(""),
            mQeuery(""), mQeuarySub(""),
            mHitStart(0), mHitLen(0), mLenA1(0), mSearchCount(0),
            mEffective(false) {}

    RH2RetValues(float pMfe, int pTargetPos, int pTargetPos0,
                 std::string pTargetSub, std::string pTarget,
                 std::string pQeuery, std::string pQeuarySub,
                 int pHitStart, int pHitLen, int pLenA1, int pSearchCount,
                 bool pEffective) :
            mMfe(pMfe), mTargetPos(pTargetPos), mTargetPos0(pTargetPos0),
            mTargetSub(pTargetSub), mTarget(pTarget),
            mQeuery(pQeuery), mQeuarySub(pQeuarySub),
            mHitStart(pHitStart), mHitLen(pHitLen), mLenA1(pLenA1),
            mSearchCount(pSearchCount),
            mEffective(pEffective) {}

    RH2RetValues(const RH2RetValues &pRetVal) :
            mMfe(pRetVal.mMfe), mTargetPos(pRetVal.mTargetPos), mTargetPos0(pRetVal.mTargetPos0),
            mTargetSub(pRetVal.mTargetSub), mTarget(pRetVal.mTarget),
            mQeuery(pRetVal.mQeuery), mQeuarySub(pRetVal.mQeuarySub),
            mHitStart(pRetVal.mHitStart), mHitLen(pRetVal.mHitLen),
            mLenA1(pRetVal.mLenA1), mSearchCount(pRetVal.mSearchCount),
            mEffective(pRetVal.mEffective) {}

    float mMfe;
    int mTargetPos;
    int mTargetPos0;
    std::string mTargetSub;
    std::string mTarget;
    std::string mQeuery;
    std::string mQeuarySub;

    int mHitStart;
    int mHitLen;
    int mLenA1;
    int mSearchCount;

    bool mEffective;
};

class RH2WorkSpace {
public:
    RH2WorkSpace(int max_target_len, int max_query_len, std::string &seed_def) :
            en(), tbl(en, max_target_len, max_query_len, seed_def),
            bt(en, tbl, max_target_len, max_query_len),
            pp(bt, max_target_len, max_query_len),
            mTargetMaxLen(max_target_len), mQueryMaxLen(max_query_len),
            mTargetLen(0), mQueryLen(0) {}

    void set_target_seq(std::vector<char> &seq, int mTargetLen);

    void set_query_seq(std::vector<char> &seq, int mQueryLen);

    void set_iloop_upper_limit(int l) { tbl.set_iloop_upper_limit(l); }

    void set_bloop_upper_limit(int l) { tbl.set_bloop_upper_limit(l); }

    void mainloop(RH2RetValues &pRetVal);

    void init_workspace(int pTargetLen, int pQueryLen);

private:
    RH2EnergyFunc en;
    RH2Table tbl;
    RH2BackTrace bt;
    RH2PrettyPrinter pp;
    std::vector<char> target_seq;
    std::vector<char> query_seq;
    int mTargetMaxLen, mQueryMaxLen;
    int mTargetLen, mQueryLen;

    bool is_effective_site();

    void create_ret_val(RH2RetValues &pRetVal, float pV2, int pSearchCount);

    void copy_ret_val(RH2RetValues &pRetVal1, RH2RetValues &pRetVal2);

    void reset_bt_work_space();

    void reset_work_space();

    void mask_hit(int idx) { target_seq[idx] = CHAR_X; }
};

} // namespace rh2

#endif /* HYBRID_CORE_HPP_ */
