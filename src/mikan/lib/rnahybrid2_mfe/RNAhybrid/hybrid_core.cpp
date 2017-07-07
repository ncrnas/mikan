/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

/* Modified by Takaya Saito 2014.01 */

/* compiled by the ADP compiler, version 0.8.342                                    */
/* source file: rnahybrid/HybridMFE3_simple_3tab.lhs                                */
/* command:                                                                         */
/* adpcompile -c rnahybrid/HybridMFE3_simple_3tab.lhs -al mfe enum -bt -bts -o rnahybrid/hybridBack4.c */
/* -------------------------------------------------------------------------------- */

#include "hybrid_core.hpp"


namespace rh2 {

/* main dynamic programming loop                                                    */
/* -------------------------------------------------------------------------------- */

void RH2WorkSpace::init_workspace(int pTargetLen, int pQueryLen) {
    tbl.table_alloc(pTargetLen, pQueryLen);
    bt.alloc_string(pTargetLen, pQueryLen);
    pp.alloc_string(pTargetLen, pQueryLen);
}

void RH2WorkSpace::set_target_seq(std::vector<char> &seq, int qTargetLen) {
    target_seq = seq;
    en.set_target_seq(&target_seq);
    tbl.set_target_seq(&target_seq, qTargetLen);
    bt.set_target_helix_start(tbl.get_target_helix_start());
    bt.set_target_helix_end(tbl.get_target_helix_end());
    pp.set_target_seq(&target_seq);
    mTargetLen = qTargetLen;
}

void RH2WorkSpace::set_query_seq(std::vector<char> &seq, int qQueryLen) {
    query_seq = seq;
    en.set_query_seq(&query_seq);
    tbl.set_query_seq(&query_seq, qQueryLen);
    bt.set_query_helix_start(tbl.get_query_helix_start());
    bt.set_query_helix_end(tbl.get_query_helix_end());
    pp.set_query_seq(&query_seq);
    mQueryLen = qQueryLen;
}

void RH2WorkSpace::create_ret_val(RH2RetValues &pRetVal, float pV2, int pSearchCount) {
    pRetVal.mMfe = pV2;
    pRetVal.mTargetPos = bt.get_gx();
    pRetVal.mTargetPos0 = bt.get_gx() - 1;
    pRetVal.mTargetSub = pp.t1.c_str();
    pRetVal.mTarget = pp.t2.c_str();
    pRetVal.mQeuery = pp.t3.c_str();
    pRetVal.mQeuarySub = pp.t4.c_str();
    pRetVal.mLenA1 = pp.get_count_to_A1();
    pRetVal.mSearchCount = pSearchCount;
    pRetVal.mEffective = false;
}

void RH2WorkSpace::copy_ret_val(RH2RetValues &pRetVal1, RH2RetValues &pRetVal2) {
    pRetVal1.mMfe = pRetVal2.mMfe;
    pRetVal1.mTargetPos = pRetVal2.mTargetPos;
    pRetVal1.mTargetPos0 = pRetVal2.mTargetPos0;
    pRetVal1.mTargetSub = pRetVal2.mTargetSub;
    pRetVal1.mTarget = pRetVal2.mTarget;
    pRetVal1.mQeuery = pRetVal2.mQeuery;
    pRetVal1.mQeuarySub = pRetVal2.mQeuarySub;
    pRetVal1.mLenA1 = pRetVal2.mLenA1;
    pRetVal1.mSearchCount = pRetVal2.mSearchCount;
    pRetVal1.mEffective = pRetVal2.mEffective;
}

void RH2WorkSpace::reset_bt_work_space() {
    bt.reset();
    pp.reset();
}

void RH2WorkSpace::reset_work_space() {
    bt.clear_processed_x();
}

void RH2WorkSpace::mainloop(RH2RetValues &pRetVal) {
    float v2;
    float min_v2;
    bool first_iter = true;
    bool loopContinue = true;
    int hit_length, hit_start;
    int iterate_counter = 0;
    std::string str_bt;
    RH2RetValues tmpRet;

    for (int targetPos = mTargetLen; targetPos >= 0; --targetPos) {
        for (int queryPos = mQueryLen; queryPos >= 0; --queryPos) {
            tbl.calc_unpaired_left_bot(targetPos, mTargetLen, queryPos, mQueryLen);
            tbl.calc_closed(targetPos, mTargetLen, queryPos, mQueryLen);
            tbl.calc_unpaired_left_top(targetPos, mTargetLen, queryPos, mQueryLen);
        }
    }

    v2 = tbl.calc_hybrid(0, mTargetLen, 0, mQueryLen);
    if (v2 == 0.0) {
        loopContinue = false;
    }

    pRetVal.mEffective = false;
    min_v2 = v2;

    while (loopContinue) {
        bt.build_bt_string(0, mTargetLen, 0, mQueryLen);
        v2 = bt.get_min_mfe();
        if (v2 == 0.0) {
            loopContinue = false;
        }

        pp.pp_str_Hybrid();
        str_bt = pp.t1.c_str();
        if (str_bt.length() == 0) {
            loopContinue = false;
        }

        create_ret_val(tmpRet, v2, iterate_counter);
        if (pp.get_count_to_A1() >= mTargetLen - 1) {
            loopContinue = false;
        }

        hit_length = pp.get_hit_length();
        hit_start = bt.get_gx() + (pp.t1[0] != ' ');
        if (loopContinue) {
            bt.add_processed_x();
        } else {
            tmpRet.mHitLen = hit_length;
            tmpRet.mHitStart = hit_start;
        }

        if (v2 < min_v2 || first_iter) {
            copy_ret_val(pRetVal, tmpRet);
            pRetVal.mEffective = true;
            min_v2 = v2;
        }

        first_iter = false;
        reset_bt_work_space();
        ++iterate_counter;
    }

    reset_work_space();

}

} // namespace rh2
