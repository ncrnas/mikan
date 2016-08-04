#include <mikan/lib/rnahybrid2_mfe/RNAhybrid/include/hybrid_backtrace.hpp>


namespace rh2 {

int RH2BackTrace::move_next_bt_string()
{
    ++mCurPos;
    return mCurPos;
}

void RH2BackTrace::build_bt_string(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen)
{
    RH2FuncID func_id;
    int i1, j1, i2, j2;

    back_hybrid(pTargetPos, pTargetLen, pQueryPos, pQueryLen, mBTStrng[mCurPos]);
    func_id = mBTStrng[mCurPos].get_func_id();
    while(func_id != FUNC_NONE)
    {
        i1 = mBTStrng[mCurPos].get_i1();
        j1 = mBTStrng[mCurPos].get_j1();
        i2 = mBTStrng[mCurPos].get_i2();
        j2 = mBTStrng[mCurPos].get_j2();

        ++mCurPos;
        if (func_id == FUNC_UNPAIRED_LEFT)
        {
            back_unpaired_left_bot(i1, j1, i2, j2, mBTStrng[mCurPos]);
        }
        else if (func_id == FUNC_CLOSED)
        {
            back_closed(i1, j1, i2, j2, mBTStrng[mCurPos]);
        }

        func_id = mBTStrng[mCurPos].get_func_id();
    }

    mEndPos = mCurPos + 1;
    mCurPos = 0;
}

void RH2BackTrace::back_hybrid(
        int pTargetPos,
        int pTargetLen,
        int pQueryPos,
        int pQueryLen,
        RH2Signature &pSig)
{
    float mfe, min_mfe;

    /* ---------------------------------- start of --------------------------------- */
    /* --------------------- v1 = nil <<< tt(uregion, uregion) --------------------- */

    RH2Signature v1(SIGID_NULL, FUNC_NONE, 65000);
    if (((pTargetLen - pTargetPos) >= 0) && ((pQueryLen - pQueryPos) >= 0))
    {
        mfe = 0.0;
        v1.set_values(SIGID_Nil, FUNC_NONE, mfe,
                pTargetPos, pTargetLen, pQueryPos, pQueryLen);
        /* No iteration neccessary! */
    }

    /* ---------------------------------- finished --------------------------------- */

    /* -------------------------- v2 = p unpaired_left_top ------------------------- */
    /* +------------------------------------------------------------------------------------ */
    /* Nonterminal unpaired_left_top is implemented as a tabulated                           */
    /* function which yields atomar results. Since we are in list context,                   */
    /* we need to wrap the result of unpaired_left_top into a single list element.           */
    /* +------------------------------------------------------------------------------------ */

    /*    Find best unpaired_left_bot without backtracking recursion to avoid */
    /*      stack overflow (hand made): */

    int k;
    int best_k = -1;
    RH2Signature v2(SIGID_NULL, FUNC_NONE, 65000);
    if (((pTargetLen - pTargetPos) >= 1) && ((pQueryLen - pQueryPos) >= 1)
            && ((pQueryPos < mQueryHelixStart) || (pQueryPos >= mQueryHelixEnd)))
    {
        min_mfe = 65000;
        for (k = pTargetPos; k < pTargetLen; k++)
        {
            mfe = tbl.get_unpaired_left_bot(k, pQueryPos);
            if ((mProcessedX.count(k+1) != 1)
                    && (mfe < min_mfe) && (mfe > mMinMfe))
            {
                best_k = k;
                min_mfe = mfe;
            }
        }
        if (best_k != -1)
        {
            v2.set_values(SIGID__NTID, FUNC_UNPAIRED_LEFT, min_mfe,
                    best_k, pTargetLen, pQueryPos, pQueryLen);
        }

    }

    /* ------------------------------- v3 = p closed ------------------------------- */
    /* +---------------------------------------------------------------------------- */
    /* Nonterminal closed is implemented as a tabulated                              */
    /* function which yields atomar results. Since we are in list context,           */
    /* we need to wrap the result of closed into a single list element.              */
    /* +---------------------------------------------------------------------------- */

    RH2Signature v3(SIGID_NULL, FUNC_NONE, 65000);
    if ((mProcessedX.count(1) != 1))
    {
        mfe = tbl.get_closed(pTargetPos, pQueryPos);
        v3.set_values(SIGID__NTID, FUNC_CLOSED, mfe,
                pTargetPos, pTargetLen, pQueryPos, pQueryLen);
    }

    /* ---------------------------- v4 = minimum(v2, v3) --------------------------- */

    bool v1_min = false;
    bool v2_min = false;
    bool v3_min = false;

    if (v2.get_mfe() < v3.get_mfe())
    {
        gx = best_k + 1;
        v2_min = true;
        min_mfe = v2.get_mfe();
    }
    else
    {
        gx = 0 + 1;
        v3_min = true;
        min_mfe = v3.get_mfe();
    }

    /*    v4 = v2.alg_mfe < v3.alg_mfe ? v2 : v3; */
    /* ---------------------------- v5 = minimum(v1, v4) --------------------------- */
    if (v1.get_mfe() < min_mfe)
    {
        gx = 0;
        v1_min = true;
        v2_min = false;
        v3_min = false;
    }

    /*    v5 = v1.alg_mfe < v4.alg_mfe ? v1 : v4; */
    /* ------------------------- build candidate structures ------------------------ */

    if (v1_min)
    {
        pSig = v1;
    }
    else if (v2_min)
    {
        pSig = v2;
    }
    else if (v3_min)
    {
        pSig = v3;
    }

    mMinMfe = pSig.get_mfe();

    //	std::cout << pSig.get_op_id() << "-"  << pSig.get_func_id() << ": " << pSig.get_i1() << ", " << pSig.get_j1() << ", ";
    //	std::cout << pSig.get_i2() << ", " << pSig.get_j2() << ", " << pSig.get_mfe() << std::endl;

}

/* table calculation for production unpaired_left_top                               */
/* -------------------------------------------------------------------------------- */

/* removed after hand made back_hybrid optimisation !*/

/* table calculation for production unpaired_left_bot                               */
/* -------------------------------------------------------------------------------- */

void RH2BackTrace::back_unpaired_left_bot(
        int pTargetPos,
        int pTargetLen,
        int pQueryPos,
        int pQueryLen,
        RH2Signature &pSig)
{
    float mfe;

    if (((pTargetLen - pTargetPos) < 1) && ((pQueryLen - pQueryPos) < 1))
    {
        return;
    }

    /* ---------------------------------- start of --------------------------------- */
    /* ---------- v1 = ulb <<< (tt(empty, lbase)) ~~~ p unpaired_left_bot ---------- */

    RH2Signature v1(SIGID_NULL, FUNC_NONE, 65000);
    if (((pQueryLen - pQueryPos) >= 2)
            && ((pQueryPos < mQueryHelixStart - 1) || (pQueryPos >= mQueryHelixEnd)))
    {
        mfe = tbl.get_unpaired_left_bot(pTargetPos, pQueryPos + 1);
        v1.set_values(SIGID_Ulb, FUNC_UNPAIRED_LEFT, mfe,
                pTargetPos, pTargetLen, pQueryPos + 1, pQueryLen,
                pTargetPos, pQueryPos + 1);
    }

    /* ---------- v1 = ulb <<< (tt(empty, lbase)) ~~~ p unpaired_left_bot ---------- */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /* ---------------- v2 = eds <<< (tt(lbase, lbase)) ~~~ p closed --------------- */

    RH2Signature v2(SIGID_NULL, FUNC_NONE, 65000);
    if (((pTargetLen - pTargetPos) >= 2) && ((pQueryLen - pQueryPos) >= 2)
            && en.is_pair(inpx(pTargetPos + 2), inpy(pQueryPos + 2))
            && ((pQueryPos < mQueryHelixStart) || (pQueryPos >= mQueryHelixEnd)))
    {
        mfe = (tbl.get_closed(pTargetPos + 1, pQueryPos + 1)
                + en.dl_energy((pTargetPos + 1) + 1, (pQueryPos + 1) + 1))
				        + en.dr_energy((pTargetPos + 1) + 1, (pQueryPos + 1) + 1);
        v2.set_values(SIGID_Eds, FUNC_CLOSED, mfe,
                pTargetPos + 1, pTargetLen, pQueryPos + 1, pQueryLen,
                pTargetPos + 1, pQueryPos + 1);
        /* No iteration neccessary! */
    }

    /* ---------------- v2 = eds <<< (tt(lbase, lbase)) ~~~ p closed --------------- */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /* ---------------- v3 = edt <<< (tt(lbase, empty)) ~~~ p closed --------------- */
    RH2Signature v3(SIGID_NULL, FUNC_NONE, 65000);
    if (((pTargetLen - pTargetPos) >= 2)
            && en.is_pair(inpx(pTargetPos + 2), inpy(pQueryPos + 1)))
    {
        mfe = tbl.get_closed(pTargetPos + 1, pQueryPos) + en.dl_energy((pTargetPos + 1) + 1, (pQueryPos) + 1);
        v3.set_values(SIGID_Edt, FUNC_CLOSED, mfe,
                pTargetPos + 1, pTargetLen, pQueryPos, pQueryLen,
                pTargetPos + 1, pQueryPos);
        /* No iteration neccessary! */
    }

    /* ---------------- v3 = edt <<< (tt(lbase, empty)) ~~~ p closed --------------- */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /* ---------------- v4 = edb <<< (tt(empty, lbase)) ~~~ p closed --------------- */
    RH2Signature v4(SIGID_NULL, FUNC_NONE, 65000);
    if (((pQueryLen - pQueryPos) >= 2)
            && en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 2))
            && ((pQueryPos < mQueryHelixStart) || (pQueryPos >= mQueryHelixEnd)))
    {
        mfe = tbl.get_closed(pTargetPos, pQueryPos + 1) + en.dr_energy((pTargetPos) + 1, (pQueryPos + 1) + 1);
        v4.set_values(SIGID_Edb, FUNC_CLOSED, mfe,
                pTargetPos, pTargetLen, pQueryPos + 1, pQueryLen,
                pTargetPos, pQueryPos + 1);
        /* No iteration neccessary! */
    }

    /* ---------------- v4 = edb <<< (tt(empty, lbase)) ~~~ p closed --------------- */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------- v5 = minimum(v3, v4) --------------------------- */
    bool v1_min = false;
    bool v2_min = false;
    bool v3_min = false;
    bool v4_min = false;
    float min_mfe = 65000;

    if (v3.get_mfe() < v4.get_mfe())
    {
        v3_min = true;
        min_mfe = v3.get_mfe();
    }
    else
    {
        v4_min = true;
        min_mfe = v4.get_mfe();
    }

    /*    v5 = v3.alg_mfe < v4.alg_mfe ? v3 : v4; */
    /* ---------------------------- v6 = minimum(v2, v5) --------------------------- */
    if (v2.get_mfe() < min_mfe)
    {
        v2_min = true;
        v3_min = false;
        v4_min = false;
        min_mfe = v2.get_mfe();
    }

    /*    v6 = v2.alg_mfe < v5.alg_mfe ? v2 : v5; */
    /* ---------------------------- v7 = minimum(v1, v6) --------------------------- */
    if (v1.get_mfe() < min_mfe)
    {
        v1_min = true;
        v2_min = false;
        v3_min = false;
        v4_min = false;
        min_mfe = v1.get_mfe();
    }

    /*    v7 = v1.alg_mfe < v6.alg_mfe ? v1 : v6; */
    /* ------------------------- build candidate structures ------------------------ */

    if (v1_min)
    {
        pSig = v1;
    }
    else if (v2_min)
    {
        pSig = v2;
    }
    else if (v3_min)
    {
        pSig = v3;
    }
    else if (v4_min)
    {
        pSig = v4;
    }

    //	std::cout << pSig.get_op_id() << "-"  << pSig.get_func_id() << ": " << pSig.get_i1() << ", " << pSig.get_j1() << ", ";
    //	std::cout << pSig.get_i2() << ", " << pSig.get_j2() << ", " << pSig.get_mfe() << std::endl;
}

/* table calculation for production closed                                          */
/* -------------------------------------------------------------------------------- */

void RH2BackTrace::back_closed(
        int pTargetPos,
        int pTargetLen,
        int pQueryPos,
        int pQueryLen,
        RH2Signature &pSig)
{
    float mfe, min_mfe;
    //str1 v1, v2, v3, v4, v5, v6, v7, v7b, v7c, v7d, v7e, v7f, v7g, v8, v9, v10, v11, v12;
    int k;
    int k2;
    int k3;
    int k4;

    if (((pTargetLen - pTargetPos) < 1) && ((pQueryLen - pQueryPos) < 1))
    {
        return;
    }

    /* ---------------------------------- start of --------------------------------- */
    /* - v1 = sr <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ p closed - */
    RH2Signature v1(SIGID_NULL, FUNC_NONE, 65000);
    if (((pTargetLen - pTargetPos) >= 2) && ((pQueryLen - pQueryPos) >= 2))
    {
        if (en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1)))
        {
            mfe = en.sr_energy(pTargetPos + 1, pQueryPos + 1)
					        + tbl.get_closed(pTargetPos + 1, pQueryPos + 1);
            v1.set_values(SIGID_Sr, FUNC_CLOSED, mfe,
                    pTargetPos + 1, pTargetLen, pQueryPos + 1, pQueryLen,
                    pTargetPos + 1, pQueryPos + 1);
            /* No iteration neccessary! */
        }
    }

    /* - v1 = sr <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ p closed - */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /*  v3 = bt <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(region, empty) `with` (sizeTT 1 15 0 0)) ~~~ p closed  */

    min_mfe = 65000;
    RH2Signature v3(SIGID_NULL, FUNC_NONE, 65000);
    if (((pTargetLen - pTargetPos) >= 3) && ((pQueryLen - pQueryPos) >= 2)
            && en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1))
            && ((pQueryPos < mQueryHelixStart) || (pQueryPos >= mQueryHelixEnd - 1)))
    {
        for (k = pTargetPos + 2; k <= std::min(pTargetPos + mBloopUpperLimit + 1, pTargetLen - 1); k++)
        {
            if (inpx(k) == CHAR_X)
            {
                break;
            }

            mfe = (tbl.get_closed(k, pQueryPos + 1)
                    + en.bl_stacking((k) - (pTargetPos + 1), 0, pTargetPos + 1, pQueryPos + 1))
					        + en.bl_ent((k) - (pTargetPos + 1));

            /* No iteration neccessary! */
            /* ------------------------- v3 = minimum(v2, v3) ------------------------ */
            if (mfe < min_mfe)
            {
                v3.set_values(SIGID_Bt, FUNC_CLOSED, mfe,
                        k, pTargetLen, pQueryPos + 1, pQueryLen,
                        pTargetPos + 1, pQueryPos + 1, pTargetPos + 1, k, pQueryPos + 1);
                min_mfe = mfe;
            }
        }
    }

    /*  v3 = bt <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(region, empty) `with` (sizeTT 1 15 0 0)) ~~~ p closed  */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /*  v5 = bb <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(empty, region) `with` (sizeTT 0 0 1 15)) ~~~ p closed  */

    min_mfe = 65000;
    RH2Signature v5(SIGID_NULL, FUNC_NONE, 65000);
    if (((pTargetLen - pTargetPos) >= 2) && ((pQueryLen - pQueryPos) >= 3)
            && en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1)))
    {
        for (k2 = pQueryPos + 2; k2 <= std::min(pQueryPos + mBloopUpperLimit + 1, pQueryLen - 1); k2++)
        {
            if ((k2 > mQueryHelixStart) && (k2 <= mQueryHelixEnd))
            {
                break;
            }
            mfe = (tbl.get_closed(pTargetPos + 1, k2)
                    + en.bl_stacking(0, (k2) - (pQueryPos + 1), pTargetPos + 1, pQueryPos + 1))
					        + en.bl_ent((k2) - (pQueryPos + 1));

            /* No iteration neccessary! */
            /* ------------------------- v5 = minimum(v4, v5) ------------------------ */
            /* v5 = v4.alg_mfe < v5.alg_mfe ? v4 : v5; */

            if (mfe < min_mfe)
            {
                v5.set_values(SIGID_Bb, FUNC_CLOSED, mfe,
                        pTargetPos + 1, pTargetLen, k2, pQueryLen,
                        pTargetPos + 1, pQueryPos + 1, pTargetPos + 1, pQueryPos + 1, k2);
                min_mfe = mfe;
            }

        }
    }

    /*  v5 = bb <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(empty, region) `with` (sizeTT 0 0 1 15)) ~~~ p closed  */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /*  v7 = il <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(region, region) `with` (sizeTT 1 15 1 15)) ~~~ p closed  */
    min_mfe = 65000;
    RH2Signature v7(SIGID_NULL, FUNC_NONE, 65000);
    if (((pTargetLen - pTargetPos) >= 3) && ((pQueryLen - pQueryPos) >= 3)
            && en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1)))
    {
        /* special internal loops: */
        for (k3 = pTargetPos + 2; k3 <= std::min(pTargetPos+std::min(3,mIloopUpperLimit+1), pTargetLen - 1); k3++)
        {
            if (inpx(k3) == CHAR_X)
            {
                break;
            }

            for (k4 = pQueryPos + 2; k4 <= std::min(pQueryPos+std::min(3,mIloopUpperLimit+1), pQueryLen - 1); k4++)
            {
                if ((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd - 1))
                {
                    break;
                }
                if (en.is_pair(inpx(k3 + 1), inpy(k4 + 1)))
                {
                    mfe = en.do_il_special(pTargetPos + 1,
                            pQueryPos + 1,
                            pTargetPos + 1,
                            k3,
                            pQueryPos + 1,
                            k4,
                            tbl.get_closed(k3, k4));
                    /* No iteration neccessary! */
                    if (mfe < min_mfe)
                    {
                        v7.set_values(SIGID_Il, FUNC_CLOSED, mfe,
                                k3, pTargetLen, k4, pQueryLen,
                                pTargetPos + 1, pQueryPos + 1, pTargetPos + 1, k3, pQueryPos + 1, k4);
                        min_mfe = mfe;
                    }
                }
            }
        }

        min_mfe = 65000;
        RH2Signature v7b(SIGID_NULL, FUNC_NONE, 65000);
        RH2Signature v7g(SIGID_NULL, FUNC_NONE, 65000);
        if ((inpx(k3) != CHAR_X))  /*  && !((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd-1))) { */
        {
            /* normal internal loops: */
            for (k3 = pTargetPos + 2; k3 <= std::min(pTargetPos + 3, pTargetLen - 1); k3++)
            {
                if (inpx(k3) == CHAR_X)
                {
                    break;
                }

                for (k4 = pQueryPos + 4; k4 <= std::min(pQueryPos + mIloopUpperLimit + 1, pQueryLen - 1); k4++)
                {
                    if ((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd - 1))
                    {
                        break;
                    }
                    if (en.is_pair(inpx(k3 + 1), inpy(k4 + 1)))
                    {
                        mfe = en.do_il(pTargetPos + 1,
                                pQueryPos + 1,
                                pTargetPos + 1,
                                k3,
                                pQueryPos + 1,
                                k4,
                                tbl.get_closed(k3, k4));
                        /* No iteration neccessary! */
                        if (mfe < min_mfe)
                        {
                            v7b.set_values(SIGID_Il, FUNC_CLOSED, mfe,
                                    k3, pTargetLen, k4, pQueryLen,
                                    pTargetPos + 1, pQueryPos + 1, pTargetPos + 1, k3, pQueryPos + 1, k4);
                            min_mfe = mfe;
                        }
                    }
                }
            }
            if (min_mfe < 65000)
            {
                v7b.set_mfe(min_mfe + en.il_stack_open(pTargetPos + 1, pQueryPos + 1));
            }

            min_mfe = 65000;
            RH2Signature v7c(SIGID_NULL, FUNC_NONE, 65000);
            if ((inpx(k3) != CHAR_X))  /*  && !((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd-1))) { */
            {
                for (k3 = pTargetPos + 4; k3 <= std::min(pTargetPos + mIloopUpperLimit + 1, pTargetLen - 1); k3++)
                {
                    if (inpx(k3) == CHAR_X)
                    {
                        break;
                    }

                    for (k4 = pQueryPos + 2; k4 <= std::min(pQueryPos + 3, pQueryLen - 1); k4++)
                    {
                        if ((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd - 1))
                        {
                            break;
                        }
                        if (en.is_pair(inpx(k3 + 1), inpy(k4 + 1)))
                        {
                            mfe = en.do_il(pTargetPos + 1,
                                    pQueryPos + 1,
                                    pTargetPos + 1,
                                    k3,
                                    pQueryPos + 1,
                                    k4,
                                    tbl.get_closed(k3, k4));
                            /* No iteration neccessary! */
                            if (mfe < min_mfe)
                            {
                                v7c.set_values(SIGID_Il, FUNC_CLOSED, mfe,
                                        k3, pTargetLen, k4, pQueryLen,
                                        pTargetPos + 1, pQueryPos + 1, pTargetPos + 1, k3, pQueryPos + 1, k4);
                                min_mfe = mfe;
                            }
                        }
                    }
                }
                if (min_mfe < 65000)
                {
                    v7c.set_mfe(min_mfe + en.il_stack_open(pTargetPos + 1, pQueryPos + 1));
                }
            }

            min_mfe = 65000;
            RH2Signature v7d(SIGID_NULL, FUNC_NONE, 65000);
            if ((inpx(k3) != CHAR_X)) /*  && !((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd-1))) { */
            {
                /* normal internal loops: */
                for (k3 = pTargetPos + 4; k3 <= std::min(pTargetPos + mIloopUpperLimit + 1, pTargetLen - 1); k3++)
                {
                    if (inpx(k3) == CHAR_X)
                    {
                        break;
                    }

                    for (k4 = pQueryPos + 4; k4 <= std::min(pQueryPos + mIloopUpperLimit + 1, pQueryLen - 1); k4++)
                    {
                        if ((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd - 1))
                        {
                            break;
                        }
                        if (en.is_pair(inpx(k3 + 1), inpy(k4 + 1)))
                        {
                            mfe = en.do_il(pTargetPos + 1,
                                    pQueryPos + 1,
                                    pTargetPos + 1,
                                    k3,
                                    pQueryPos + 1,
                                    k4,
                                    tbl.get_closed(k3, k4));
                            /* No iteration neccessary! */
                            if (mfe < min_mfe)
                            {
                                v7d.set_values(SIGID_Il, FUNC_CLOSED, mfe,
                                        k3, pTargetLen, k4, pQueryLen,
                                        pTargetPos + 1, pQueryPos + 1, pTargetPos + 1, k3, pQueryPos + 1, k4);
                                min_mfe = mfe;
                            }
                        }
                    }
                }
                if (min_mfe < 65000)
                {
                    v7d.set_mfe(min_mfe + en.il_stack_open(pTargetPos + 1, pQueryPos + 1));
                }

            }

            bool v7b_min = false;
            bool v7c_min = false;
            bool v7d_min = false;

            if (v7b.get_mfe() < v7c.get_mfe())
            {
                v7b_min = true;
                min_mfe = v7b.get_mfe();
            }
            else
            {
                v7c_min = true;
                min_mfe = v7c.get_mfe();
            }

            if (v7d.get_mfe() < min_mfe)
            {
                v7b_min = false;
                v7c_min = false;
                v7d_min = true;
                min_mfe = v7d.get_mfe();
            }

            if (v7g.get_mfe() > min_mfe)
            {
                if (v7b_min)
                {
                    v7g = v7b;
                }
                else if (v7c_min)
                {
                    v7g = v7c;
                }
                else if (v7d_min)
                {
                    v7g = v7d;
                }
            }

        }

        if (v7g.get_mfe() < v7.get_mfe())
        {
            v7 = v7g;
        }

    }

    /*  v7 = il <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(region, region) `with` (sizeTT 1 15 1 15)) ~~~ p closed  */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /*  v8 = el <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(uregion, uregion))  */
    min_mfe = 65000;
    RH2Signature v8(SIGID_NULL, FUNC_NONE, 65000);
    if ((pTargetLen == pTargetPos + 1 || inpx(pTargetPos + 2) != 'X')
            && ((pQueryPos >= mQueryHelixEnd - 1) || (mQueryHelixEnd > pQueryLen)))
    {
        if (en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1)))
        {
            mfe = ((((pTargetLen) - (pTargetPos + 1)) > 0) ? en.dli_energy(pTargetPos + 1, pQueryPos + 1) : 0)
				          + ((((pQueryLen) - (pQueryPos + 1)) > 0) ? en.dri_energy(pTargetPos + 1, pQueryPos + 1) : 0);
            v8.set_values(SIGID_El, FUNC_NONE, mfe,
                    0, 0, 0, 0,
                    pTargetPos + 1, pQueryPos + 1, pTargetPos + 1, pTargetLen, pQueryPos + 1, pQueryLen);
            /* No iteration neccessary! */
        }

    }

    /*  v8 = el <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(uregion, uregion))  */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------- v9 = minimum(v7, v8) --------------------------- */
    bool v8_min = false;
    bool v7_min = false;
    bool v5_min = false;
    bool v3_min = false;
    bool v1_min = false;

    if (v7.get_mfe() < v8.get_mfe())
    {
        v7_min = true;
        min_mfe = v7.get_mfe();
    }
    else
    {
        v8_min = true;
        min_mfe = v8.get_mfe();
    }
    /*    v9 = v7.alg_mfe < v8.alg_mfe ? v7 : v8; */

    /* --------------------------- v10 = minimum(v5, v9) --------------------------- */
    if (v5.get_mfe() < min_mfe)
    {
        v5_min = true;
        v7_min = false;
        v8_min = false;
        min_mfe = v5.get_mfe();
    }

    /*    v10 = v5.alg_mfe < v9.alg_mfe ? v5 : v9; */

    /* --------------------------- v11 = minimum(v3, v10) -------------------------- */
    if (v3.get_mfe() < min_mfe)
    {
        v3_min = true;
        v5_min = false;
        v7_min = false;
        v8_min = false;
        min_mfe = v3.get_mfe();
    }

    /*    v11 = v3.alg_mfe < v10.alg_mfe ? v3 : v10; */

    /* --------------------------- v12 = minimum(v1, v11) -------------------------- */
    if (v1.get_mfe() < min_mfe)
    {
        v1_min = true;
        v3_min = false;
        v5_min = false;
        v7_min = false;
        v8_min = false;
        min_mfe = v1.get_mfe();
    }

    /*    v12 = v1.alg_mfe < v11.alg_mfe ? v1 : v11; */

    /* ------------------------- build candidate structures ------------------------ */

    if (v1_min)
    {
        pSig = v1;
    }
    else if (v3_min)
    {
        pSig = v3;
    }
    else if (v5_min)
    {
        pSig = v5;
    }
    else if (v7_min)
    {
        pSig = v7;
    }
    else if (v8_min)
    {
        pSig = v8;
    }
    //	std::cout << pSig.get_op_id() << "-"  << pSig.get_func_id() << ": " << pSig.get_i1() << ", " << pSig.get_j1() << ", ";
    //	std::cout << pSig.get_i2() << ", " << pSig.get_j2() << ", " << pSig.get_mfe() << std::endl;
}

} // namespace rh2
