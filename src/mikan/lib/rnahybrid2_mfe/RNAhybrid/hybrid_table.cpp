#include <mikan/lib/rnahybrid2_mfe/RNAhybrid/include/hybrid_energy.hpp>
#include <mikan/lib/rnahybrid2_mfe/RNAhybrid/include/hybrid_table.hpp>

namespace rh2 {

void RH2Table::set_target_seq(std::vector<char> *seq, int)
{
    target_seq = seq;
}
void RH2Table::set_query_seq(std::vector<char> *seq, int)
{
    query_seq = seq;
    //	set_query_seed_pos(qQueryLen);
}

void RH2Table::set_target_seed_pos(int qTargetLen)
{
    mTargetHelixStart = qTargetLen - seed_end - 1;
    mTargetHelixEnd = qTargetLen - seed_start;
}

void RH2Table::set_query_seed_pos(int qQueryLen)
{
    mQueryHelixStart = qQueryLen - seed_end;
    mQueryHelixEnd = qQueryLen - seed_start + 1;
}

/* table memory allocation                                                          */
/* -------------------------------------------------------------------------------- */

void RH2Table::table_alloc(int pTargetMaxLen, int pQueryMaxLen)
{
    tbl_unpaired_left_bot.resize((unsigned)pTargetMaxLen);
    for (int i = 0; i < pTargetMaxLen; ++i)
    {
        tbl_unpaired_left_bot[i].resize((unsigned)pQueryMaxLen, 0);
    }

    tbl_closed.resize((unsigned)pTargetMaxLen);
    for (int i = 0; i < pTargetMaxLen; ++i)
    {
        tbl_closed[i].resize((unsigned)pQueryMaxLen, 0);
    }

    tbl_unpaired_left_top.resize((unsigned)pTargetMaxLen);
    for (int i = 0; i < pTargetMaxLen; ++i)
    {
        tbl_unpaired_left_top[i].resize((unsigned)pQueryMaxLen, 0);
    }

    if (seed_def[0] == '7')
    {
        seed_start = 2;
        seed_end = 8;
    }
    else if (seed_def[0] == '6')
    {
        seed_start = 2;
        seed_end = 7;
    }
    mQueryHelixStart = pQueryMaxLen - seed_end;
    mQueryHelixEnd = pQueryMaxLen - seed_start + 1;
}

/* table calculations                                                               */
/* -------------------------------------------------------------------------------- */

/* table calculation for production hybrid                                          */
/* -------------------------------------------------------------------------------- */

float RH2Table::calc_hybrid(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen)
{
    float v1, v2, v3, v4, v5;

    /* ---------------------------------- start of --------------------------------- */
    /* --------------------- v1 = nil <<< tt(uregion, uregion) --------------------- */
    if (((pTargetLen - pTargetPos) >= 0) && ((pQueryLen - pQueryPos) >= 0))
    {
        v1 = 0;
        /* No iteration neccessary! */
    }
    else
    {
        v1 = 65000;
    }
    /* --------------------- v1 = nil <<< tt(uregion, uregion) --------------------- */
    /* ---------------------------------- finished --------------------------------- */

    /* -------------------------- v2 = p unpaired_left_top ------------------------- */
    if (((pTargetLen - pTargetPos) >= 1) && ((pQueryLen - pQueryPos) >= 1))
    {
        v2 = tbl_unpaired_left_top[pTargetPos][pQueryPos];
    }
    else
    {
        v2 = 65000;
    }
    /* ------------------------------- v3 = p closed ------------------------------- */
    if (((pTargetLen - pTargetPos) >= 1) && ((pQueryLen - pQueryPos) >= 1))
    {
        v3 = tbl_closed[pTargetPos][pQueryPos];
    }
    else
    {
        v3 = 65000;
    }
    v4 = v2 < v3 ? v2 : v3;
    v5 = v1 < v4 ? v1 : v4;

    return (v5);
}

/* table calculation for production unpaired_left_top                               */
/* -------------------------------------------------------------------------------- */

void RH2Table::calc_unpaired_left_top(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen)
{
    float v1, v2, v3;

    if (((pTargetLen - pTargetPos) < 1) || ((pQueryLen - pQueryPos) < 1))
    {
        return;
    }

    /* ---------------------------------- start of --------------------------------- */
    /* ---------- v1 = ult <<< (tt(lbase, empty)) ~~~ p unpaired_left_top ---------- */
    if (((pTargetLen - pTargetPos) >= 2))
    {
        v1 = tbl_unpaired_left_top[pTargetPos + 1][pQueryPos];
        /* No iteration neccessary! */
    }
    else
    {
        v1 = 65000;
    }
    /* ---------- v1 = ult <<< (tt(lbase, empty)) ~~~ p unpaired_left_top ---------- */
    /* ---------------------------------- finished --------------------------------- */

    /* -------------------------- v2 = p unpaired_left_bot ------------------------- */
    if ((pQueryPos < mQueryHelixStart) || (pQueryPos >= mQueryHelixEnd))
    {
        v2 = tbl_unpaired_left_bot[pTargetPos][pQueryPos];
    }
    else
    {
        v2 = 65000;
    }
    v3 = v1 < v2 ? v1 : v2;
    /* ------------------------- assign table entry result ------------------------- */

    tbl_unpaired_left_top[pTargetPos][pQueryPos] = v3;

}

/* table calculation for production unpaired_left_bot                               */
/* -------------------------------------------------------------------------------- */

void RH2Table::calc_unpaired_left_bot(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen)
{
    float v1, v2, v3, v4, v5, v6, v7;

    if (((pTargetLen - pTargetPos) < 1) || ((pQueryLen - pQueryPos) < 1))
    {
        return;
    }

    /* ---------------------------------- start of --------------------------------- */
    /* ---------- v1 = ulb <<< (tt(empty, lbase)) ~~~ p unpaired_left_bot ---------- */
    if (((pQueryLen - pQueryPos) >= 2)
          && ((pQueryPos < mQueryHelixStart - 1) || (pQueryPos >= mQueryHelixEnd)))
    {
        v1 = tbl_unpaired_left_bot[pTargetPos][pQueryPos + 1];
        /* No iteration neccessary! */
    }
    else
    {
        v1 = 65000;
    }
    /* ---------- v1 = ulb <<< (tt(empty, lbase)) ~~~ p unpaired_left_bot ---------- */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /* ---------------- v2 = eds <<< (tt(lbase, lbase)) ~~~ p closed --------------- */
    if (((pTargetLen - pTargetPos) >= 2) && ((pQueryLen - pQueryPos) >= 2)
          && en.is_pair(inpx(pTargetPos + 2), inpy(pQueryPos + 2))
          && ((pQueryPos < mQueryHelixStart) || (pQueryPos >= mQueryHelixEnd)))
    {
        v2 = (tbl_closed[pTargetPos + 1][pQueryPos + 1]
             + en.dl_energy((pTargetPos + 1) + 1, (pQueryPos + 1) + 1))
             + en.dr_energy((pTargetPos + 1) + 1, (pQueryPos + 1) + 1);
        /* No iteration neccessary! */
    }
    else
    {
        v2 = 65000;
    }
    /* ---------------- v2 = eds <<< (tt(lbase, lbase)) ~~~ p closed --------------- */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /* ---------------- v3 = edt <<< (tt(lbase, empty)) ~~~ p closed --------------- */
    if (((pTargetLen - pTargetPos) >= 2)
          && en.is_pair(inpx(pTargetPos + 2), inpy(pQueryPos + 1)))
    {
        v3 = tbl_closed[pTargetPos + 1][pQueryPos]
             + en.dl_energy((pTargetPos + 1) + 1, (pQueryPos) + 1);
        /* No iteration neccessary! */
    }
    else
    {
        v3 = 65000;
    }
    /* ---------------- v3 = edt <<< (tt(lbase, empty)) ~~~ p closed --------------- */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /* ---------------- v4 = edb <<< (tt(empty, lbase)) ~~~ p closed --------------- */

    if (((pQueryLen - pQueryPos) >= 2)
          && en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 2))
          && ((pQueryPos < mQueryHelixStart) || (pQueryPos >= mQueryHelixEnd)))
    {
        v4 = tbl_closed[pTargetPos][pQueryPos + 1]
             + en.dr_energy((pTargetPos) + 1, (pQueryPos + 1) + 1);
        /* No iteration neccessary! */
    }
    else
    {
        v4 = 65000;
    }
    /* ---------------- v4 = edb <<< (tt(empty, lbase)) ~~~ p closed --------------- */
    /* ---------------------------------- finished --------------------------------- */

    v5 = v3 < v4 ? v3 : v4;
    v6 = v2 < v5 ? v2 : v5;
    v7 = v1 < v6 ? v1 : v6;
    /* ------------------------- assign table entry result ------------------------- */

    tbl_unpaired_left_bot[pTargetPos][pQueryPos] = v7;

}

/* table calculation for production closed                                          */
/* -------------------------------------------------------------------------------- */
void RH2Table::calc_closed(int pTargetPos, int pTargetLen, int pQueryPos, int pQueryLen)
{
    float v1, v2, v3, v4, v5, v6, v7, v7b, v7c, v7d, v7e, v7f, v7g, v8, v9, v10, v11, v12;
    int k;
    int k2;
    int k3;
    int k4;
    int k_min, k_min2;

    if (((pTargetLen - pTargetPos) < 1) || ((pQueryLen - pQueryPos) < 1))
    {
        return;
    }

    v7d = 65000;
    /* ---------------------------------- start of --------------------------------- */
    /* - v1 = sr <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ p closed - */
    if (((pTargetLen - pTargetPos) >= 2) && ((pQueryLen - pQueryPos) >= 2))
    {
        if (en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1)))
        {
            v1 = en.sr_energy(pTargetPos + 1, pQueryPos + 1)
               + tbl_closed[pTargetPos + 1][pQueryPos + 1];
            /* No iteration neccessary! */
        }
        else
        {
            v1 = 65000;
        }
    }
    else
    {
        v1 = 65000;
    }
    /* - v1 = sr <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ p closed - */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /*  v3 = bt <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(region, empty) `with` (sizeTT 1 15 0 0)) ~~~ p closed  */
    if (((pTargetLen - pTargetPos) >= 3) && ((pQueryLen - pQueryPos) >= 2)
          && en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1))
          && ((pQueryPos < mQueryHelixStart) || (pQueryPos >= mQueryHelixEnd - 1)))
    {
        v3 = 65000;
        k_min = std::min(pTargetPos + bloop_upper_limit + 1, pTargetLen - 1);
        for (k = pTargetPos + 2; k <= k_min; ++k)
        {
            if (inpx(k) == TX)
            {
                break;
            }

            v2 = (tbl_closed[k][pQueryPos + 1]
                 + en.bl_stacking((k) - (pTargetPos + 1), 0, pTargetPos + 1, pQueryPos + 1))
                 + en.bl_ent((k) - (pTargetPos + 1));
            /* No iteration neccessary! */
            v3 = v2 < v3 ? v2 : v3;
        }
    }
    else
    {
        v3 = 65000;
    }
    /*  v3 = bt <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(region, empty) `with` (sizeTT 1 15 0 0)) ~~~ p closed  */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /*  v5 = bb <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(empty, region) `with` (sizeTT 0 0 1 15)) ~~~ p closed  */
    if (((pTargetLen - pTargetPos) >= 2) && ((pQueryLen - pQueryPos) >= 3)
          && en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1)))
    {
        v5 = 65000;
        k_min = std::min(pQueryPos + bloop_upper_limit + 1, pQueryLen - 1);
        for (k2 = pQueryPos + 2; k2 <= k_min; ++k2)
        {
            if ((k2 > mQueryHelixStart) && (k2 <= mQueryHelixEnd))
            {
                break;
            }
            v4 = (tbl_closed[pTargetPos + 1][k2]
                 + en.bl_stacking(0, (k2) - (pQueryPos + 1), pTargetPos + 1, pQueryPos + 1))
                 + en.bl_ent((k2) - (pQueryPos + 1));
            /* No iteration neccessary! */
            v5 = v4 < v5 ? v4 : v5;

        }
    }
    else
    {
        v5 = 65000;
    }
    /*  v5 = bb <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(empty, region) `with` (sizeTT 0 0 1 15)) ~~~ p closed  */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /*  v7 = il <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(region, region) `with` (sizeTT 1 15 1 15)) ~~~ p closed  */

    if (((pTargetLen - pTargetPos) >= 3) && ((pQueryLen - pQueryPos) >= 3)
          && en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1)))
    {
        v7 = 65000;
        /* special internal loops: */
        k_min = std::min(pTargetPos+std::min(3,iloop_upper_limit+1), pTargetLen - 1);
        for (k3 = pTargetPos + 2; k3 <= k_min; ++k3)
        {
            if (inpx(k3) == TX)
            {
                break;
            }

            k_min2 = std::min(pQueryPos+std::min(3,iloop_upper_limit+1), pQueryLen - 1);
            for (k4 = pQueryPos + 2; k4 <= k_min2; ++k4)
            {
                if ((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd - 1))
                {
                    break;
                }
                if (en.is_pair(inpx(k3 + 1), inpy(k4 + 1)))
                {
                    v6 = en.do_il_special(pTargetPos + 1,
                            pQueryPos + 1,
                            pTargetPos + 1,
                            k3,
                            pQueryPos + 1,
                            k4,
                            tbl_closed[k3][k4]);
                    /* No iteration neccessary! */
                    v7 = v6 < v7 ? v6 : v7;

                }
            }
        }
        v7g = 65000;
        v7b = 65000;
        if ((inpx(k3) != TX))  /*  && !((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd-1))) { */
        {
            /* normal internal loops: */
            k_min =  std::min(pTargetPos + 3, pTargetLen - 1);
            for (k3 = pTargetPos + 2; k3 <= k_min; ++k3)
            {
                if (inpx(k3) == TX)
                {
                    break;
                }

                k_min2 = std::min(pQueryPos + iloop_upper_limit + 1, pQueryLen - 1);
                for (k4 = pQueryPos + 4; k4 <= k_min2; ++k4)
                {
                    if ((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd - 1))
                    {
                        break;
                    }
                    if (en.is_pair(inpx(k3 + 1), inpy(k4 + 1)))
                    {
                        v6 = en.do_il(pTargetPos + 1,
                                pQueryPos + 1,
                                pTargetPos + 1,
                                k3,
                                pQueryPos + 1,
                                k4,
                                tbl_closed[k3][k4]);
                        /* No iteration neccessary! */
                        v7b = v6 < v7b ? v6 : v7b;

                    }
                }
            }
            if (v7b < 65000)
            {
                v7b += en.il_stack_open(pTargetPos + 1, pQueryPos + 1);
            }
            v7c = 65000;
            /* normal internal loops: */
            if ((inpx(k3) != TX))  /*  && !((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd-1))) { */
            {
                k_min = std::min(pTargetPos + iloop_upper_limit + 1, pTargetLen - 1);
                for (k3 = pTargetPos + 4; k3 <= k_min; ++k3)
                {
                    if (inpx(k3) == TX)
                    {
                        break;
                    }

                    k_min2 = std::min(pQueryPos + 3, pQueryLen - 1);
                    for (k4 = pQueryPos + 2; k4 <= k_min2; ++k4)
                    {
                        if ((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd - 1))
                        {
                            break;
                        }
                        if (en.is_pair(inpx(k3 + 1), inpy(k4 + 1)))
                        {
                            v6 = en.do_il(pTargetPos + 1,
                                    pQueryPos + 1,
                                    pTargetPos + 1,
                                    k3,
                                    pQueryPos + 1,
                                    k4,
                                    tbl_closed[k3][k4]);
                            /* No iteration neccessary! */
                            v7c = v6 < v7c ? v6 : v7c;

                        }
                    }
                }
                if (v7c < 65000)
                {
                    v7c += en.il_stack_open(pTargetPos + 1, pQueryPos + 1);
                }
                v7d = 65000;
                /* normal internal loops: */
                if ((inpx(k3) != TX)) /*  && !((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd-1))) { */
                {
                    k_min = std::min(pTargetPos + iloop_upper_limit + 1, pTargetLen - 1);
                    for (k3 = pTargetPos + 4; k3 <= k_min; ++k3)
                    {
                        if (inpx(k3) == TX)
                        {
                            break;
                        }

                        k_min2 = std::min(pQueryPos + iloop_upper_limit + 1, pQueryLen - 1);
                        for (k4 = pQueryPos + 4; k4 <= k_min2; ++k4)
                        {
                            if ((k4 > mQueryHelixStart) && (pQueryPos < mQueryHelixEnd - 1))
                            {
                                break;
                            }
                            if (en.is_pair(inpx(k3 + 1), inpy(k4 + 1)))
                            {
                                v6 = en.do_il(pTargetPos + 1,
                                        pQueryPos + 1,
                                        pTargetPos + 1,
                                        k3,
                                        pQueryPos + 1,
                                        k4,
                                        tbl_closed[k3][k4]);
                                /* No iteration neccessary! */
                                v7d = v6 < v7d ? v6 : v7d;

                            }
                        }
                    }
                    if (v7d < 65000)
                    {
                        v7d += en.il_stack_open(pTargetPos + 1, pQueryPos + 1);
                    }
                }

            }
            v7e = v7b < v7c ? v7b : v7c;
            v7f = v7d < v7e ? v7d : v7e;
            v7g = v7g < v7f ? v7g : v7f;
        }
        v7 = v7g < v7 ? v7g : v7;
    }
    else
    {
        v7 = 65000;
    }
    /*  v7 = il <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(region, region) `with` (sizeTT 1 15 1 15)) ~~~ p closed  */
    /* ---------------------------------- finished --------------------------------- */

    /* ---------------------------------- start of --------------------------------- */
    /*  v8 = el <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(uregion, uregion))  */
    if ((pTargetLen == pTargetPos + 1 || inpx(pTargetPos + 2) != 'X')
         && ((pQueryPos >= mQueryHelixEnd - 1) || (mQueryHelixEnd > pQueryLen)))
    {
        if (en.is_pair(inpx(pTargetPos + 1), inpy(pQueryPos + 1)))
        {
            v8 = ((((pTargetLen) - (pTargetPos + 1)) > 0) ? en.dli_energy(pTargetPos + 1, pQueryPos + 1) : 0)
                 + ((((pQueryLen) - (pQueryPos + 1)) > 0) ? en.dri_energy(pTargetPos + 1, pQueryPos + 1) : 0);
            /* No iteration neccessary! */
        }
        else
        {
            v8 = 65000;
        }
    }
    else
    {
        v8 = 65000;
    }
    /*  v8 = el <<< (tt(lbase, lbase) `with` (pairingTTcross compl_ij)) ~~~ (tt(uregion, uregion))  */
    /* ---------------------------------- finished --------------------------------- */

    v9 = v7 < v8 ? v7 : v8;
    v10 = v5 < v9 ? v5 : v9;
    v11 = v3 < v10 ? v3 : v10;
    v12 = v1 < v11 ? v1 : v11;
    /* ------------------------- assign table entry result ------------------------- */

    tbl_closed[pTargetPos][pQueryPos] = v12;

}

} // namespace rh2
