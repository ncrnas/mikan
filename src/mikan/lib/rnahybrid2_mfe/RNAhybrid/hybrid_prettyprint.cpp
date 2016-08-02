#include <mikan/lib/rnahybrid2_mfe/RNAhybrid/include/hybrid_backtrace.hpp>
#include <mikan/lib/rnahybrid2_mfe/RNAhybrid/include/hybrid_prettyprint.hpp>

namespace rh2 {

void RH2PrettyPrinter::sx(std::string::iterator& t, int i)
{
    switch(inpx(i))
    {
    case CHAR_A:
        *t = 'A';
        break;
    case CHAR_C:
        *t = 'C';
        break;
    case CHAR_G:
        *t = 'G';
        break;
    case CHAR_U:
        *t = 'U';
        break;
    default:
        *t = 'N';
    }
    *(++t) = 0;
}

void RH2PrettyPrinter::shift(std::string::iterator& t, int u, int v)
{
    if ((v - u) == 0)
    {
        *t = ' ';
        *(++t) = ' ';
        *(++t) = 0;
    }
    else
    {
        *t = ' ';
        *(++t) = 0;
    }
}

void RH2PrettyPrinter::one_space(std::string::iterator& t)
{
    *t = ' ';
    *(++t) = 0;
}

void RH2PrettyPrinter::ssx(std::string::iterator& t, int i, int j)
{
    for (int k = i + 1; k <= j; k++)
    {
        sx(t, k);
    }
}

void RH2PrettyPrinter::sy(std::string::iterator& t, int i)
{
    switch(inpy(i))
    {
    case CHAR_A:
        *t = 'A';
        break;
    case CHAR_C:
        *t = 'C';
        break;
    case CHAR_G:
        *t = 'G';
        break;
    case CHAR_U:
        *t = 'U';
        break;
    default:
        *t = 'N';
    }
    *(++t) = 0;
}

void RH2PrettyPrinter::ssy(std::string::iterator& t, int i, int j)
{
    for (int k = i + 1; k <= j; k++)
    {
        sy(t, k);
    }
}

void RH2PrettyPrinter::blanks(std::string::iterator& t, int i, int j)
{
    int k, l;

    l = j - i;
    for (k = 0; k <= l - 1; ++k, ++t)
    {
        *t = ' ';
    }

    *t = 0;
}

void RH2PrettyPrinter::blanks2(std::string::iterator& t, int i, int j)
{
    int k, l;

    l = std::max(j - i - 1, 0);
    for (k = 0; k <= l - 1; ++k, ++t)
    {
        *t = ' ';
    }

    *t = 0;
}

void RH2PrettyPrinter::asym1(std::string::iterator& t, int l, int r, int u, int v)
{
    int k, ln;

    ln = std::max((v - u) - (r - l), 0);
    for (k = 0; k <= ln - 1; ++k, ++t)
    {
        *t = ' ';
    }

    *t = 0;
}

void RH2PrettyPrinter::asym2(std::string::iterator& t, int l, int r, int u, int v)
{
    int k, ln;

    ln = std::max((r - l) - (v - u), 0);
    for (k = 0; k <= ln - 1; ++k, ++t)
    {
        *t = ' ';
    }

    *t = 0;
}

void RH2PrettyPrinter::pp_str_Hybrid()
{
    RH2SigID opType;
    int curPos = 0;
    int endPos = bt.get_end_pos();
    int a1, a2, a3, a4, a5, a6;

    while (curPos < endPos)
    {
        opType = bt.get_op_type(curPos);
        a1 = bt.get_a1(curPos);
        a2 = bt.get_a2(curPos);
        a3 = bt.get_a3(curPos);
        a4 = bt.get_a4(curPos);
        a5 = bt.get_a5(curPos);
        a6 = bt.get_a6(curPos);

        if (opType == SIGID_Ulb)
        {
            one_space(r1);
            one_space(r2);
            one_space(r3);
            sy(r4, a2);
        }
        else if (opType == SIGID_Eds)
        {
            sx(r1, a1);
            one_space(r2);
            one_space(r3);
            sy(r4, a2);

            mCountToA1 = a1;
        }
        else if (opType == SIGID_Edt)
        {
            sx(r1, a1);
            one_space(r2);
            one_space(r3);
            one_space(r4);

            mCountToA1 = a1;
        }
        else if (opType == SIGID_Edb)
        {
            one_space(r1);
            one_space(r2);
            one_space(r3);
            sy(r4, a2);
        }
        else if (opType == SIGID_Sr)
        {
            one_space(r1);
            sx(r2, a1);
            sy(r3, a2);
            one_space(r4);

            mCountToA1 = a1;
        }
        else if (opType == SIGID_Bt)
        {
            one_space(r1);
            ssx(r1, a3, a4);
            sx(r2, a1);
            blanks(r2, a3, a4);
            sy(r3, a2);
            blanks(r3, a3, a4);
            one_space(r4);
            blanks(r4, a3, a4);

            mCountToA1 = a4;
        }
        else if (opType == SIGID_Bb)
        {
            one_space(r1);
            blanks(r1, a4, a5);
            sx(r2, a1);
            blanks(r2, a4, a5);
            sy(r3, a2);
            blanks(r3, a4, a5);
            one_space(r4);
            ssy(r4, a4, a5);

            mCountToA1 = a1;
        }
        else if (opType == SIGID_Il)
        {
            one_space(r1);
            ssx(r1, a3, a4);
            asym1(r1, a3, a4, a5, a6);

            sx(r2, a1);
            blanks(r2, a3, a4);
            asym1(r2, a3, a4, a5, a6);

            sy(r3, a2);
            blanks(r3, a5, a6);
            asym2(r3, a3, a4, a5, a6);

            one_space(r4);
            ssy(r4, a5, a6);
            asym2(r4, a3, a4, a5, a6);

            mCountToA1 = a4;
        }
        else if (opType == SIGID_El)
        {
            if ((a4 - a3) == 0)
            {
                one_space(r1);
                blanks(r1, a5, a6);

                sx(r2, a3);
                blanks(r2, a5, a6);

                sy(r3, a5);
                blanks(r3, a5, a6);

                one_space(r4);
                ssy(r4, a5, a6);

                mCountToA1 = a3;
            }
            else
            {
                one_space(r1);
                sx(r1, a3 + 1);
                blanks2(r1, a5, a6);

                sx(r2, a3);
                one_space(r2);
                blanks2(r2, a5, a6);

                sy(r3, a5);
                one_space(r3);
                blanks2(r3, a5, a6);

                shift(r4, a5, a6);
                ssy(r4, a5, a6);

                mCountToA1 = a3 + 1;
            }
        }
        else if (opType == SIGID_Nil)
        {
            *r1 = '\0';
            *r2 = '\0';
            *r3 = '\0';
            *r4 = '\0';
        }

        ++curPos;
    }
}

int RH2PrettyPrinter::get_hit_length()
{
    int k;
    int hit_length = 0;

    std::string str_t1(t1.c_str());
    int t_len = (int)str_t1.length();

    for (k = 0; k < t_len && t1[k] == ' ' && t2[k] == ' '; ++k) {}
    if (k < t_len)
    {
        t1[k] = ' '; /* hide left dangle */
    }

    for (k = t_len - 1; k >= 0 && t1[k] == ' ' && t2[k] == ' '; --k) {}
    if (k >= 0)
    {
        t1[k] = ' '; /* hide right dangle */
    }

    for (k = 0; k < t_len; k++)
    {
        if (t1[k] != ' ')
        {
            hit_length++;
        }
    }

    for (k = 0; k < t_len; k++)
    {
        if (t2[k] != ' ')
        {
            hit_length++;
        }
    }

    return hit_length;
}

} // namespace rh2
