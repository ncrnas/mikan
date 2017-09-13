#ifndef DP_ALIGN_HPP_
#define DP_ALIGN_HPP_

#include <string>
#include <deque>
#include <iostream>

namespace mr3dp {

class MR3DPAlign {

public:

    // Constructor
    MR3DPAlign () {}

    // Create alignment from path
    void create_align(std::string &pQSeq, std::string &pDSeq, std::deque<int> &pPathI, std::deque<int> &pPathJ) {
        mQAlignSeq = "";
        mAlignBar = "";
        mDAlignSeq = "";

        int prevI = pPathI[0];
        int prevJ = pPathJ[0];

        for (unsigned k = 1; k < pPathI.size(); k++) {
            if (prevI == pPathI[k]) {
                mQAlignSeq += "-";
                mDAlignSeq += pDSeq[pPathJ[k] - 1];
                mAlignBar += " ";
            } else if (prevJ == pPathJ[k]) {
                mQAlignSeq += pQSeq[pPathI[k] - 1];
                mDAlignSeq += "-";
                mAlignBar += " ";
            } else {
                mQAlignSeq += pQSeq[pPathI[k] - 1];
                mDAlignSeq += pDSeq[pPathJ[k] - 1];
                if ((pQSeq[pPathI[k] - 1] == 'G' && pDSeq[pPathJ[k] - 1] == 'U')
                        || (pQSeq[pPathI[k] - 1] == 'U' && pDSeq[pPathJ[k] - 1] == 'G') ) {
                    mAlignBar += ":";
                } else {
                    mAlignBar += "|";
                }

            }
            prevI = pPathI[k];
            prevJ = pPathJ[k];
        }

        for (unsigned i = prevI; i < pQSeq.size(); ++i) {

            mQAlignSeq += pQSeq[i];
            if (i < pDSeq.size()) {
                mDAlignSeq += pDSeq[i];
            } else {
                mDAlignSeq += " ";
            }
            mAlignBar += " ";

        }

    }

    // Get query sequence of alignment
    std::string &get_q_align() {
        return mQAlignSeq;
    }

    // Get database sequence of alignment
    std::string &get_d_align() {
        return mDAlignSeq;
    }

    // Print alignment
    void print_align() {
        std::cout << mQAlignSeq << std::endl;
        std::cout << mAlignBar << std::endl;
        std::cout << mDAlignSeq << std::endl;
    }


private:
    // Alignment
    std::string mQAlignSeq;
    std::string mAlignBar;
    std::string mDAlignSeq;
};

} // namespace mr3dp

#endif /* DP_ALIGN_HPP_ */