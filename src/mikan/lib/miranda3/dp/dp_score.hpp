#ifndef DP_SCORE_HPP_
#define DP_SCORE_HPP_

namespace mr3dp {

// The user defined data table.
//      A   C   G   U
//  A  -3  -3  -3   5
//  C  -3  -3   5  -3
//  G  -3   5  -3   1
//  U   5  -3   1  -3

class MR3DPScore {
public:
    int gap_open() {
        return 9;
    }

    int gap_extend() {
        return 4;
    }

    int score_ab(char a, char b) {
        if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')
            || (a == 'C' && b == 'G') || (a == 'G' && b == 'C')) {
            return 5;
        } else if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) {
            return 1;
        } else if ((a == 'A' && b == 'A')
                   || (a == 'A' && b == 'C') || (a == 'C' && b == 'A')
                   || (a == 'A' && b == 'G') || (a == 'G' && b == 'A')
                   || (a == 'C' && b == 'C')
                   || (a == 'C' && b == 'U') || (a == 'U' && b == 'C')
                   || (a == 'G' && b == 'G')
                   || (a == 'U' && b == 'U')) {
            return -3;
        } else {
            return -1;
        }
    }

};

} // namespace mr3dp

#endif /* DP_SCORE_HPP_ */
