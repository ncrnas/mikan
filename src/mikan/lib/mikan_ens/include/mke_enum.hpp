#ifndef MKE_ENUM_HPP_
#define MKE_ENUM_HPP_

namespace mkens {

struct ToolIdx {
    enum idx {
        MR, PT, RH, TM, TS, SV
    };
    static const unsigned Count = 6;
};

struct SiteIdx {
    enum idx {
        MRAlg, MREng,
        PTDdg, PTDpx, PTOpn,
        RHMfe, RHNrm,
        TSCtx,
        SVSvm
    };
    static const unsigned Count = 9;
};

struct RNAIdx {
    enum idx {
        MRAlg, MREng,
        PTDdg,
        RHMfe, RHNrm,
        TMSvm,
        TSCtx,
        SVSvm
    };
    static const unsigned Count = 8;
};

} // namespace mkens

#endif //MKE_ENUM_HPP_
