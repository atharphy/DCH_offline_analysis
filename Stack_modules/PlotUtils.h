#pragma once

#include <string>
#include <vector>

#include "TH1D.h"
#include "TROOT.h"

#include "StackConfig.h"

inline void divideByBinWidth(TH1D* h) {
    if (!h) return;
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
        double width = h->GetXaxis()->GetBinWidth(ib);
        if (width <= 0.0) continue;
        h->SetBinContent(ib, h->GetBinContent(ib) / width);
        h->SetBinError(ib, h->GetBinError(ib) / width);
    }
}

inline std::string channelFromHname(const std::string& hname) {
    std::string base = hname;

    std::vector<std::string> suffixes = {"_OS_Zwin", "_OS_Zveto", "_SS_Zwin", "_SS_Zveto"};
    for (const auto& suf : suffixes) {
        if (base.size() > suf.size() && base.compare(base.size() - suf.size(), suf.size(), suf) == 0) {
            base = base.substr(0, base.size() - suf.size());
            break;
        }
    }

    for (const auto& reg : regions) {
        std::string tag = "_" + reg;
        size_t pos = base.find(tag);
        if (pos != std::string::npos) {
            base = base.substr(0, pos);
            break;
        }
    }

    size_t last = base.find_last_of('_');
    if (last == std::string::npos) return "unknown";
    return base.substr(last + 1);
}

inline void CleanUpROOTMemory() {
    gROOT->cd();
    TH1::AddDirectory(kFALSE);
}
