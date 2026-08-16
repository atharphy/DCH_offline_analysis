#pragma once

#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TMemFile.h"
#include "TObject.h"

#include "HistCache.h"
#include "Labels.h"
#include "SystematicSources.h"

inline void mergeInputFilesByProcess(std::map<std::string, std::vector<TFile*>>& handles, const std::unordered_set<std::string>& selectedVariables, const std::vector<SystSource>& activeSources, const std::string& namePrefix = "merged_") {
    for (auto& processEntry : handles) {
        std::cout << "  merging " << processEntry.first << " from " << processEntry.second.size() << " files" << std::endl;
        std::unordered_map<std::string, TH1D*> sums;

        for (TFile* input : processEntry.second) {
            if (!input || input->IsZombie()) continue;

            for (const std::string& name : keysInFile(input)) {
                if (name.rfind("h_", 0) != 0) continue;
                if (!selectedVariables.count(getVariableFromHname(name))) continue;
                if (!DRAW_BANDS && isSystematicVariantName(name, activeSources)) continue;

                TH1D* source = dynamic_cast<TH1D*>(input->Get(name.c_str()));
                if (!source) continue;

                auto found = sums.find(name);
                if (found == sums.end()) {
                    TH1D* sum = dynamic_cast<TH1D*>(source->Clone(name.c_str()));
                    sum->SetDirectory(nullptr);
                    sum->Sumw2();
                    sums[name] = sum;
                } else {
                    found->second->Add(source);
                }
            }
        }

        const std::size_t nCached = sums.size();
        TMemFile* merged = new TMemFile((namePrefix + processEntry.first + ".root").c_str(), "RECREATE");
        merged->cd();
        for (auto& item : sums) {
            item.second->Write(item.first.c_str(), TObject::kOverwrite);
            delete item.second;
        }
        merged->Write();

        for (TFile* input : processEntry.second) {
            if (!input) continue;
            fileKeys.erase(input);
            input->Close();
            delete input;
        }

        processEntry.second.clear();
        processEntry.second.push_back(merged);
        std::cout << "    cached " << nCached << " histograms" << std::endl;
    }
}
