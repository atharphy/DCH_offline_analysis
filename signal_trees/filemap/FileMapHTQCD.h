#ifndef FILEMAP_HTQCD_H
#define FILEMAP_HTQCD_H

#include "FileMap.h"

inline std::map<std::string, std::vector<std::string>> getFileMapHTQCD(const std::string& year) {
    auto fm = getFileMap(year);
    std::vector<std::string> names;
    for (const char* b : {"50to100", "100to200", "200to300", "300to500", "500to700", "700to1000", "1000to1500", "1500to2000", "2000toInf"})

        if (year == "2018" && (std::string(b) == "700to1000" || std::string(b) == "1000to1500")) continue;
        else names.push_back(std::string("QCD_HT") + b + "_" + year + ".root");
    fm["QCD"] = mcList(year, names);
    return fm;
}

#endif
