#ifndef CAT_H
#define CAT_H

#include <string>

const char* numberToCat(int number) {
    static const char* cat_names[] = {
        "",

        "eeee", "eeem", "eeet", "eemm", "eemt", "eett",
        "emem", "emet", "emmm", "emmt", "emtt",
        "etet", "etmm", "etmt", "ettt",
        "mmmm", "mmmt", "mmtt",
        "mtmt", "mttt",
        "tttt",

        "eee", "eem", "eet",
        "eme", "emm", "emt",
        "ete", "etm", "ett",
        "mme", "mmm", "mmt",
        "mte", "mtm", "mtt",
        "tte", "ttm", "ttt",

        "ee", "em", "et", "mm", "mt", "tt",

        "e", "m", "t"
    };

    int total_cats = sizeof(cat_names) / sizeof(cat_names[0]);

    if (number < 0 || number >= total_cats) {
        return "Invalid";
    }

    return cat_names[number];
}

int cat_lepCount(const std::string& cat, char lep1, char lep2) {
    int count = 0;
    for (size_t i = 0; i < cat.length(); i++) {
        if (cat[i] == lep1 || cat[i] == lep2)
            count++;
    }
    return count;
}

#endif
