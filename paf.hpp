#pragma once

#include <string>
#include <vector>
#include <iostream>

using namespace std;

struct PafLine {
    string query_name;
    int64_t query_len;
    int64_t query_start;
    int64_t query_end;
    char strand;
    string target_name;
    int64_t target_len;
    int64_t target_start;
    int64_t target_end;
    int64_t num_matching;
    int64_t num_bases;
    int64_t mapq;
    string cigar;
};

inline vector<string> split_delims(const string &s, const string& delims, vector<string> &elems) {
    size_t start = string::npos;
    for (size_t i = 0; i < s.size(); ++i) {
        if (delims.find(s[i]) != string::npos) {
            if (start != string::npos && i > start) {
                elems.push_back(s.substr(start, i - start));
            }
            start = string::npos;
        } else if (start == string::npos) {
            start = i;
        }
    }
    if (start != string::npos && start < s.size()) {
        elems.push_back(s.substr(start, s.size() - start));
    }
    return elems;
}

inline PafLine parse_paf_line(const string& paf_line) {
    vector<string> toks;
    split_delims(paf_line, "\t\n", toks);
    assert(toks.size() > 12);

    PafLine paf;
    paf.query_name = toks[0];
    paf.query_len = stol(toks[1]);
    paf.query_start = stol(toks[2]);
    paf.query_end = stol(toks[3]);
    assert(toks[4] == "+" || toks[4] == "-");
    paf.strand = toks[4][0];
    paf.target_name = toks[5];
    paf.target_len = stol(toks[6]);
    paf.target_start = stol(toks[7]);
    paf.target_end = stol(toks[8]);
    paf.num_matching = stol(toks[9]);
    paf.num_bases = stol(toks[10]);
    paf.mapq = stol(toks[11]);

    for (size_t i = 12; i < toks.size(); ++i) {
        if (toks[i].compare(0, 3, "cg:Z:") == 0) {
            paf.cigar = toks[i].substr(5);
            break;
        }
    }

    return paf;
}

inline ostream& operator<<(ostream& os, const PafLine& paf) {
    os << paf.query_name << "\t" << paf.query_len << "\t" << paf.query_start << "\t" << paf.query_end << "\t"
       << string(1, paf.strand) << "\t"
       << paf.target_name << "\t" << paf.target_len << "\t" << paf.target_start << "\t" << paf.target_end << "\t"
       << paf.num_matching << "\t" << paf.num_bases << "\t" << paf.mapq;
    return os;
}
