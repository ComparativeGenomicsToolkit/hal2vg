#pragma once

#include <string>
#include <sstream>
#include <cassert>

inline std::string parse_subpath_name(const std::string& path_name, int64_t* out_start = nullptr, int64_t* out_end = nullptr) {

    std::string base_name = path_name;
    if (out_start) {
        *out_start = -1;
    }
    if (out_end) {
        *out_end = -1;
    }
    
    size_t first_length = 0;
    size_t start_offset = 0;
    while (true) {
        size_t sp = base_name.rfind("_sub_");
        if (sp != std::string::npos) {
            size_t up = base_name.rfind("_");
            if (up != std::string::npos && up > sp + 1) {
                int64_t start;
                int64_t end;
                try {
                    start = stol(base_name.substr(sp + 5, up - sp - 5));
                    end = stol(base_name.substr(up + 1));
                } catch (...) {
                    return base_name;
                }
                std::stringstream new_name;
                start_offset += start; // final offset is sum of all nested offsets
                if (first_length == 0) {
                    first_length = end - start;
                    assert(first_length > 0);
                } else {
                    // in the case of nested subpaths, the end coordinate will always
                    // be derived from the start, plus the length of the "top" path
                    end = start_offset + first_length;
                }
                if (out_start) {
                    *out_start = start_offset;
                }
                if (out_end) {
                    *out_end = end;
                }
                base_name = base_name.substr(0, sp);
            }
        } else {
            break;
        }
    }
    return base_name;
}

inline void resolve_subpath_naming(std::string& path_name) {
    int64_t sub_start;
    int64_t sub_end;
    std::string new_name = parse_subpath_name(path_name, &sub_start, &sub_end);
    if (sub_start != -1) {
        assert(new_name != path_name);
        path_name = new_name + "[" + std::to_string(sub_start) + "-" + std::to_string(sub_end) + "]";
    }
}

