#ifndef FTAB_HPP
#define FTAB_HPP

#include <string>
#include <cinttypes>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>


class FTab : public std::map<std::string, std::pair<size_t,size_t>> {
    public:

    void load(std::ifstream& ifs) {
        (*this).clear();
        std::string line;
        std::string kmer;
        size_t s, e;
        while (std::getline(ifs, line)) {
            std::stringstream ss(line);
            ss >> kmer;
            ss >> s;
            ss >> e;
            (*this)[kmer] = std::make_pair(s, e);
        }
        k = kmer.size();
    }

    void serialize(std::ofstream& ofs) const {
        for (const auto& kv: *this) {
            ofs << kv.first << " " << kv.second.first << " " << kv.second.second << std::endl;
        }
    }

    size_t get_k() const { return k; }

    private:
    
    size_t k = 10;
};

#endif