#ifndef DOC_SET_HPP
#define DOC_SET_HPP

#include <vector>
#include <cstdint>

const uint64_t MOD = 1000000000000000003ll;
const uint32_t ARR_SIZE = (1 << 16);
extern uint64_t pow2[ARR_SIZE];

class DocSet {
    public:
        std::vector<uint16_t> docs;
    
        DocSet() {}
        DocSet(std::vector<uint16_t> &_docs) : docs(_docs) {}
        DocSet(std::vector<uint16_t> &&_docs) : docs(_docs) {}
    
        bool operator==(const DocSet &o) const {
            return docs == o.docs;
        }
};
    
template<>
struct std::hash<DocSet> {
    std::size_t operator()(const DocSet &dc) const {
        size_t hash = 0;
        for (int doc : dc.docs) {
            hash += pow2[doc];
            if (hash >= MOD) {
                hash -= MOD;
            }
        }
        return hash;
    }
};

#endif