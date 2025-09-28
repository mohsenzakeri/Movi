#ifndef MOVE_INTERVALS_HPP
#define MOVE_INTERVALS_HPP

#include <iostream>
#include <vector>
#include <cstdint>

#include "move_row.hpp"

struct MoveInterval {
    MoveInterval() {}

    MoveInterval& operator =(const MoveInterval& interval) {
        // Check self assignment
        if (this == &interval)
            return *this;
        run_start = interval.run_start;
        offset_start = interval.offset_start;
        run_end = interval.run_end;
        offset_end = interval.offset_end;
        return *this;
    }

    uint64_t run_start;
    uint64_t offset_start;
    uint64_t run_end;
    uint64_t offset_end;

    MoveInterval(uint64_t run_start_, uint64_t offset_start_, uint64_t run_end_, uint64_t offset_end_) {
        run_start = run_start_;
        offset_start = offset_start_;
        run_end = run_end_;
        offset_end = offset_end_;
    }

    void make_empty() {
        run_start = 1;
        offset_start = 0;
        run_end = 0;
        offset_end = 0;
    }

    bool is_empty() {
        return !((run_start < run_end) or (run_start == run_end and offset_start <= offset_end));
    }

    uint64_t count(std::vector<MoveRow>& rlbwt) {
        uint64_t row_count = 0;
        if (run_start == run_end) {
            row_count = offset_end - offset_start + 1;
        } else {
            row_count = (rlbwt[run_start].get_n() - offset_start) + (offset_end + 1);
            for (uint64_t k = run_start + 1; k < run_end; k ++) {
                row_count += rlbwt[k].get_n();
            }
        }
        return row_count;
    }

    friend std::ostream& operator<<(std::ostream& output, const MoveInterval& mi) {
        output << mi.run_start << ":" << mi.offset_start << " --- " << mi.run_end << ":" << mi.offset_end;
        return output;
    }

    bool operator==(const MoveInterval &m) {
        if (run_start == m.run_start and offset_start == m.offset_start and
            run_end == m.run_end and offset_end == m.offset_end)
            return true;
        return false;
    }
};

struct MoveBiInterval {
    MoveInterval fw_interval;
    MoveInterval rc_interval;
    uint64_t match_len;

    MoveBiInterval() {
        fw_interval = MoveInterval(1, 0, 0, 0);
        rc_interval = MoveInterval(1, 0, 0, 0);
        match_len = 0;
    }

    friend std::ostream& operator<<(std::ostream& output, const MoveBiInterval& mi_bi) {
        output << mi_bi.match_len << ": " << mi_bi.fw_interval << "\t" << mi_bi.rc_interval;
        return output;
    }
};

#endif