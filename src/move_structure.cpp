#include "move_structure.hpp"
#include "classifier.hpp"

uint64_t pow2[ARR_SIZE];

MoveStructure::MoveStructure(MoviOptions* movi_options_) {
    movi_options = movi_options_;
    onebit = false; // This is not used any more as the onebit modes is deprecated
    no_ftab = 0;
    all_initializations = 0;
}

MoveStructure::MoveStructure(MoviOptions* movi_options_, uint16_t splitting_, bool constant_) {
    movi_options = movi_options_;
    onebit = false; // This is not used any more as the onebit modes is deprecated
    splitting = splitting_;
    constant = constant_;
    no_ftab = 0;
    all_initializations = 0;
    reconstructed = false;

    std::string bwt_filename = movi_options->get_ref_file() + std::string(".bwt");
    // std::ifstream bwt_file(bwt_filename);

    if (movi_options->is_thresholds()) {
        std::string thr_filename = movi_options->get_ref_file() + std::string(".thr_pos");
        read_thresholds(thr_filename, thresholds);
    }
    build();
}

std::vector<MoveRow> MoveStructure::get_rlbwt() {
    return rlbwt;
}

/*char MoveStructure::compute_char(uint64_t idx) {
    if (verbose) {
        std::cerr << idx << " bit1_begin: " << bit1_begin << "\n";
        std::cerr << idx << " bit1_after_eof: " << bit1_after_eof << "\n";
        std::cerr << "alphabet.size(): " << alphabet.size() << "\n";
    }
    char row_c;
    if (idx < eof_row) {
        row_c = idx%2 == 0 ? alphabet[bit1_begin] : alphabet[(1 + bit1_begin) % 2];
    } else if (idx == eof_row) {
        row_c = END_CHARACTER;
    } else{
        row_c = (idx - eof_row + 1)%2 == 0 ? alphabet[bit1_after_eof] : alphabet[(1 + bit1_after_eof) % 2];
    }
    return row_c;
}*/

uint64_t MoveStructure::LF(uint64_t row_number, uint64_t alphabet_index) {
    uint64_t lf = 0;
    lf += 1;
    for (uint64_t i = 0; i < alphabet_index; i++) {
        lf += counts[i];
    }
    auto& occ_rank = *occs_rank[alphabet_index];
    lf += static_cast<uint64_t>(occ_rank(row_number));
    return lf;
}

/*std::string MoveStructure::reconstruct() {
    if (bwt_string == "") {
        for (uint32_t i = 0; i < r; i++) {
            for (uint32_t j = 0; j < get_n(i); j++)
                bwt_string += rlbwt[i].get_c();
        }
    }
    if (!reconstructed) {
        orig_string = "";
        orig_string += static_cast<char>(END_CHARACTER);
        uint64_t alphabet_index = rlbwt[idx_].get_c();
        for (uint64_t bwt_row = 0; bwt_row != end_bwt_row; bwt_row = LF(bwt_row, alphabet_index)) {
            orig_string = bwt_string[bwt_row] + orig_string;
            uint64_t alphabet_index = alphamap[static_cast<uint64_t>(bwt_string[bwt_row])];
        }
        reconstructed = true;
    }
    return orig_string;
}*/

/*uint64_t MoveStructure::naive_lcp(uint64_t row1, uint64_t row2) {
    if (row1 >= length or row2 >= length)
        return 0;

    if (!reconstructed)
        orig_string = reconstruct();
    uint64_t lcp = 0;
    uint64_t sa1 = naive_sa(row1);
    uint64_t sa2 = naive_sa(row2);

    while (orig_string[sa1 + lcp] == orig_string[sa2 + lcp]) {
        lcp += 1;
    }

    return lcp;
}*/

/*uint64_t MoveStructure::naive_sa(uint64_t bwt_row) {
    uint64_t sa = 0;
    for (; bwt_row != end_bwt_row; bwt_row = LF(bwt_row, alphabet_index)) {
        sa += 1;
        uint64_t alphabet_index = alphamap[static_cast<uint64_t>(bwt_string[bwt_row])];
    }
    return sa;
}*/

uint32_t MoveStructure::compute_index(char row_char, char lookup_char) {
    uint32_t alpha_index = alphamap[lookup_char];
    return alpha_index;
    /*if (lookup_char < row_char)
        return alpha_index;
    else
        return alpha_index+1;*/
}

uint16_t MoveStructure::LF_move(uint64_t& offset, uint64_t& i, uint64_t id) {
    /* if (movi_options->is_verbose()) {
        std::cerr << "\t in LF:\n";
        std::cerr << "\t \t i: " << i << " offset: " << offset << "\n";
    } */
    auto& row = rlbwt[i];
    auto idx = (id == std::numeric_limits<uint64_t>::max()) ? get_id(i) : id;
    if (idx >= r) {
        std::cerr << "This should not happen.\n";
        std::cerr << "idx:::" << idx << " i:" << i << "\n";
        exit(0);
    }
    offset = get_offset(i) + offset;
    uint16_t ff_count = 0;
    /* if (movi_options->is_verbose()) {
        std::cerr << "\t \t i: " << i << " offset: " << offset << " idx: " << idx << "\n";
    } */

    if (idx < r - 1 && offset >= get_n(idx)) {
        uint64_t idx_ = fast_forward(offset, idx, 0);
        idx += idx_;
        if (idx_ >= std::numeric_limits<uint16_t>::max()) {
            std::cerr << "Number of fast forwards for a query was greater than 2^16: " << idx_ << "\n";
            std::cerr << "offset: " << offset << "\n";
            std::cerr << "idx: " << idx << "\n";
            exit(0);
        }
        ff_count = static_cast<uint16_t>(idx_);
    }

    /* if (movi_options->is_verbose()) {
        std::cerr << "\t \t after fast forward:\n";
        std::cerr << "\t \t i: " << i << " offset: " << offset << " idx: " << idx << "\n";
    } */

    if (movi_options->is_logs()) {
        if (ff_counts.find(ff_count) != ff_counts.end())
            ff_counts[ff_count] += 1;
        else
            ff_counts[ff_count] = 1;
    }
    i = idx;
    return ff_count;
}

// Finds documents corresponding to rows in BWT
void MoveStructure::build_doc_pats() {
    doc_pats.resize(length);
    uint64_t offset = 0;
    uint64_t index = 0;
    uint64_t SA_val = length;
    uint32_t doc_offset_ind = num_docs - 1;
    for (uint64_t i = 0; i < length; i++) {
        if (i % 1000000 == 0)
            std::cerr << "Finding suffix array entries: " << i << "\r";
        SA_val--;
        if (doc_offset_ind > 0 && SA_val < doc_offsets[doc_offset_ind - 1]) {
            doc_offset_ind--;
        }
        uint64_t row_ind = run_offsets[index] + offset;
        doc_pats[row_ind] = doc_ids[doc_offset_ind];
        LF_move(offset, index);
    }
    std::cerr << "\n";
}

void MoveStructure::fill_run_offsets() {
    // Fill run offset information for finding all SA entries                                
    uint64_t run_offset = 0;
    run_offsets.resize(r);
    for (uint64_t i = 0; i < r; i++) {
        run_offsets[i] = run_offset;
        run_offset += rlbwt[i].get_n();
    }
}

// Finds document sets for each run in the rlbwt.
void MoveStructure::build_doc_sets() {
    // Stores set of documents appearing in rows in each run.
    std::unordered_map<DocSet, uint32_t> unique; 
    doc_set_inds.resize(r);

    uint64_t unique_cnt = 0;
    uint16_t temp[MAX_RUN_LENGTH];
    for (uint64_t i = 0; i < r; i++) {
        if (i % 100000000 == 0) {
            std::cerr << "Processed " << i << " runs, " << unique_cnt << " unique doc sets so far\r";
        }

        uint64_t n = rlbwt[i].get_n();
        uint16_t ptr = 0;
        for (uint64_t j = 0; j < n; j++) {
            uint64_t row_ind = run_offsets[i] + j;
            uint16_t doc = doc_pats[row_ind];
            temp[ptr++] = doc;
        }

        // Sort documents and initialize unique elements as DocSet
        std::sort(temp, temp + ptr);
        std::vector<uint16_t> docs;
        for (size_t j = 0; j < ptr; j++) {
            if (j == 0 || temp[j] != temp[j - 1]) {
                docs.push_back(temp[j]);
            }
        }
        DocSet cur(std::move(docs));

        auto it = unique.find(cur);
        if (it != unique.end()) {
            doc_set_inds[i] = it->second;
        } else {
            doc_set_inds[i] = unique_cnt;
            unique[cur] = unique_cnt++;
        }    
    }

    std::cerr << std::endl;
    std::cerr << "Done building and hashing doc sets" << std::endl;
    unique_doc_sets.resize(unique.size());
    for (auto it = unique.begin(); it != unique.end(); it++) {
        unique_doc_sets[it->second] = it->first.docs;
    }
}

uint32_t MoveStructure::hash_collapse(std::unordered_map<DocSet, uint32_t> &keep_set, DocSet &bv) {
    // Remove documents with highest multiplicity one at a time
    // until doc set is one of the kept ones.
    /*std::vector<std::pair<uint16_t, uint16_t>> cnts(num_docs);
    for (size_t i = 0; i < num_docs; i++) {
        cnts[i].second = i;
    }
    
    DocSet cur(bv);
    //std::sort(cnts.begin(), cnts.end());
    std::random_shuffle(cnts.begin(), cnts.end());
    for (size_t i = 0; i <= num_docs; i++) {
        if (keep_set.count(cur)) {
            return keep_set[cur];
        }
        if (i < num_docs) {
            cur.unset(cnts[i].second);
        }
    }
    assert(false);*/
    
    /*std::vector<int> min_inds;
    int min_diff = INT_MAX;
    for (auto it = keep_set.begin(); it != keep_set.end(); it++) {
        int diff = 0;
        for (int i = 0; i < num_docs; i++) {
            if (bv[i] != it->first.bv[i]) {
                diff++;
            }
        }
        if (diff <= min_diff) {
            if (diff < min_diff) min_inds.clear();
            min_inds.push_back(it->second);
            min_diff = diff;
        }
    }
    return min_inds[rand() % min_inds.size()];*/
    return 0;
}

void MoveStructure::compress_doc_sets() {
    // How many doc sets to keep.
    int take = (1 << 16);

    // Get doc set counts.
    doc_set_cnts.resize(unique_doc_sets.size());
    for (size_t i = 0; i < r; i++) {
        doc_set_cnts[doc_set_inds[i]]++;
    }

    /*std::vector<uint64_t> set_dist(num_species + 1);
    for (size_t i = 0; i < unique_doc_sets.size(); i++) {
        set_dist[unique_doc_sets[i].size()] += doc_set_cnts[i];
    }
    */

    // Sort document sets by their frequency (and ensuring singletons are put first).
    std::vector<std::tuple<bool, uint64_t, uint32_t>> sorted(doc_set_cnts.size());
    for (size_t i = 0; i < doc_set_cnts.size(); i++) {
        sorted[i] = {unique_doc_sets[i].size() == 1, doc_set_cnts[i], i};
    }
    std::sort(sorted.begin(), sorted.end(), std::greater<>());
    std::cerr << "Sorted document sets by frequency (ensuring singletons first)" << std::endl;

    /*std::ofstream sort_freq_out("../indices/pseudomonadota/sort_freq_out.txt");
    for (size_t i = 0; i < sorted.size(); i++) {
        auto [singleton, freq, ind] = sorted[i];
        sort_freq_out << freq << " ";
        for (uint16_t doc : unique_doc_sets[ind]) {
            sort_freq_out << doc << " ";
        }
        sort_freq_out << "\n";
    }
    sort_freq_out.close();*/
    
    // Only keep the most frequent document sets.
    std::vector<std::vector<uint16_t>> keep(take);
    std::vector<bool> in_keep(unique_doc_sets.size());
    std::vector<uint32_t> compress_to(unique_doc_sets.size(), take);
    for (size_t i = 0; i < take; i++) {
        int ind = std::get<2>(sorted[i]);
        keep[i] = unique_doc_sets[ind];
        in_keep[ind] = true;
        compress_to[ind] = i;
    }
    
    uint64_t missing_cnt = 0;
    for (size_t i = 0; i < r; i++) {
        if (!in_keep[doc_set_inds[i]]) {
            missing_cnt++;
        }
        doc_set_inds[i] = compress_to[doc_set_inds[i]];
    }
    unique_doc_sets = keep;
    std::cerr << "Fraction of runs without doc set: " << (double) missing_cnt / r << std::endl;
}

// Returns if x is an ancestor of y.
bool MoveStructure::is_ancestor(uint16_t x, uint16_t y) {
    return t_in[x] <= t_in[y] && t_out[x] >= t_out[y];
}

uint16_t MoveStructure::LCA(uint16_t x, uint16_t y) {
    if (is_ancestor(x, y)) return x;
    int cur = x;
    for (int i = 15; i >= 0; i--) {
        if (!is_ancestor(bin_lift[i][cur], y)) {
            cur = bin_lift[i][cur];
        }
    }
    return bin_lift[0][cur];
}

void MoveStructure::dfs_times(uint16_t cur, uint16_t &t) {
    t_in[cur] = t++;
    for (int child : tree[cur]) {
        dfs_times(child, t);
    }
    t_out[cur] = t++;
}

void MoveStructure::build_tree_doc_sets() {
    std::ifstream fin(movi_options->get_index_dir() + "/doc_set_similarities.txt");
    double distmat[num_species * (num_species - 1) / 2];
    int ind = 0;
    for (int i = 0; i < num_species; i++) {
        for (int j = 0; j < num_species; j++) {
            double cur;
            fin >> cur;
            if (j > i) {
                distmat[ind++] = cur;
            }
        }
    }
    for (int i = 0; i < num_species * (num_species - 1) / 2; i++) {
        distmat[i] = 1 - distmat[i] / r;
    }
    fin.close();

    std::cerr << "Read in similarities matrix" << std::endl;

    // Set up tree for hierarchical clustering
    int nodes = num_species * 2 - 1;
    tree.resize(nodes);
    tree_doc_sets.resize(nodes);
    bin_lift.resize(16, std::vector<uint16_t>(nodes, nodes - 1));

    for (int i = 0; i < num_species; i++) {
        tree_doc_sets[i].emplace_back(i);
    }

    // Cluster documents based on doc set similarities.
    int *merge = new int[2 * (num_species - 1)];
    double *height = new double[num_species - 1];
    hclust_fast(num_species, distmat, HCLUST_METHOD_AVERAGE, merge, height);

    std::cerr << "Clustered documents using hierarchical clustering algorithm" << std::endl;

    // Go through merges in hierarchical clustering, construct doc sets at each node.
    std::vector<int> last_merge(num_species, 0);
    for (int i = 0; i < num_species; i++) last_merge[i] = i;
    for (int i = 0; i < num_species - 1; i++) {
        int node1 = merge[i];
        if (node1 < 0) node1 = -node1 - 1;
        else node1 += num_species - 1;

        int node2 = merge[num_species - 1 + i];
        if (node2 < 0) node2 = -node2 - 1;
        else node2 += num_species - 1;
        
        for (int j = 0; j < num_species; j++) {
            if (last_merge[j] == node1 || last_merge[j] == node2) {
                last_merge[j] = num_species + i;
                tree_doc_sets[num_species + i].push_back(j);
            }
        }
        
        tree[num_species + i].push_back(node1);
        tree[num_species + i].push_back(node2);
        bin_lift[0][node1] = num_species + i;
        bin_lift[0][node2] = num_species + i;
    }

    for (int i = 1; i < 16; i++) {
        for (int j = 0; j < nodes; j++) {
            bin_lift[i][j] = bin_lift[i - 1][bin_lift[i - 1][j]];
        }
    }
    uint16_t timer = 0;
    t_in.resize(nodes); t_out.resize(nodes);
    dfs_times(nodes - 1, timer);

    std::cerr << "Built compression tree" << std::endl;

    // Compress doc sets by LCA in tree
    std::vector<uint32_t> compress_to(unique_doc_sets.size());
    std::vector<bool> in_keep(unique_doc_sets.size());
    for (size_t i = 0; i < unique_doc_sets.size(); i++) {
        int lca = -1;
        for (uint16_t doc : unique_doc_sets[i]) {
            if (lca == -1) lca = doc;
            else lca = LCA(lca, doc);
        }
        assert(lca != -1);
        compress_to[i] = lca;
        // debug_out << unique_doc_sets[i] << " " << tree_doc_sets[lca] << "\n";
    }
    
    for (size_t i = 0; i < r; i++) {
        doc_set_inds[i] = compress_to[doc_set_inds[i]];
    }
    unique_doc_sets = tree_doc_sets;
    std::cerr << "Completed tree compression of colors" << std::endl;
}

void MoveStructure::build_doc_set_similarities() {
    // Get doc set counts.
    doc_set_cnts.resize(unique_doc_sets.size());
    for (size_t i = 0; i < r; i++) {
        doc_set_cnts[doc_set_inds[i]]++;
    }

    std::vector<std::vector<uint64_t>> similarities(num_species, std::vector<uint64_t>(num_species));
    for (size_t i = 0; i < unique_doc_sets.size(); i++) {
        std::vector<uint16_t> &docs = unique_doc_sets[i];
        for (size_t j = 0; j < docs.size(); j++) {
            for (size_t k = j + 1; k < docs.size(); k++) {
                similarities[docs[j]][docs[k]] += doc_set_cnts[i];
            }
        }
    }

    std::string fname = movi_options->get_index_dir() + "/doc_set_similarities.txt";
    std::ofstream fout(fname);
    for (int i = 0; i < num_species; i++) {
        for (int j = 0; j < num_species; j++) {
            fout << similarities[i][j] << " ";
        }
        fout << "\n";
    }
    fout.close();
}

std::string MoveStructure::reconstruct_lf() {
    orig_string = "";
    uint64_t offset = 0;
    uint64_t run_index = 0;
    uint64_t i = 0;
    uint64_t ff_count_tot = 0;

    uint64_t total_elapsed = 0;
    for (; run_index != end_bwt_idx; ) {
        auto begin = std::chrono::system_clock::now();
        ff_count_tot += LF_move(offset, run_index);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        total_elapsed += elapsed.count();

        if (i % 10000 == 0)
            std::cerr << i << "\r";
        i += 1;
        // orig_string = rlbwt[run_index].get_c() + orig_string;
    }
    std::printf("Time measured for reconstructing the original text: %.3f seconds.\n", total_elapsed * 1e-9);

    std::cerr << "Finished reconstructing the original string.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
    return orig_string;
}

void MoveStructure::sequential_lf() {
    uint64_t line_index = 0;
    uint64_t row_index = 0;
    uint64_t ff_count_tot = 0;

    uint64_t total_elapsed = 0;
    for (uint64_t row_index = 0; row_index < r; row_index++) {
        auto& current = rlbwt[row_index];
        for (uint64_t j = 0; j < current.get_n(); j ++) {
            if (line_index % 10000 == 0)
                std::cerr << line_index << "\r";
            uint64_t offset = j;
            uint64_t i = row_index;
            auto begin = std::chrono::system_clock::now();
            ff_count_tot += LF_move(offset, i);
            auto end = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            total_elapsed += elapsed.count();
            line_index += 1;
        }

    }
    std::printf("Time measured for LF-mapping of all the characters: %.3f seconds.\n", total_elapsed * 1e-9);

    std::cerr << line_index << "\n";
    std::cerr << "Finished performing LF-mapping for all the BWT characters.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
}

void MoveStructure::random_lf() {
    uint64_t ff_count_tot = 0;

    // generate the random order from 1 to length
    std::vector<uint64_t> random_order(length);
    for (uint64_t i = 0; i < length; i++) {
        random_order[i] = i;
    }
    std::srand(time(0));
    std::random_shuffle(random_order.begin(), random_order.end());

    // find the n and id for each random BWT row
    std::vector<uint64_t> n_to_id(length);
    std::vector<uint64_t> id_to_p(r);
    uint64_t current_n = 0;
    for (uint64_t i = 0; i < r; i++) {
        id_to_p[i] = current_n;
        for (uint64_t j = 0; j < get_n(i); j++) {
            n_to_id[current_n] = i;
            current_n += 1;
        }
    }


    uint64_t total_elapsed = 0;
    for (uint64_t i = 0; i < length; i++) {
        if (i % 10000 == 0)
            std::cerr << i << "\r";

        // uint64_t n = std::rand() % r;
        uint64_t n = random_order[i];
        uint64_t id = n_to_id[n];
        uint64_t offset = n - id_to_p[id];
        auto begin = std::chrono::system_clock::now();
        ff_count_tot += LF_move(offset, id);
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        total_elapsed += elapsed.count();
    }
    std::printf("Time measured for LF-mapping of all the characters in the random order: %.3f seconds.\n", total_elapsed * 1e-9);

    std::cerr << length << "\n";
    std::cerr << "Finished performing LF query for all the BWT characters.\n";
    std::cerr << "Total fast forward: " << ff_count_tot << "\n";
}

#define my_prefetch_rr(address) __builtin_prefetch((void *)address, 0, 1)

uint64_t MoveStructure::get_id(uint64_t idx) {
#if NO_EXTRA_TABLE
    return rlbwt[idx].get_id();
#endif
#if BLOCKED_MODE
    if (idx != end_bwt_idx) {
        uint64_t block_number = idx / BLOCK_SIZE;
        return rlbwt[idx].get_id() + static_cast<uint64_t>(id_blocks[rlbwt[idx].get_c()][block_number]) + first_runs[rlbwt[idx].get_c() + 1];
    }
    else
        return rlbwt[idx].get_id();
#endif
#if TALLY_MODE
    uint64_t id = 0;
    // The $ always goes to row/run 0
    if (idx == end_bwt_idx) {
        return 0;
    }

    uint8_t char_index = rlbwt[idx].get_c();
    uint64_t tally_a = idx / tally_checkpoints;

    // The id of the last run is always stored at the last tally_id row
    if (idx == r - 1) {
        uint64_t tally_ids_len = tally_ids[char_index].size();
        return tally_ids[char_index][tally_ids_len - 1].get();
    }


    // If we are at a checkpoint, simply return the stored id
    if (idx % tally_checkpoints == 0) {
        id = tally_ids[char_index][tally_a].get();
        return tally_ids[char_index][tally_a].get();
    }

    // find the closest check point
    uint64_t tally_b = tally_a + 1;

    // Sanity Check to see if we fall off the tally vector
    // if (tally_b >= tally_ids[char_index].size()) {
    //     dbg << "ta: " << tally_a << " tb: " << tally_b <<
    //            " tsize: " << tally_ids[char_index].size() << "\n";
    // }

    uint64_t next_check_point = tally_b * tally_checkpoints;
    // If we are at the end, just look at the last run
    if (tally_b * tally_checkpoints >= r) {
        next_check_point = r - 1;
    }

    int64_t prev_check_point = tally_a * tally_checkpoints;

    // Find the closest checkpoint which is either before or after the current row
    // For the first run at 0, the previous tallies are not stored, so we have to look below
    // Since the above direction never helps, for now always set forward_direction to true
    bool forward_direciton = true;
    if (forward_direciton or (tally_a == 0 or tally_b - idx <= idx - tally_a)) {
        uint64_t rows_until_tally = 0;

        // The id for the run at the next_check_point
        id = tally_ids[char_index][tally_b].get();

        uint64_t last_count = 0;
        uint64_t last_id = r;

        // TODO: Could prefetching help with performance?
        // Prelimninary results were not hopeful
        // for (uint64_t i = 0; i < tally_checkpoints; i++) {
        //     my_prefetch_rr((void*)(&(rlbwt[0]) + idx + i));
        //     my_prefetch_rr((void*)(&(rlbwt[0]) + id - i));
        // }


        // Look for the rows between idx and the next_check_point
        // Count how many rows with the same character exists between
        // the current run head and the next run head for which we know the id
        // dbg << "from: " << idx << " to: " << next_check_point << "\n";
        for (uint64_t i = idx; i < next_check_point; i++) {
            // Only count the rows with the same character
            if (get_char(i) == get_char(idx)) {
                rows_until_tally += get_n(i);
                last_id = i;
            }
        }

        // The id stored at the checkpoint is actually the id for the current run
        // because there was no run with the same character between the current run and the next_check_point
        if (last_id == idx and get_char(idx) != get_char(next_check_point)) {
            return id;
        }
        if (last_id == r) {
            std::cerr << "last_id should never be equal to r.\n";
            exit(0);
        }

        // We should know what will be the offset of the run head in the destination run (id)
        uint16_t offset = get_offset(next_check_point);

        // At the checkpoint, one id is stored for each character
        // For the character of the checkpoint, the id of the checkpoint is stored
        // For other characters, the id of the last run before the checkpoint with that character is stored
        // So, we have to decrease the number of rows of the last run as they are not between the current run
        // and the run for which we know the id
        // Also, for other characters, offset should be stored based on that run (not the checkpoint run)
        if (get_char(idx) != get_char(next_check_point)) {
            rows_until_tally -= get_n(last_id);
            offset = get_offset(last_id);
        }

        // Sanity check, the offset in the destination run should be always smaller than the length of that run
        if (offset >= get_n(id)) {
            std::cout << "offset: " << offset << " n: " << get_n(id) << "\n"; // << dbg.str() << "\n";
            exit(0);
        }
        if (offset >= rows_until_tally) {
            return id;
        } else {
            rows_until_tally -= (offset + 1);
            id -= 1;
        }

        while (rows_until_tally != 0) {
            if (rows_until_tally >= get_n(id)) {
                rows_until_tally -= get_n(id);
                id -= 1;
                my_prefetch_rr((void*)(&(rlbwt[0]) + id));
            } else {
                rows_until_tally = 0;
            }
        }
    } else {
        uint64_t rows_until_tally = 0;

        // The id for the run at the previous_check_point
        id = tally_ids[char_index][tally_a].get();

        uint64_t last_count = 0;
        uint64_t last_id = r;

        // Look for the rows between idx and the prev_check_point
        // Count how many rows with the same character exists between
        // the current run head and the prev run head for which we know the id
        // dbg << "from: " << idx << " to: " << prev_check_point << "\n";
        for (int64_t i = idx - 1; i >= prev_check_point; i--) {
            // Only count the rows with the same character
            if (get_char(i) == get_char(idx)) {
                rows_until_tally += get_n(i);
                last_id = i;
            }
        }

        // The last checkpoint is storing an id for a run above, so we have to count the
        // number of rows in that run too
        // We have to first find that run with the same character
        if (get_char(idx) != get_char(prev_check_point)) {
            last_id = prev_check_point;
            while (get_char(idx) != get_char(last_id)) {
                last_id -= 1;
            }
            rows_until_tally += get_n(last_id);
        }

        if (last_id == r) {
            std::cerr << idx << " " << tally_a << "\n";
            std::cerr << "last_id should never be equal to r.\n";
            exit(0);
        }

        uint16_t offset = get_offset(last_id);

        // Sanity check, the offset in the destination run should be always smaller than the length of that run
        if (offset >= get_n(id)) {
            std::cout << "idx: " << idx << " last_id: " << last_id << " id: " << id << " prev_check_point: " << prev_check_point << "\n";
            std::cout << "offset: " << offset << " n: " << get_n(id) << "\n"; // << dbg.str() << "\n";
            exit(0);
        }
        if ((get_n(id) - offset - 1) >= rows_until_tally) {
            return id;
        } else {
            rows_until_tally -= (get_n(id) - offset);
            id += 1;
        }

        while (rows_until_tally != 0) {
            if (rows_until_tally >= get_n(id)) {
                rows_until_tally -= get_n(id);
                id += 1;
                my_prefetch_rr((void*)(&(rlbwt[0]) + id));
            } else {
                rows_until_tally = 0;
            }
        }
    }

    return id;
#endif
}

char MoveStructure::get_char(uint64_t idx) {
    if (idx == end_bwt_idx)
        return '$';
    else
        return alphabet[rlbwt[idx].get_c()];
}

#if SPLIT_MAX_RUN
uint64_t MoveStructure::get_n(uint64_t idx) {
    return rlbwt[idx].get_n();
}

uint64_t MoveStructure::get_offset(uint64_t idx) {
    return rlbwt[idx].get_offset();
}
#endif

#if SPLIT_THRESHOLDS_TRUE
uint64_t MoveStructure::get_thresholds(uint64_t idx, uint32_t alphabet_index) {
    return rlbwt[idx].get_threshold(alphabet_index) == 0 ? 0 : get_n(idx);
}
#endif

#if SPLIT_THRESHOLDS_FALSE
uint64_t MoveStructure::get_n(uint64_t idx) {
    if (rlbwt[idx].is_overflow_n()) {
        return n_overflow[rlbwt[idx].get_n()];
    } else {
        return rlbwt[idx].get_n();
    }
}

uint64_t MoveStructure::get_offset(uint64_t idx) {
    if (rlbwt[idx].is_overflow_offset()) {
        return offset_overflow[rlbwt[idx].get_offset()];
    } else {
        return rlbwt[idx].get_offset();
    }
}

// for all the threshold related functions
uint64_t MoveStructure::get_thresholds(uint64_t idx, uint32_t alphabet_index) {
    if (rlbwt[idx].is_overflow_thresholds()) {
        return thresholds_overflow[get_rlbwt_thresholds(idx, alphabet_index)][alphabet_index];
    } else {
        return get_rlbwt_thresholds(idx, alphabet_index);
    }
}

uint16_t MoveStructure::get_rlbwt_thresholds(uint64_t idx, uint16_t i) {
    if (i >= alphabet.size() - 1) {
        std::cerr << "get_thresholds: " << i << " is greater than or equal to " << alphabet.size() - 1 << "\n";
        exit(0);
    }

    uint8_t status = rlbwt[idx].get_threshold_status(i);
    switch (status) {
        case 0: return 0; break;
        case 1: return rlbwt[idx].get_threshold(); break;
        case 3: return get_n(idx); break;
        default:
            std::cerr << "Undefined status for thresholds status: " << status << "\n";
            exit(0);
    }

    std::cerr << "Undefined behavior!\n";
    exit(0);
}

void MoveStructure::set_rlbwt_thresholds(uint64_t idx, uint16_t i, uint16_t value) {
    if (i >= alphabet.size() - 1) {
        std::cerr << "get_thresholds: " << i << " is greater than or equal to " << alphabet.size() - 1 << "\n";
        exit(0);
    }

    uint8_t status = 0;
    if (value == 0) {
        status = 0;
    } else if (value == get_n(idx)) {
        status = 3;
    } else {
        status = 1;
        // [TODO] Not all the states where the multiple non-trivial thresholds exists are checked here
        if (rlbwt[idx].get_threshold() != value and
            rlbwt[idx].get_threshold() != 0 and
            rlbwt[idx].get_threshold() != get_n(idx) and
            !rlbwt[idx].is_overflow_thresholds()) {
            // std::cerr << "idx: " << idx << " i: " << i << " value: " << value << "\n";
            // std::cerr << rlbwt[idx].get_threshold() << " " << !rlbwt[i].is_overflow_thresholds() << "\n";
            // std::cerr << "There are more than 1 non-trivial threshold values.\n";
            // exit(0);
            rlbwt[i].set_overflow_thresholds();
            return;
        }
        rlbwt[idx].set_threshold(value);
    }
    rlbwt[idx].set_threshold_status(i, status);
}
#endif // for all the threshold related functions

void MoveStructure::build_rlbwt() {
    std::ifstream bwt_file(movi_options->get_bwt_file());
    bwt_file.clear();
    bwt_file.seekg(0,std::ios_base::end);
    std::streampos end_pos = bwt_file.tellg();
    if (movi_options->is_verbose())
        std::cerr << "end_pos: " << end_pos << "\n";
    std::cerr << static_cast<uint64_t>(end_pos) << "\n";
    bwt_file.seekg(0);
    char current_char = bwt_file.get();
    char last_char = current_char;
    r = 0;
    size_t len = 0;
    
    std::ofstream len_file(movi_options->get_bwt_file() + ".len", std::ios::out | std::ios::binary);
    std::ofstream heads_file(movi_options->get_bwt_file() + ".heads");
    while (current_char != EOF) {
        if (r % 10000 == 0)
            std::cerr << r << "\r";
        if (current_char != last_char) {
            r += 1;
            // write output
            heads_file << last_char;
            len_file.write(reinterpret_cast<char*>(&len), 5);
            len = 0;
        } 
        len += 1;
        last_char = current_char;
        current_char = bwt_file.get();
        
    }

    // write output
    heads_file << last_char;
    len_file.write(reinterpret_cast<char*>(&len), 5);

    heads_file.close();
    len_file.close();
}

void MoveStructure::build() {
#if TALLY_MODE
    tally_checkpoints = movi_options->get_tally_checkpoints();
#endif
    std::string bwt_filename = movi_options->get_ref_file() + std::string(".bwt");
    std::ifstream bwt_file(bwt_filename);
    uint64_t end_pos = 0;
    if (movi_options->is_preprocessed()) {
        std::ifstream heads_file(bwt_filename + ".heads", std::ios::in | std::ios::binary);
        std::vector<char> heads_((std::istreambuf_iterator<char>(heads_file)), std::istreambuf_iterator<char>());
        std::cerr << "Number of BWT runs: " << heads_.size() << "\n";
        original_r  = heads_.size();
        std::ifstream len_file(bwt_filename + ".len", std::ios::in | std::ios::binary);
        for (uint64_t i = 0; i < original_r; i++) {
            if (i>0 && i % 100000 == 0)
                std::cerr << "original_r: " << i << "\r";
            size_t len = 0;
            len_file.read(reinterpret_cast<char*>(&len), 5);
            end_pos += len;
        }
    } else {
        bwt_file.clear();
        bwt_file.seekg(0, std::ios_base::end);
        std::streampos end_pos_ = bwt_file.tellg();
        if (movi_options->is_verbose())
            std::cerr << "end_pos: " << end_pos_ << "\n";
        bwt_file.seekg(0);
        end_pos = static_cast<uint64_t>(end_pos_);
    }
    std::cerr << "building.. \n";
    if (splitting) {
        std::string splitting_filename = movi_options->get_ref_file() + std::string(".d_col");
        std::ifstream splitting_file(splitting_filename);

        bits.load(splitting_file);
        std::cerr << "bits.size after loading the d_col file: " << bits.size() << "\n";
        rbits = sdsl::rank_support_v<>(&bits);
        std::cerr << "The main bit vector (bits) is loaded from the d_col file.\n";
    } else {
        bits = sdsl::bit_vector(static_cast<uint64_t>(end_pos) + 1, 0); // 5137858051
        bits[0] = 1;
        std::cerr << "The main bit vector (bits) is initialized.\n";
    }
    bwt_string = "";
    bwt_string.resize(static_cast<uint64_t>(end_pos) + 1);
    uint64_t all_chars_count = 256;
    alphamap.resize(all_chars_count);
    std::fill(alphamap.begin(), alphamap.end(), alphamap.size());
    std::vector<uint64_t> all_chars(all_chars_count, 0);
    uint64_t bwt_curr_length = 0;
    // if (!splitting)
    r = 1;
    original_r = 1;
    // TODO Use a size based on the input size

    uint64_t split_by_max_run = 0;
    uint64_t split_by_thresholds = 0;

    // A any mode with splitting (modes 1 and 4) does not work with the preprocessed build
    if (movi_options->is_preprocessed()) {
        length = static_cast<uint64_t>(end_pos);
        // We know only 4 characters (A, C, G, T) exists in the alphabet:
        char chars[4] = {'A', 'C', 'G', 'T'};
        for (size_t i = 0; i < 4; i++) {
            alphabet.push_back(chars[i]);
            alphamap[static_cast<size_t>(chars[i])] = i;
            sdsl::bit_vector* new_bit_vector = new sdsl::bit_vector(length, 0);
            occs.emplace_back(std::unique_ptr<sdsl::bit_vector>(new_bit_vector));
            std::cerr << i + 1 << " alphabet bitvectors " << (i==0 ? "is" : "are") << " built.\r";
        }
        std::cerr << "\n";

        std::ifstream len_file(bwt_filename + ".len", std::ios::in | std::ios::binary);
        std::ifstream heads_file(bwt_filename + ".heads", std::ios::in | std::ios::binary);
        std::vector<char> heads_((std::istreambuf_iterator<char>(heads_file)), std::istreambuf_iterator<char>());
        std::cerr << "Number of BWT runs: " << heads_.size() << "\n";
        /* heads_file.clear();
        heads_file.seekg(0, std::ios_base::end);
        std::streampos heads_end_pos = heads_file.tellg();
        std::cerr << "heads_end_pos:\t" << heads_end_pos << "\n"; */
        original_r  = heads_.size();
        std::vector<size_t> lens;
        // std::vector<char> heads;
        for (uint64_t i = 0; i < original_r; i++) {
            if (i>0 && i % 100000 == 0)
                std::cerr << "original_r: " << i << "\r";
            size_t len = 0;
            len_file.read(reinterpret_cast<char*>(&len), 5);
#if SPLIT_MAX_RUN
            size_t remaining_length = len;
            while (remaining_length > MAX_RUN_LENGTH) {
                split_by_max_run += 1;
                lens.push_back(MAX_RUN_LENGTH);
                heads.push_back(heads_[i]);
                remaining_length -= MAX_RUN_LENGTH;
            }
            lens.push_back(remaining_length);
            heads.push_back(heads_[i]);
#endif
#if NO_SPLIT
            lens.push_back(len);
            heads.push_back(heads_[i]);
#endif
        }
        original_r  = heads.size();

        uint64_t all_bits_size = splitting ? rbits(length) + 1 : original_r + 1;
        r =  splitting ? rbits(length) : 1;

        std::cerr << "all_bits_size: " << all_bits_size << "\n";
        if (splitting)
            std::cerr << "rbits(length): " << rbits(length) << "\n";
        all_p.resize(all_bits_size);
        all_p[all_bits_size - 1] = length;
        bwt_curr_length = 0;
        for (uint64_t i = 0; i < original_r; i++) {
            if (!splitting) bits[bwt_curr_length] = 1;
            if (i>0 && i % 100000 == 0) {
                std::cerr << "original_r: " << i << "\t";
                std::cerr << "bwt_curr_length: " << bwt_curr_length << "\r";
            }
            // We know only 4 characters (A, C, G, T) exists in the alphabet:
            if (heads[i] != 'A' and heads[i] != 'C' and heads[i] != 'G' and heads[i] != 'T' and heads[i] != END_CHARACTER) {
                std::cerr << "The preprocessed file includes non A/C/G/T character.\n";
                exit(0);
            }
            size_t len = lens[i];
            /* size_t len = 0;
            len_file.read(reinterpret_cast<char*>(&len), 5); */
            all_p[i] = bwt_curr_length;
            if (heads[i] == END_CHARACTER) {
                end_bwt_idx = i;

                bwt_string[bwt_curr_length] = heads[i];
                bwt_curr_length++;
            } else {
                all_chars[static_cast<size_t>(heads[i])] += len;
                auto& bit_vec = *occs[alphamap[static_cast<uint64_t>(heads[i])]];
                for (size_t j = 0; j < len; j ++) {
                    /* if (splitting && bwt_curr_length > 0 && bits[bwt_curr_length]) {
                        r += 1;
                    } */
                    bwt_string[bwt_curr_length] = heads[i];
                    bit_vec[bwt_curr_length] = 1;
                    bwt_curr_length++;
                }
            }
        }
        for (size_t i = 0; i < 4; i++) {
            counts.push_back(all_chars[static_cast<size_t>(alphabet[i])]);
        }
        if (!splitting) r = original_r;
    } else {
        if (movi_options->is_thresholds() and (SPLIT_THRESHOLDS_TRUE))
            fill_bits_by_thresholds();
        // Reading the BWT from the file
        uint64_t current_char = bwt_file.get();
        uint16_t run_length = 0;
        while (current_char != EOF) { // && current_char != 10
            uint64_t current_char_ = static_cast<uint64_t>(current_char); // Is this line important?!
            run_length += 1;
            // if (current_char != 'A' and current_char != 'C' and current_char != 'G' and current_char != 'T')
            //    std::cerr << "\ncurrent_char:" << current_char << "---" << static_cast<uint64_t>(current_char) << "---\n";
            if (original_r % 100000 == 0) {
                std::cerr << "r: " << r << "\t";
                std::cerr << "original_r: " << original_r << "\t";
                std::cerr << "bwt_curr_length: " << bwt_curr_length << "\r";
            }
#if SPLIT_MAX_RUN
            // The first row is already set and accounted for, so we skip
            if (bwt_curr_length > 0 && current_char != bwt_string[bwt_curr_length - 1]) {
                // 1) A new run is detected if the next character is different
                original_r += 1;
                r += 1;
                run_length = 0;
                bits[bwt_curr_length] = 1;
            } else if (movi_options->is_thresholds() and bwt_curr_length > 0 and bits[bwt_curr_length] == 1) {
                // 2) A new run is detected if there is a non-trivial threshold at the next offset
                // The bit was already set by one of the threshold values
                // So, we have found a new run, and reset the run length
                r += 1;
                run_length = 0;
                split_by_thresholds += 1;
            } else if (run_length == MAX_RUN_LENGTH) {
                // 3) A new run is detected if the length of the run is greater than MAX_RUN_LENGTH
                r += 1;
                run_length = 0;
                split_by_max_run += 1;
                bits[bwt_curr_length] = 1;
            }
#endif
#if SPLIT_THRESHOLDS_FALSE
            if (bwt_curr_length > 0 && current_char != bwt_string[bwt_curr_length - 1]) {
                // 1) A new run is detected if the next character is different
                original_r += 1;
                r += 1;
                run_length = 0;
                if (splitting and !bits[bwt_curr_length]) {
                    std::cerr << "There is something wrong with the splitting vector.\n";
                    std::cerr << "The run boundaries should have been set to 1 since a new character was detected.\n";
                    exit(0);
                }
                bits[bwt_curr_length] = 1;
            } else if (splitting && bwt_curr_length > 0 && bits[bwt_curr_length]) {
                // 2) A new run is detected based on Nishimoto-Tabei splitting
                r += 1;
                run_length = 0;
            } else if (run_length == MAX_RUN_LENGTH) {
                // 3) A new run is detected if the length of the run is greater than MAX_RUN_LENGTH
                split_by_max_run += 1;
                r += 1;
                run_length = 0;
                bits[bwt_curr_length] = 1;
            }
#endif
            bwt_string[bwt_curr_length] = current_char;
            bwt_curr_length++;
            all_chars[current_char] += 1;

            current_char = bwt_file.get();
        }
// #if SPLIT_THRESHOLDS_FALSE
// r is laready set to r' which might be larger than r because of the length splitting
// so, we don't have to set it back to original_r anymore (the next line is commented)
//         if (!splitting) r = original_r;
// #endif
        length = bwt_curr_length; // bwt_string.length();
        std::cerr << "\n";

        // Building the auxilary structures
        uint64_t alphabet_index = 0;
        // END_CHARACTER in the bwt created by pfp is 0
        for (uint64_t i = 1; i < all_chars_count; i++) {
            if (all_chars[i] != 0) {
                auto current_char = static_cast<unsigned char>(i);
                if (movi_options->is_verbose())
                    std::cerr << "i is " << i << "\t" << current_char
                            << "\t" << all_chars[i] << " alphabet_index: " << alphabet_index << "\n";

                alphabet.push_back(current_char);
                counts.push_back(all_chars[i]);
                alphamap[i] = alphabet_index;
                alphabet_index += 1;

                sdsl::bit_vector* new_bit_vector = new sdsl::bit_vector(length, 0);
                occs.emplace_back(std::unique_ptr<sdsl::bit_vector>(new_bit_vector));
            }
        }
        std::cerr << "All the characters are indexed.\n";

        all_p.resize(r + 1);
        all_p[0] = 0;
        heads.resize(r);
        heads[0] = bwt_string[0];

        uint64_t r_idx = 0;
        for (uint64_t i = 0; i < length; i++) {
            if (i % 10000 == 0)
                std::cerr <<"length processed: " << i << "\r";
            if (i == length - 1 or bwt_string[i] != bwt_string[i+1] or bits[i+1]) {
                all_p[r_idx + 1] = i + 1;
                heads[r_idx + 1] = bwt_string[i + 1];
                r_idx += 1;
            }
            if (static_cast<uint64_t>(bwt_string[i]) == END_CHARACTER) {
                end_bwt_idx = r_idx - 1;
                continue;
            }
            auto& bit_vec = *occs[alphamap[static_cast<uint64_t>(bwt_string[i])]];
            bit_vec[i] = 1;
        }
        std::cerr << "\nAll the Occ bit vectors are built.\n";
    }
    // We don't need the original BWT anymore, the run-head characters are stored in heads.
    bwt_string.clear();
    bwt_string.shrink_to_fit();

    std::cerr << "\nsplit_by_max_run: " << split_by_max_run << "\n";
    std::cerr << "split_by_thresholds: " << split_by_thresholds << "\n";
    std::cerr << "n: " << length << "\n";
    std::cerr << "r: " << r << "\n";
    std::cerr << "original_r: " << original_r << "\n";
    rlbwt.resize(r);

    if (movi_options->is_verbose() and bits.size() < 1000)
        std::cerr << "bits: " << bits << "\n";
    if (!splitting) // The rank vector is already built if it was the splitting mode
        rbits = sdsl::rank_support_v<>(&bits);

    if (alphabet.size() > 4) {
        std::cerr << "Warning: There are more than 4 characters, the index expexts only A, C, T and G in the reference.\n";
    }

    for (auto& occ: occs) {
        std::cerr << occs_rank.size() << "\n";
        if (movi_options->is_verbose() and (*occ).size() < 1000)
            std::cerr << *occ << "\n";
        occs_rank.emplace_back(std::unique_ptr<sdsl::rank_support_v<> >(new sdsl::rank_support_v<>(occ.get())));
    }
    if (movi_options->is_verbose()) {
        std::cerr << "size occs_rank:" << occs_rank.size() << "\n";
        std::cerr << "All Occ rank vectors are built.\n";
    }


    // Building the move structure rows with O(r) loop
    uint64_t offset = 0;
    uint64_t max_len = 0;
    if (movi_options->is_verbose()) {
        std::cerr << "bits.size(): " << bits.size() << "\n";
        std::cerr << "rank_support_v<>(&bits)(bits.size()): " << sdsl::rank_support_v<>(&bits)(bits.size()) << "\n";
    }

#if TALLY_MODE
    tally_ids.resize(alphabet.size());
    uint64_t tally_ids_rows_count = r / tally_checkpoints + 2;
    std::vector<uint64_t> current_tally_ids;
    for (uint32_t alphabet_ind = 0; alphabet_ind < alphabet.size(); alphabet_ind++) {
        current_tally_ids.push_back(r);
        tally_ids[alphabet_ind].resize(tally_ids_rows_count);
    }
#endif

    //if (splitting)
    //    sbits = sdsl::select_support_mcl<>(&bits);
    std::vector<uint64_t> raw_ids;
    raw_ids.resize(r);
    std::cerr << r << "\t" << rlbwt.size() << "\t" << all_p.size() << "\n";
    for (uint64_t r_idx = 0; r_idx < r; r_idx++) {
        if (r_idx % 10000 == 0)
            std::cerr << r_idx << "\r";
        uint64_t lf  = 0;
        if (r_idx != end_bwt_idx) {
            uint64_t alphabet_index = alphamap[static_cast<uint64_t>(heads[r_idx])];
            lf = LF(all_p[r_idx], alphabet_index);
        }
        else
            lf = 0;
        uint64_t pp_id = rbits(lf) - 1;
        if (bits[lf] == 1)
            pp_id += 1;
        uint64_t len = all_p[r_idx + 1] - all_p[r_idx];
        // check the boundaries before performing select
        if (pp_id >= r) {
            std::cerr << "pp_id: " << pp_id << " r: " << r << " r_idx: " << r_idx << " lf: " << lf << "\n";
            exit(0); // TODO: add error handling
        }
        //TODO: Can we really not use sbits?
        // sbits(pp_id + 1) is replaced by all_p[pp_id]
        // if (lf < sbits(pp_id + 1)) {
        if (lf < all_p[pp_id]) {
            // std::cerr << lf << " " << sbits(pp_id + 1);
            std::cerr << pp_id << " " << lf << " " << all_p[pp_id] << "\n";
            exit(0); // TODO: add error handling
        }
        // if (splitting)
        //     offset = lf - sbits(pp_id + 1);
        // else
        offset = lf - all_p[pp_id];
        /* if (sbits(pp_id + 1) != all_p[pp_id])
            exit(0); */

        if (movi_options->is_verbose() and r_idx == 0) // or any run to be inspected
            std::cerr << "r_idx: " << r_idx
                        << " len: " << len
                        << " lf: " << lf
                        << " offset: " << offset
                        << " pp_id: " << pp_id
                        /*<< " sbits(pp_id): " << all_p[pp_id - 1]
                        << " sbits(pp_id + 1): " << all_p[pp_id]
                        << " sbits(pp_id - 1): " << all_p[pp_id - 2]*/
                        << "\n";

#if TALLY_MODE
        rlbwt[r_idx].init(len, offset);

        if (r_idx != end_bwt_idx) {
            // Only update the tally corresponding to the current run's character
            uint16_t char_index = alphamap[heads[r_idx]];
            // The first tally id might have been set to the wrong value since no run with that character was observed
            // So, we fix the values once we the first run with that character
            if (current_tally_ids[char_index] == r) {
                uint64_t tally_ids_index = r_idx / tally_checkpoints;
                for (int tid = 0; tid <= tally_ids_index; tid++) {
                    tally_ids[char_index][tid].set_value(pp_id);
                }
            }
            current_tally_ids[char_index] = pp_id;
        }

        if (r_idx % tally_checkpoints == 0) {
            uint64_t tally_ids_index = r_idx / tally_checkpoints;
            // We reached a check_point, store the tally ids for all the characters
            for (uint32_t alphabet_ind = 0; alphabet_ind < alphabet.size(); alphabet_ind++) {
                tally_ids[alphabet_ind][tally_ids_index].set_value(current_tally_ids[alphabet_ind]);
            }
        }
#else
        // rlbwt[r_idx].init(bwt_row, len, lf, offset, pp_id);
        rlbwt[r_idx].init(len, offset, pp_id);
#endif

        raw_ids[r_idx] = pp_id;
        // To take care of cases where length of the run
        // does not fit in uint16_t
#if SPLIT_MAX_RUN
        if (offset > MAX_RUN_LENGTH or len > MAX_RUN_LENGTH) {
            // Should not get here in the regular mode: MODE = 3
            std::cerr << "The length or the offset are too large.\n";
            std::cerr << "offset: " << offset << "\tlength: " << length << "\n";
            exit(0);
        }
#endif
#if SPLIT_THRESHOLDS_FALSE
        if (len > MAX_RUN_LENGTH) {
            std::cerr << "This shouldn't happen anymore, because the run splitting should be applied in every mode now.\n";
            exit(0);
            n_overflow.push_back(len);
            if (n_overflow.size() - 1 >= MAX_RUN_LENGTH) {
                std::cerr << "Warning: the number of runs with overflow n is beyond " << MAX_RUN_LENGTH<< "! " << n_overflow.size() - 1 << "\n";
                exit(0);
            }
            rlbwt[r_idx].set_n(n_overflow.size() - 1);
            rlbwt[r_idx].set_overflow_n();
        }
        if (offset > MAX_RUN_LENGTH) {
            std::cerr << "This shouldn't happen anymore, because the run splitting should be applied in every mode now.\n";
            exit(0);
            offset_overflow.push_back(offset);
            if (offset_overflow.size() - 1 >= MAX_RUN_LENGTH) {
                std::cerr << "Warning: the number of runs with overflow offset is beyond " << MAX_RUN_LENGTH<< "! " << offset_overflow.size() - 1 << "\n";
                exit(0);
            }
            rlbwt[r_idx].set_offset(offset_overflow.size() - 1);
            rlbwt[r_idx].set_overflow_offset();
        }
#endif
        if (len > max_len)
            max_len = len;
        if (movi_options->is_logs()) {
            if (run_lengths.find(len) != run_lengths.end())
                run_lengths[len] += 1;
            else
                run_lengths[len] = 1;
        }
        // rlbwt[r_idx].set_c(bwt_string[all_p[r_idx]], alphamap);
        rlbwt[r_idx].set_c(heads[r_idx], alphamap);
        // bit1_after_eof = alphamap[bwt_string[i+1]];
    }
#if TALLY_MODE
    // Set the last tally_id using the current_tally_ids
    std::cout << "\n";
    for (uint32_t alphabet_ind = 0; alphabet_ind < alphabet.size(); alphabet_ind++) {
        tally_ids[alphabet_ind][tally_ids[alphabet_ind].size()-1].set_value(current_tally_ids[alphabet_ind]);
    }
#endif

    std::cerr << "All the move rows are built.\n";
    std::cerr << "Max run length: " << max_len << "\n";

#if USE_THRESHOLDS
    // compute the thresholds
    compute_thresholds();
#endif

    first_runs.push_back(0);
    first_offsets.push_back(0);
    last_runs.push_back(0);
    last_offsets.push_back(0);
    uint64_t char_count = 1;
    for (uint64_t i = 0; i < counts.size(); i++) {
        uint64_t last_run = last_runs.back();
        uint64_t last_offset = last_offsets.back();
        if (last_offset + 1 >= rlbwt[last_run].get_n()) {
            first_runs.push_back(last_run + 1);
            first_offsets.push_back(0);
        } else {
            first_runs.push_back(last_run);
            first_offsets.push_back(last_offset + 1);
        }
        char_count += counts[i];
        // char_count points to the position of the first row with character i in the BWT
        if (movi_options->is_verbose())
            std::cerr << i << " char_count: " << char_count << "\n";
        auto occ_rank = rbits(char_count);
        // occ_rank is the number of set bits until and including the last row with character i - 1
        if (movi_options->is_verbose())
            std::cerr << i << " occ_rank: " << occ_rank << "\n";
        // The index of the run is number of set bits - 1 (because we count from 0):
        last_runs.push_back(static_cast<uint64_t>(occ_rank - 1));
        if (movi_options->is_verbose())
            std::cerr << char_count << " " << all_p[last_runs.back()] << " " << get_n(last_runs.back()) <<  "\n";
        last_offsets.push_back(char_count - all_p[last_runs.back()] - 1);
    }
    if (movi_options->is_verbose()) {
        for (uint64_t i = 0; i < first_runs.size(); i++) {
            std::cerr << "<--- " << first_runs[i] << "\t" << first_offsets[i] << "\n";
            std::cerr << "    -\n    -\n    -\n";
            std::cerr << ">--- " << last_runs[i] << "\t" << last_offsets[i] << "\n";
        }
    }

#if BLOCKED_MODE
    uint64_t max_raw_id = 0;
    uint64_t max_blocked_id = 0;
    std::vector<uint32_t> block_start_id;
    block_start_id.resize(alphabet.size(), 0);
    uint64_t block_count = 0;
    id_blocks.resize(alphabet.size());
    uint32_t max_diff = 0;
    for (uint64_t i = 0; i < rlbwt.size(); i++) {
        if (i % 10000 == 0)
            std::cerr << "i: " << i << "\r";
        uint64_t block_boundary = (BLOCK_SIZE) * block_count; // BLOCK_SIZE = 1 << 22 - 1
        if (i >= block_boundary) {
            std::cout << "\n\n" << block_count << "\t";
            // A new block is being initiated
            block_count += 1;
            // Store the largest (last) observed id for each character in the last block as a check point for the new block
            for (uint64_t j = 0; j < alphabet.size(); j++) {
                id_blocks[j].push_back(block_start_id[j]);
                std::cout << id_blocks[j][block_count - 1] << "\t";
                if (block_count > 1)
                    max_diff = std::max(max_diff, id_blocks[j][block_count - 1] - id_blocks[j][block_count - 2] );
            }
            std::cout << "\n";
        }
        if (i != end_bwt_idx) {
            uint64_t id = raw_ids[i];
            max_raw_id = std::max(id, max_raw_id);

            uint64_t block_number = i / BLOCK_SIZE;
            if (block_number != block_count - 1) {
                std::cerr << "The block calculation is incorrect.\n";
                std::cerr << "block_count: " << block_count << " block_number: " << block_number << "\n";
                exit(0);
            }

            // First Calcuate the id with respect to the first run with the same character -> adjusted_id
            // The first entry of the first _runs stores the ids for the global end run, so +1 is required
            uint64_t adjusted_id = id - first_runs[rlbwt[i].get_c() + 1];
            // Calcuated the distance of the ajustedt_id from the check point of that character
            uint64_t blocked_id = adjusted_id - static_cast<uint64_t>(id_blocks[rlbwt[i].get_c()][block_number]);
            if (blocked_id > MAX_BLOCKED_ID) {
                std::cerr << "The number of bits in the runs are not enough for storing the blocked_id.\n";
                std::cerr << "adjusted_id: " << adjusted_id << " id: " << id << " first_runs[rlbwt[i].get_c() + 1]: " << first_runs[rlbwt[i].get_c() + 1] << "\n";
                std::cerr << "id_blocks[rlbwt[i].get_c()][block_number]: " << id_blocks[rlbwt[i].get_c()][block_number] << "\n";
                std::cerr << "rlbwt[i].get_c(): " << static_cast<uint32_t>(rlbwt[i].get_c()) << "\n";
                std::cerr << "blocked_id: " << blocked_id << " MAX_BLOCKED_ID: " << MAX_BLOCKED_ID << "\n";
                std::cerr << "block_count: " << block_count << " block_number: " << block_number << "\n";
                std::cerr << "i: " << i << " BLOCK_SIZE: " << BLOCK_SIZE << "\n";
                exit(0);
            }
            rlbwt[i].set_id(blocked_id);
            if (rlbwt[i].get_id() != blocked_id) {
                std::cerr << rlbwt[i].get_id() << "-------" << blocked_id << "\n";
                exit(0);
            }

            max_blocked_id = std::max(blocked_id, max_blocked_id);

            // always holds the last id seen for each character
            block_start_id[rlbwt[i].get_c()] = static_cast<uint32_t>(id - first_runs[rlbwt[i].get_c() + 1]);
            if (id - first_runs[rlbwt[i].get_c() + 1] > std::numeric_limits<uint32_t>::max()) {
                std::cerr << "id - first_runs[rlbwt[i].get_c() + 1]: " << id - first_runs[rlbwt[i].get_c() + 1] << "\n";
                std::cerr << "The block_start_id does not fit in uint32_t\n";
                exit(0);
            }
        }
    }
    std::cerr << "max raw id: " << max_raw_id << "\t max blocked id: " << max_blocked_id << "\n";
    std::cerr << "max allowed blocked id: " << MAX_BLOCKED_ID << "\n";
    std::cerr << "Maximum distance between the check points: " << max_diff << "\n";
#endif

#if USE_NEXT_POINTERS
    if (constant) {
        std::cerr << "Computing the next ups and downs.\n";
        compute_nexts();
    }
#endif
    std::cerr << "The move structure building is done.\n";
}

void MoveStructure::fill_bits_by_thresholds() {
    for (int i = 0; i < thresholds.size(); i++) {
        bits[thresholds[i]] = 1;
    }
    std::cerr << "The bits vector is updated by thresholds.\n";
}

void MoveStructure::compute_run_lcs() {
    std::cout << "run,length,lcs\n";
    for (int i = 0; i < r; i++) {
        if (i % 1000 == 0) std::cerr << i << "\r";
        if (get_n(i) > 1) {
            uint64_t id_top = i;
            uint64_t id_bottom = i;
            uint64_t offset_top = 0;
            uint64_t offset_bottom = get_n(i) - 1;

            uint64_t lcs = 0;
            char c_top = get_char(id_top);
            char c_bottom = get_char(id_bottom);
            while (c_top == c_bottom and lcs <= 32) {
                lcs += 1;
                LF_move(offset_top, id_top);
                LF_move(offset_bottom, id_bottom);
                c_top = get_char(id_top);
                c_bottom = get_char(id_bottom);
            }
            std::cout << i << "," << get_n(i) << "," << lcs << "\n";
        } else {
            std::cout << i << "," << 1 << "," << -1 << "\n";
        }
    }
}

void MoveStructure::compute_ftab() {
    size_t ftab_k = movi_options->get_ftab_k();

    uint64_t ftab_size = std::pow(4, ftab_k);
    std::cerr << "ftab_size: " << ftab_size*sizeof(MoveInterval)*std::pow(10, -6) << " MB \n";
    ftab.clear();
    ftab.resize(ftab_size);
    std::cerr << "Number of ftab entries (ftab-k=" << ftab_k << "): " << ftab_size << "\n";
    for (uint64_t i = 0; i < ftab_size; i++) {
        std::string kmer = number_to_kmer(i, 2*(ftab_k), alphabet, alphamap);
        uint64_t kmer_code = kmer_to_number(ftab_k, kmer, 0, alphamap);
        if (kmer_code != i) {
            std::cerr << kmer << " " << i << " " << kmer_code << "\n";
        }
        int32_t pos_on_kmer = ftab_k - 1;
        MoveInterval interval(
            first_runs[alphamap[kmer[pos_on_kmer]] + 1],
            first_offsets[alphamap[kmer[pos_on_kmer]] + 1],
            last_runs[alphamap[kmer[pos_on_kmer]] + 1],
            last_offsets[alphamap[kmer[pos_on_kmer]] + 1]
        );
        MoveInterval kmer_interval = backward_search(kmer, pos_on_kmer, interval, std::numeric_limits<int32_t>::max());
        uint64_t match_count = kmer_interval.count(rlbwt);
        if (match_count >= 0 and pos_on_kmer == 0) {
            // std::cerr << kmer << " " << match_count << " " << kmer_interval << "\n";
            // if (!(ftab[i] == kmer_interval)) {
            //     std::cerr << "ftab[i] is different from kmer_interval: " << i << "\n";
            //     std::cerr << ftab[i] << "\n" << kmer_interval << "\n";
            //     exit(0);
            // }
            ftab[i] = kmer_interval;
        } else {
            // std::cerr << kmer << " not found!\n";
            // if (!ftab[i].is_empty()) {
            //     std::cerr << "ftab[i] is non-empty and different from kmer_interval: " << i << "\n";
            //     std::cerr << ftab[i].is_empty() << "\n" << match_count << " " << (pos_on_kmer == 0) << "\n";
            //     exit(0);
            // }
            ftab[i].make_empty();
        }
        // if (pos_on_kmer != 0) pos_on_kmer += 1;
        // std::cerr << kmer << "\n";
        // std::cerr << kmer.length() - pos_on_kmer << "/" << kmer.length() << "\t" << match_count << "\n";
    }
}

void MoveStructure::write_ftab() {
    size_t ftab_k = movi_options->get_ftab_k();

    uint64_t ftab_size = std::pow(4, ftab_k);
    if (ftab_size != ftab.size()) {
        std::cerr << "The size of the ftab is not correct: " << ftab_size << " != " << ftab.size() << "\n";
        exit(0);
    }
    std::string fname = movi_options->get_index_dir() + "/ftab." + std::to_string(ftab_k) + ".bin";
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    fout.write(reinterpret_cast<char*>(&ftab_k), sizeof(ftab_k));
    fout.write(reinterpret_cast<char*>(&ftab_size), sizeof(ftab_size));
    fout.write(reinterpret_cast<char*>(&ftab[0]), ftab_size*sizeof(ftab[0]));
    fout.close();
}

void MoveStructure::read_ftab() {
    size_t ftab_k = movi_options->get_ftab_k();
    if (movi_options->is_multi_ftab()) {
        ftabs.resize(ftab_k);
        while (ftab_k > 1) {
            std::string fname = movi_options->get_index_dir() + "/ftab." + std::to_string(ftab_k) + ".bin";
            std::ifstream fin(fname, std::ios::in | std::ios::binary);
            fin.read(reinterpret_cast<char*>(&ftab_k), sizeof(ftab_k));

            uint64_t ftab_size = 0;
            fin.read(reinterpret_cast<char*>(&ftab_size), sizeof(ftab_size));
            if (ftab_size != std::pow(4, ftab_k)) {
                std::cerr << "The size of the ftab is not correct: " << ftab_size << " != " << std::pow(4, ftab_k) << "\n";
                exit(0);
            }
            std::vector<MoveInterval> new_ftab;
            new_ftab.resize(ftab_size);
            fin.read(reinterpret_cast<char*>(&new_ftab[0]), ftab_size*sizeof(new_ftab[0]));
            fin.close();
            ftabs[ftab_k - 1] = new_ftab;
            ftab_k -= 1;
        }
    } else {
        std::string fname = movi_options->get_index_dir() + "/ftab." + std::to_string(ftab_k) + ".bin";
        std::ifstream fin(fname, std::ios::in | std::ios::binary);
        fin.read(reinterpret_cast<char*>(&ftab_k), sizeof(ftab_k));
        uint64_t ftab_size = 0;
        fin.read(reinterpret_cast<char*>(&ftab_size), sizeof(ftab_size));
        std::cerr << "ftab_k: " << ftab_k << "\t" << "ftab_size: " << ftab_size << "\n";
        if (ftab_size != std::pow(4, ftab_k)) {
            std::cerr << "The size of the ftab is not correct: " << ftab_size << " != " << std::pow(4, ftab_k) << "\n";
            exit(0);
        }
        ftab.resize(ftab_size);
        fin.read(reinterpret_cast<char*>(&ftab[0]), ftab_size*sizeof(ftab[0]));
        fin.close();
    }
}

#if USE_THRESHOLDS
void MoveStructure::compute_thresholds() {
    // initialize the start threshold at the last row
    std::vector<uint64_t> alphabet_thresholds(alphabet.size(), length);
    uint64_t thr_i = original_r - 1;
    uint64_t run_p = 0;
    if (movi_options->is_verbose()) {
        std::cerr << "thresholds.size():" << thresholds.size() << " length: "
                  << length << " r: " << r <<  " original_r: " << original_r << "\n";
        std::cerr << "thresholds[r]: " << thresholds[original_r-1] << " r-1: "
                  << thresholds[original_r - 2] << " r-2: " << thresholds[original_r - 3] << "\n";
    }
    for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
        if (i % 10000 == 0)
            std::cerr << "i: " << i << "\r";
        char rlbwt_c = alphabet[rlbwt[i].get_c()];
        /* if (thr_i != i) {
            std::cerr << "thr_i: " << thr_i << " i: " << i << "\n";
            exit(0);
        } */
        if (movi_options->is_verbose() and i >= rlbwt.size() - 10)
            std::cerr << "i: " << i << "\n"
                << "rlbwt[i].get_offset(): " << get_offset(i) << "\n "
                << "get_n(i): " << get_n(i) << "\n"
                << "thresholds[thr_i]: " << thresholds[thr_i] << " "
                << "rlbwt_c: " << rlbwt_c << "\n";
        std::vector<uint64_t> current_thresholds;
        current_thresholds.resize(alphabet.size() - 1);
        for (uint64_t j = 0; j < alphabet.size(); j++) {
            if (alphabet[j] == rlbwt_c) {
                if (thr_i >= thresholds.size()) {
                    std::cerr << " thr_i = " << thr_i << " is out of bound:\n";
                    std::cerr << " thresholds.size = " << thresholds.size() << "\n";
                    exit(0); // TODO: add error handling
                }
                alphabet_thresholds[j] = thresholds[thr_i];
            } else {
                if (alphamap_3[alphamap[rlbwt_c]][j] >= alphabet.size() - 1) {
                    std::cerr << "alphamap_3 is not working in general:\n"
                                << "alphabet.size() - 1 = " << alphabet.size() - 1 << "\n"
                                << "alphamap_3[alphamap[rlbwt_c]][j] = " << alphamap_3[alphamap[rlbwt_c]][j] << "\n";
                    exit(0); // TODO: add error handling
                }
                if (alphabet_thresholds[j] >= all_p[i] + get_n(i)) {
                    // rlbwt[i].thresholds[j] = get_n(i);
                    if (i == end_bwt_idx) {
                        end_bwt_idx_thresholds[j] = get_n(i);
                        std::cerr << "condition 1: end_bwt_idx_thresholds[" << j << "]:" << end_bwt_idx_thresholds[j] << "\n";
                        continue;
                    }
#if SPLIT_THRESHOLDS_FALSE
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], get_n(i));
#endif
#if SPLIT_THRESHOLDS_TRUE
                    rlbwt[i].set_threshold(alphamap_3[alphamap[rlbwt_c]][j], 1);
#endif
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = get_n(i);
                } else if (alphabet_thresholds[j] <= all_p[i]) {
                    if (i == end_bwt_idx) {
                        end_bwt_idx_thresholds[j] = 0;
                        std::cerr << "condition 2: end_bwt_idx_thresholds[" << j << "]:" << end_bwt_idx_thresholds[j] << "\n";
                        continue;
                    }
#if SPLIT_THRESHOLDS_FALSE
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], 0);
#endif
#if SPLIT_THRESHOLDS_TRUE
                    rlbwt[i].set_threshold(alphamap_3[alphamap[rlbwt_c]][j], 0);
#endif
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = 0;
                } else {
                    if (i == end_bwt_idx) {
                        end_bwt_idx_thresholds[j] = alphabet_thresholds[j] - all_p[i];
                        std::cerr << "condition 3: end_bwt_idx_thresholds[" << j << "]:" << end_bwt_idx_thresholds[j] << "\n";
                        continue;
                    }
#if SPLIT_THRESHOLDS_FALSE
                    if (alphabet_thresholds[j] - all_p[i] >= std::numeric_limits<uint16_t>::max()) {
                        rlbwt[i].set_overflow_thresholds();
                    }
                    set_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j], alphabet_thresholds[j] - all_p[i]);
#endif
#if SPLIT_THRESHOLDS_TRUE
                    std::cerr << "j: " << j << " i:" << i << " n:" << get_n(i) << " atj:"
                              << alphabet_thresholds[j] << " api:" << all_p[i] << "\n";
                    std::cerr << "This should never happen, since the runs are split at threshold boundaries.\n";
                    exit(0);
#endif
                    current_thresholds[alphamap_3[alphamap[rlbwt_c]][j]] = alphabet_thresholds[j] - all_p[i];
                }
                // printing the values for last 10 runs to debug
                if (movi_options->is_verbose() and i >= rlbwt.size() - 10) {
                    std::cerr << "\t j: \t" << j << " "
                        << "alphabet[j]: " << alphabet[j] << "  "
                        << "alphamap_3[alphamap[rlbwt_c]][j]: " << alphamap_3[alphamap[rlbwt_c]][j] << " "
                        << "alphabet_thresholds[j]: " << alphabet_thresholds[j] << " "
#if SPLIT_THRESHOLDS_FALSE
                        << "rlbwt[i].thresholds[j]:" << get_rlbwt_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j]) << "\n";
#endif
#if SPLIT_THRESHOLDS_TRUE
                        << "rlbwt[i].thresholds[j]:" << get_thresholds(i, alphamap_3[alphamap[rlbwt_c]][j]) << "\n";
#endif
                }
            }
        }
        if (i > 0 && (rlbwt[i].get_c() != rlbwt[i - 1].get_c()  || i == end_bwt_idx || i-1 == end_bwt_idx)) {
            thr_i--;
        }
#if SPLIT_THRESHOLDS_FALSE
        if (rlbwt[i].is_overflow_thresholds()) {
            if (thresholds_overflow.size() >= std::numeric_limits<uint16_t>::max()) {
                std::cerr << "Undefined behaviour: the number of runs with overflow thresholds is beyond uint16_t!"
                            << thresholds_overflow.size() << "\n";
                exit(0);
            }
            for (uint64_t k = 0; k < alphabet.size() - 1; k++) {
                set_rlbwt_thresholds(i, k, thresholds_overflow.size());
            }
            thresholds_overflow.push_back(current_thresholds);
        }
#endif
        run_p += get_n(i);
    } // end of the main for
    // since the thresholds for the first run was not calculated in the for
    for (uint64_t j = 0; j < alphabet.size() - 1; j++) {
#if SPLIT_THRESHOLDS_FALSE
        set_rlbwt_thresholds(0, j, 0);
#endif
#if SPLIT_THRESHOLDS_TRUE
        rlbwt[0].set_threshold(j, 0);
#endif
    }
}
#endif

uint64_t scan_count;
#if USE_NEXT_POINTERS
void MoveStructure::compute_nexts() {
    for (uint64_t i = rlbwt.size() - 1; i > 0; --i) {
        if (i % 100000 == 0)
            std::cerr << i << "\r";

        char rlbwt_c = alphabet[rlbwt[i].get_c()];
        for (uint64_t j = 0; j < alphabet.size(); j++) {
            if (i == end_bwt_idx) {
                auto idx = reposition_up(i, alphabet[j], scan_count);
                end_bwt_idx_next_up[j] = (idx == r) ? std::numeric_limits<uint16_t>::max() : i - idx;
                idx = reposition_down(i, alphabet[j], scan_count);
                end_bwt_idx_next_down[j] = (idx == r) ? std::numeric_limits<uint16_t>::max() : idx - i;
                continue;
            }
            if (alphabet[j] != rlbwt_c) {
                auto alphabet_idx = alphamap_3[alphamap[rlbwt_c]][j];
                auto idx = reposition_up(i, alphabet[j], scan_count);
                if (idx == r) {
                    rlbwt[i].set_next_up(alphabet_idx, std::numeric_limits<uint16_t>::max());
                } else {
                    if (i - idx > std::numeric_limits<uint16_t>::max())
                        std::cerr << "Warning - reposition up " << i - idx << " does not fit in 16 bits.\n";
                    rlbwt[i].set_next_up(alphabet_idx, i - idx);
                }

                idx = reposition_down(i, alphabet[j], scan_count);
                if (idx == r) {
                    rlbwt[i].set_next_down(alphabet_idx, std::numeric_limits<uint16_t>::max());
                } else {
                    if (idx - i > std::numeric_limits<uint16_t>::max())
                        std::cerr << "Warning - reposition down " << idx - i << " does not fit in 16 bits.\n";
                    rlbwt[i].set_next_down(alphabet_idx, idx - i);
                }
            }
        }
    }
}
#endif

uint64_t MoveStructure::fast_forward(uint64_t& offset, uint64_t idx, uint64_t x) {
    uint64_t idx_ = idx;
    /* if (movi_options->is_verbose()) {
        std::cerr << "\t \t fast forwarding:\n";
        std::cerr << " \t \t idx: " << idx << " offset: " << offset << " n:" << get_n(idx) << "\n";
    } */
    while (idx < r - 1 && offset >= get_n(idx)) {
        offset -= get_n(idx);
        idx += 1;
        // if (movi_options->is_verbose()) std::cerr << "\t \t ff offset based: +" << idx - idx_ << "\n";
    }
    /* if (movi_options->is_verbose())
        std::cerr << " \t \t idx: " << idx << " offset: " << offset << " n:" << get_n(idx) << "\n"; */
    return idx - idx_;
}

uint64_t MoveStructure::reposition_up(uint64_t idx, char c, uint64_t& scan_count) {
    if (idx == 0)
        return r;
    char row_c = alphabet[rlbwt[idx].get_c()];

    while (idx > 0 and row_c != c) {
        scan_count += 1;
        idx -= 1;
        row_c = alphabet[rlbwt[idx].get_c()];
    }
    /* if (logs) {
        if (repositions.find(scan_count) != repositions.end())
            repositions[scan_count] += 1;
        else
            repositions[scan_count] = 1;
    } */
    /* if (movi_options->is_verbose())
        std::cerr << "\t \t \t \t idx after the while in the reposition up" << idx << "\n";*/
    return (row_c == c) ? idx : r;
}

uint64_t MoveStructure::reposition_down(uint64_t idx, char c, uint64_t& scan_count) {
    if (idx == r - 1)
        return r;
    char row_c = alphabet[rlbwt[idx].get_c()];

    while (idx < r - 1 && row_c != c) {
        scan_count += 1;
        idx += 1;
        row_c = alphabet[rlbwt[idx].get_c()];
    }
    /* if (logs) {
        if (repositions.find(scan_count) != repositions.end())
            repositions[scan_count] += 1;
        else
            repositions[scan_count] = 1;
    } */
    /*if (movi_options->is_verbose())
        std::cerr << "\t \t \t \t idx after the while in the reposition down: " << idx << " " << c << " " << row_c << "\n";*/
    return (row_c == c) ? idx : r;
}

void MoveStructure::update_interval(MoveInterval& interval, char next_char) {
    if (!check_alphabet(next_char)) {
        std::cerr << "This should not happen! The character should have been checked before.\n";
        exit(0);
    }
#if USE_NEXT_POINTERS
    // if (movi_options->is_debug())
    //     dbg << alphabet[rlbwt[interval.run_start].get_c()] << " " << alphabet[rlbwt[interval.run_end].get_c()] << " " << next_char << "\n";
    uint64_t read_alphabet_index = alphamap[static_cast<uint64_t>(next_char)];
    if ((interval.run_start <= interval.run_end) and (alphabet[rlbwt[interval.run_start].get_c()] != next_char)) {
        if (interval.run_start == 0) {
            // To check if this case ever happens. If not, we should get rid of this condition.
            // if (movi_options->is_debug())
            //     dbg << "run_start is 0 before updating the interval!\n";
            while ((interval.run_start <= interval.run_end) and (alphabet[rlbwt[interval.run_start].get_c()] != next_char)) {
                interval.run_start += 1;
                interval.offset_start = 0;
                if (interval.run_start >= r) {
                    break;
                }
            }
        } else {
            char rlbwt_char = alphabet[rlbwt[interval.run_start].get_c()];
            uint64_t alphabet_index = alphamap_3[alphamap[rlbwt_char]][read_alphabet_index];
            if (rlbwt[interval.run_start].get_next_down(alphabet_index) == std::numeric_limits<uint16_t>::max()) {
                interval.run_start = r;
            } else {
                uint64_t run_start_ = interval.run_start + rlbwt[interval.run_start].get_next_down(alphabet_index);
                interval.run_start = run_start_;
                interval.offset_start = 0;
            }
        }
    }
    if ((interval.run_end >= interval.run_start) and (alphabet[rlbwt[interval.run_end].get_c()] != next_char)) {
        char rlbwt_char = alphabet[rlbwt[interval.run_end].get_c()];
        uint64_t alphabet_index = alphamap_3[alphamap[rlbwt_char]][read_alphabet_index];
        if (rlbwt[interval.run_end].get_next_up(alphabet_index) == std::numeric_limits<uint16_t>::max()) {
            interval.run_end = r;
        } else {
            uint64_t run_end_ = interval.run_end - rlbwt[interval.run_end].get_next_up(alphabet_index);
            interval.run_end = run_end_;
            interval.offset_end = rlbwt[interval.run_end].get_n() - 1;
        }
    }
#else
    while (interval.run_start <= interval.run_end and get_char(interval.run_start) != next_char) { //  >= or >
        interval.run_start += 1;
        interval.offset_start = 0;
        if (interval.run_start >= r) {
            break;
        }
    }
    while (interval.run_end >= interval.run_start and get_char(interval.run_end) != next_char) { //  >= or >
        interval.run_end -= 1;
        interval.offset_end = rlbwt[interval.run_end].get_n() - 1;
        if (interval.run_end == 0) {
            break;
        }
    }
#endif

}

bool MoveStructure::extend_bidirectional(char c_, MoveInterval& fw_interval, MoveInterval& rc_interval) {
    MoveInterval fw_interval_before_extension = fw_interval;
    char c_comp = complement(c_);

    bool res = backward_search_step(c_, fw_interval);
    if (res) {
        // The alphabet is already checked to be legal (ACGT)
        uint64_t skip = 0;
        uint64_t current_run = fw_interval_before_extension.run_start;
        uint64_t current_offset = fw_interval_before_extension.offset_start;
        while (current_run <= fw_interval_before_extension.run_end ) {
            if (current_run != end_bwt_idx) {
                if (complement(get_char(current_run)) < c_comp) {
                    uint64_t char_count = current_run != fw_interval_before_extension.run_end ?
                                        get_n(current_run) - current_offset : fw_interval_before_extension.offset_end - current_offset + 1;
                    skip += char_count;
                }
            } else {
                skip += 1;
            }
            current_run += 1;
            current_offset = 0;
        }

        while (skip != 0) {
            int rows_after = get_n(rc_interval.run_start) - 1 - rc_interval.offset_start;
            if (rows_after >= skip) {
                rc_interval.offset_start += skip;
                skip = 0;
            } else {
                rc_interval.run_start += 1;
                rc_interval.offset_start = 0;
                skip -= rows_after + 1;
            }
        }
        // Compute the run end for the rc interval
        skip = fw_interval.count(rlbwt) - 1;
        rc_interval.run_end = rc_interval.run_start;
        rc_interval.offset_end = rc_interval.offset_start;
        while (skip != 0) {
            int rows_after = get_n(rc_interval.run_end) - 1 - rc_interval.offset_end;
            if (rows_after >= skip) {
                rc_interval.offset_end += skip;
                skip = 0;
            } else {
                rc_interval.run_end += 1;
                rc_interval.offset_end = 0;
                skip -= rows_after + 1;
            }
        }
        return true;
    } else {
        return false;
    }
}

bool MoveStructure::extend_left(char c, MoveBiInterval& bi_interval) {
    char c_ = c;
    if (extend_bidirectional(c_, bi_interval.fw_interval, bi_interval.rc_interval)) {
        bi_interval.match_len += 1;
        return true;
    } else {
        return false;
    }
}

bool MoveStructure::extend_right(char c, MoveBiInterval& bi_interval) {
    char c_ = complement(c);
    if (extend_bidirectional(c_, bi_interval.rc_interval, bi_interval.fw_interval)) {
        bi_interval.match_len += 1;
        return true;
    } else {
        return false;
    }
}

MoveBiInterval MoveStructure::backward_search_bidirectional(std::string& R, int32_t& pos_on_r, MoveBiInterval interval, int32_t max_length) {
    // If the pattern is found, the pos_on_r will be equal to 0 and the interval will be non-empty
    // Otherwise the interval corresponding to match from the end until and including the updated pos_on_r will be returned
    // The input interval is non-empty and corresponds to the interval that matches the read (R) at pos_on_r
    MoveBiInterval prev_interval = interval;
    int32_t pos_on_r_saved = pos_on_r;
    while (pos_on_r > 0 and !interval.fw_interval.is_empty()) {
        /*if (!check_alphabet(R[pos_on_r - 1])) {
            return interval;
        }*/
        prev_interval = interval;
        bool res = extend_left(R[pos_on_r - 1], interval);
        if (!interval.fw_interval.is_empty()) {
            pos_on_r -= 1;
        }

        // The following is only for the backward_search is called for "look ahead" for kmer skipping
        if (pos_on_r_saved - pos_on_r > max_length)
            break;
    }
    if (interval.fw_interval.is_empty()) {
        return prev_interval;
    } else {
        return interval;
    }
}

MoveInterval MoveStructure::backward_search(std::string& R, int32_t& pos_on_r, MoveInterval interval, int32_t max_length) {
    // If the pattern is found, the pos_on_r will be equal to 0 and the interval will be non-empty
    // Otherwise the interval corresponding to match from the end until and including the updated pos_on_r will be returned
    // The input interval is non-empty and corresponds to the interval that matches the read (R) at pos_on_r
    MoveInterval prev_interval = interval;
    // uint64_t match_len = 1; --- unused variable
    int32_t pos_on_r_saved = pos_on_r;
    while (pos_on_r > 0 and !interval.is_empty()) {
        /*if (!check_alphabet(R[pos_on_r - 1])) {
            return interval;
        }*/
        prev_interval = interval;
        backward_search_step(R, pos_on_r, interval);
        if (!interval.is_empty()) {
            pos_on_r -= 1;
        }
        /*update_interval(interval, R[pos_on_r - 1]);
        if (!interval.is_empty()) {
            LF_move(interval.offset_start, interval.run_start);
            LF_move(interval.offset_end, interval.run_end);
            pos_on_r -= 1;
        }*/

        // The following is only for the backward_search is called for "look ahead" for kmer skipping
        if (pos_on_r_saved - pos_on_r > max_length)
            break;
    }
    if (interval.is_empty()) {
        return prev_interval;
    } else {
        return interval;
    }
}

MoveInterval MoveStructure::try_ftab(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, size_t ftab_k, bool rc) {
    auto& query_seq = mq.query();
    if (ftab_k > 1 and pos_on_r >= ftab_k - 1) {
        // uint64_t kmer_code = kmer_to_number(ftab_k, query_seq.substr(pos_on_r - ftab_k + 1, ftab_k), alphamap);
        uint64_t kmer_code = kmer_to_number(ftab_k, query_seq, pos_on_r - ftab_k + 1, alphamap, rc);
        if (kmer_code != std::numeric_limits<uint64_t>::max()) {
            auto& current_ftab = movi_options->is_multi_ftab() ? ftabs[ftab_k - 1] : ftab;
            if (!current_ftab[kmer_code].is_empty()) {
                // Add the skipped matching length, e.g., for zml computation
                if (movi_options->is_zml()) {
                    for (size_t i = 0; i < ftab_k - 1; i++) {
                        mq.add_ml(i, movi_options->is_stdout());
                    }
                }
                if (movi_options->is_zml() or movi_options->is_kmer()) {
                    match_len = ftab_k - 1;
                }
                pos_on_r = pos_on_r - ftab_k + 1;
                return current_ftab[kmer_code];
            }
        }
    }
    // If we reach here, we know the ftab could not be used, so we return an empty interval
    MoveInterval empty_interval;
    empty_interval.make_empty();
    return empty_interval;
}

MoveBiInterval MoveStructure::initialize_bidirectional_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len) {
    MoveBiInterval bi_interval;
    int32_t pos_on_r_before = pos_on_r;
    bi_interval.fw_interval = initialize_backward_search(mq, pos_on_r, match_len);
    // If the initialization was unsuccessfull, match_len of the bidirectional interval should be set to 0
    // Without the following condition, the match_len will be increamented in the following line
    if (match_len == 0) {
        bi_interval.match_len = match_len;
        return bi_interval;
    }
    // This is needed because of the way match_len is off by one
    match_len += 1;
    bi_interval.match_len = match_len;

    int32_t pos_on_r_rc = pos_on_r_before;
    uint64_t match_len_rc = 0;
    bi_interval.rc_interval = initialize_backward_search(mq, pos_on_r_rc, match_len_rc, true);

    // match_len_rc will be equal to match_len if both fw and rc are present int he reference
    if (match_len - 1 != match_len_rc) {
        std::cerr << "The reverse complement might not be present in the reference.\n";
        exit(0);
    }

    return bi_interval;
}

MoveInterval MoveStructure::initialize_backward_search(MoveQuery& mq, int32_t& pos_on_r, uint64_t& match_len, bool rc) {
    all_initializations += 1;
    // Initialize assuming that the character at pos_on_r exists in the alphabet
    size_t ftab_k = movi_options->get_ftab_k();
    if (movi_options->is_multi_ftab()) {
        while (ftab_k > 1 and pos_on_r >= ftab_k - 1) {
            //std::cerr << ftab_k << "\n";
            MoveInterval ftab_res = try_ftab(mq, pos_on_r, match_len, ftab_k, rc);
            if (!ftab_res.is_empty())
                return ftab_res;
            ftab_k -= 2;
        }
    } else if (ftab_k > 1) {
        MoveInterval ftab_res = try_ftab(mq, pos_on_r, match_len, ftab_k, rc);
        if (!ftab_res.is_empty())
            return ftab_res;
    }
    // If we reach here, we know the ftab could not be used, so we initialize regularly
    if (movi_options->is_multi_ftab() and ftab_k < movi_options->get_ftab_k()) {
        ftab_k += 2;
    }
    if (pos_on_r >= ftab_k - 1) {
        no_ftab += 1;
    }
    auto& query_seq = mq.query();
    auto first_char_index = alphamap[rc ? complement(query_seq[pos_on_r]): query_seq[pos_on_r]] + 1;
    MoveInterval interval(
        first_runs[first_char_index],
        first_offsets[first_char_index],
        last_runs[first_char_index],
        last_offsets[first_char_index]
    );
    return interval;
}

bool MoveStructure::backward_search_step(char c, MoveInterval& interval) {
    if (!check_alphabet(c)) {
        interval.make_empty();
        return false;
    }

    update_interval(interval, c);
    if (!interval.is_empty()) {
        LF_move(interval.offset_start, interval.run_start);
        LF_move(interval.offset_end, interval.run_end);
        return true;
    } else {
        return false;
    }
}

uint64_t MoveStructure::backward_search_step(std::string& R, int32_t& pos_on_r, MoveInterval& interval) {
    // It is assumed that the interval represents a match until and including the position pos_on_r
    // Then we try to see if the match can be extended to position pos_on_r - 1
    // The interval becomes empty if the match cannot be extended
    // otherwise it is updated according to the character on pos_on_r - 1
    uint64_t ff_count = 0;
    if (pos_on_r <= 0) {
        std::cerr << "The backward search step never be called on position 0 of the read.\n";
        exit(0);
    }

    if (!check_alphabet(R[pos_on_r - 1])) {
        interval.make_empty();
        return ff_count;
    }

    update_interval(interval, R[pos_on_r - 1]);
    if (!interval.is_empty()) {
        ff_count += LF_move(interval.offset_start, interval.run_start);
        ff_count += LF_move(interval.offset_end, interval.run_end);
    }

    return ff_count;
}

uint64_t MoveStructure::query_zml(MoveQuery& mq) {
    auto& query_seq = mq.query();
    int32_t pos_on_r = query_seq.length() - 1;
    uint64_t match_len = 0;
    uint64_t ff_count_tot = 0;

    while (!check_alphabet(query_seq[pos_on_r]) and pos_on_r >= 0) {
        mq.add_ml(0, movi_options->is_stdout());
        pos_on_r -= 1;
    }

    if (pos_on_r < 0) {
        // Special case where no character in the read exists in the index.
        return 0;
    }

    // Multi-class classification
    if (movi_options->is_multi_classify()) {
        for (uint16_t i = 0; i < num_species; i++) {
            classify_cnts[i] = 0;
        }
    }

    MoveInterval interval = initialize_backward_search(mq, pos_on_r, match_len);
    while (pos_on_r > 0) {
        MoveInterval prev_interval = interval;
        ff_count_tot += backward_search_step(query_seq, pos_on_r, interval);
        if (!interval.is_empty()) {
            mq.add_ml(match_len, movi_options->is_stdout());
            pos_on_r -= 1;
            match_len += 1;
        } else {
            // Classification based on maximal matching
            if (movi_options->is_multi_classify() && match_len >= movi_options->get_min_match_len()) {
                if (prev_interval.run_start == prev_interval.run_end) {
                    for (uint64_t i = prev_interval.offset_start; i <= prev_interval.offset_end; i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_start] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }   
                } else {
                    for (uint64_t i = prev_interval.offset_start; i < get_n(prev_interval.run_start); i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_start] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }    
                    for (uint64_t r_ind = prev_interval.run_start + 1; r_ind < prev_interval.run_end; r_ind++) {
                        for (uint64_t i = 0; i < get_n(r_ind); i++) {
                            uint64_t full_ind = run_offsets[r_ind] + i;
                            uint16_t cur_doc = doc_pats[full_ind];
                            classify_cnts[cur_doc] += match_len;
                        }  
                    } 
                    for (uint64_t i = 0; i <= prev_interval.offset_end; i++) {
                        uint64_t full_ind = run_offsets[prev_interval.run_end] + i;
                        uint16_t cur_doc = doc_pats[full_ind];
                        classify_cnts[cur_doc] += match_len;
                    }           
                }
            }

            mq.add_ml(match_len, movi_options->is_stdout());
            pos_on_r -= 1;
            match_len = 0;
            while (!check_alphabet(query_seq[pos_on_r]) and pos_on_r > 0) {
                mq.add_ml(match_len, movi_options->is_stdout());
                pos_on_r -= 1;
            }
            // Special case where the character at position 0 of the read does not exist in the index.
            if (check_alphabet(query_seq[pos_on_r]))
                interval = initialize_backward_search(mq, pos_on_r, match_len);
        }
    }
    if (interval.is_empty()) {
        match_len = 0;
    }
    mq.add_ml(match_len, movi_options->is_stdout());

    // Document occuring the most is the genotype we think the query is from.
    if (movi_options->is_multi_classify()) {
        uint16_t best_doc = 0;
        for (uint16_t i = 1; i < num_species; i++) {
            if (classify_cnts[i] > classify_cnts[best_doc]) {
                best_doc = i;
            }
        }
        
        // Document occuring the most is the genotype we think the query is from.
        out_file << to_taxon_id[best_doc] << " ";
        for (uint16_t i = 0; i < num_species; i++) {
        //    out_file << classify_cnts[i] << " ";
        }
        out_file << "\n";
    }

    return ff_count_tot;
}

uint64_t MoveStructure::query_backward_search(MoveQuery& mq, int32_t& pos_on_r) {
    auto& query_seq = mq.query();
    // Check the special case of non-existing character at the end of the read
    // before initializing the interval based on that character
    if (!check_alphabet(query_seq[pos_on_r])) {
        pos_on_r += 1; // even the character at the end_pos was not found.
        return 0;
    }
    // Initial the interval by matching the character at the end of the read (pos_on_r)
    uint64_t not_used = 0;
    MoveInterval initial_interval = initialize_backward_search(mq, pos_on_r, not_used);
    return backward_search(query_seq, pos_on_r, initial_interval, std::numeric_limits<int32_t>::max()).count(rlbwt);
}

bool MoveStructure::look_ahead_ftab(MoveQuery& mq, uint32_t pos_on_r, int32_t& step) {
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();
    // int32_t pos_on_r_ahead = pos_on_r - static_cast<int32_t>(k/2);
    for (step = 0; step <= 19 ; step += 1) {
        int32_t pos_on_r_ahead = pos_on_r - k + ftab_k + step;
        uint64_t kmer_code = kmer_to_number(ftab_k, query_seq, pos_on_r_ahead - ftab_k, alphamap);
        if (kmer_code != std::numeric_limits<uint64_t>::max() and !ftab[kmer_code].is_empty()) {
            // return true;
        } else {
            return false;
        }
    }
    return true;
}

bool MoveStructure::look_ahead_backward_search(MoveQuery& mq, uint32_t pos_on_r, int32_t step) {
    size_t ftab_k = movi_options->get_ftab_k();
    size_t k = movi_options->get_k();
    auto& query_seq = mq.query();

    uint64_t match_len = 0;
    int32_t pos_on_r_ahead = pos_on_r - step;
    MoveInterval initial_interval = initialize_backward_search(mq, pos_on_r_ahead, match_len);
    backward_search(query_seq, pos_on_r_ahead, initial_interval, k - step - match_len);
    if (pos_on_r - pos_on_r_ahead >= k - 1) {
        return true;
    } else {
        return false;
    }
}

void MoveStructure::flat_and_serialize_colors_vectors() {
    std::unordered_map<uint32_t, uint64_t> flat_index;
    for (uint64_t i = 0; i < unique_doc_sets.size(); i++) {
        if (i % 100000 == 0) {
            std::cerr << "i: " << i << "\r";
        }
        flat_index[i] = flat_colors.size();
        flat_colors.push_back(unique_doc_sets[i].size());
        for (uint64_t j = 0; j < unique_doc_sets[i].size(); j++) {
            flat_colors.push_back(unique_doc_sets[i][j]);
        }
    }
    std::cerr << "\n";

    doc_set_flat_inds.resize(doc_set_inds.size());
    for (uint64_t i = 0; i < doc_set_inds.size(); i++) {
        if (i % 100000 == 0) {
            std::cerr << "i: " << i << "\r";
        }
        doc_set_flat_inds[i].set_value(flat_index[doc_set_inds[i]]);
    }
    std::cerr << "\n";

    std::cerr << "flat_colors.size(): " << flat_colors.size() << "\n";

    std::string fname = movi_options->get_index_dir() + "/doc_sets_flat.bin";
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    if (!fout) {
        throw std::runtime_error("Failed to open the index file at: " + movi_options->get_index_dir());
    }
    uint64_t flat_colors_size = flat_colors.size();
    fout.write(reinterpret_cast<char*>(&flat_colors_size), sizeof(uint64_t));
    fout.write(reinterpret_cast<char*>(&flat_colors[0]), flat_colors.size() * sizeof(uint16_t));
    fout.write(reinterpret_cast<char*>(&doc_set_flat_inds[0]), doc_set_flat_inds.size() * sizeof(doc_set_flat_inds[0]));
    fout.close();
}

void MoveStructure::add_colors_to_rlbwt() {
#if MODE == 3 or MODE == 6
    rlbwt_colored.resize(rlbwt.size());
    for (uint64_t i = 0; i < rlbwt.size(); i++) {
        if (i % 100000 == 0) {
            std::cerr << "i: " << i << "\r";
        }
        rlbwt_colored[i].id = rlbwt[i].id;
        rlbwt_colored[i].n = rlbwt[i].n;
        rlbwt_colored[i].offset = rlbwt[i].offset;
        if (movi_options->is_compressed()) {
            if (doc_set_inds[i] >= unique_doc_sets.size()) {
                rlbwt_colored[i].color_id = std::numeric_limits<uint16_t>::max();
            } else {
                rlbwt_colored[i].color_id = static_cast<uint16_t>(doc_set_inds[i]);
            }
        } else {
            rlbwt_colored[i].color_id = doc_set_inds[i];
        }
    }
    std::cerr << "\n";
#endif
}
uint64_t MoveStructure::query_pml(MoveQuery& mq) {
    auto& R = mq.query();
    int32_t pos_on_r = R.length() - 1;
    uint64_t idx = r - 1; // or we can start from a random position in the rlbwt std::rand() % r
    uint64_t offset = get_n(idx) - 1;

    uint64_t match_len = 0;
    uint16_t ff_count = 0;
    uint64_t ff_count_tot = 0;
    uint64_t scan_count = 0;
    auto t1 = std::chrono::high_resolution_clock::now();

    if (movi_options->is_verbose()) {
        std::cerr << "beginning of the search \ton query: " << mq.query() << "\t";
        std::cerr << "and on BWT, idx(r-1): " << idx << " offset: " << offset << "\n";
    }

    // Multi-class classification
    if (movi_options->is_multi_classify()) {
        for (uint16_t i = 0; i < num_species; i++) {
            if (!movi_options->is_pvalue_scoring()) {
                classify_cnts[i] = 0;
            } else {
                doc_scores[i] = 0;
            }
        }
    }

    uint16_t best_doc = std::numeric_limits<uint16_t>::max(); // for multi-class classification
    uint16_t second_best_doc = std::numeric_limits<uint16_t>::max();
    uint64_t iteration_count = 0;
    uint32_t sum_matching_lengths = 0;
    while (pos_on_r > -1) {
        iteration_count += 1;
        if (movi_options->is_logs() and (iteration_count-1)%200 == 0) {
            t1 = std::chrono::high_resolution_clock::now();
        }

        if (movi_options->is_verbose())
            std::cerr << "Searching position " << pos_on_r << " of the read:\n";

        auto& row = rlbwt[idx];
        uint64_t row_idx = idx;
        char row_c = alphabet[row.get_c()];
        if (!check_alphabet(R[pos_on_r])) {
            // The character from the read does not exist in the reference
            match_len = 0;
            scan_count = 0;

            if (movi_options->is_verbose())
                std::cerr << "\t The character " << R[pos_on_r] << " does not exist.\n";
        } else if (row_c == R[pos_on_r]) {
            // Case 1
            match_len += 1;
            scan_count = 0;

            if (movi_options->is_verbose()) {
                std::cerr << "\t Cas1: It was a match. \n" << "\t Continue the search...\n";
                std::cerr << "\t match_len: " << match_len << "\n";
                std::cerr << "\t current_id: " << idx << "\t row.id: " << get_id(row_idx) << "\n"
                          << "\t row.get_n: " << get_n(row_idx) << " rlbwt[idx].get_n: " << get_n(get_id(row_idx)) << "\n"
                          << "\t offset: " << offset << "\t row.get_offset(): " << get_offset(row_idx) << "\n";
            }
        } else {
            // Case 2
            // Repositioning up or down (randomly or with thresholds)
            if (movi_options->is_verbose())
                std::cerr << "\t Case2: Not a match, looking for a match either up or down...\n";

            uint64_t idx_before_reposition = idx;
#if USE_THRESHOLDS
            bool up = movi_options->is_random_repositioning() ?
                               reposition_randomly(idx, offset, R[pos_on_r], scan_count) :
                               reposition_thresholds(idx, offset, R[pos_on_r], scan_count);
#else
            // When there is no threshold, reposition randomly
            bool up = reposition_randomly(idx, offset, R[pos_on_r], scan_count);
#endif

            match_len = 0;
            // scan_count = (!constant) ? std::abs((int)idx - (int)idx_before_reposition) : 0;

            char c = alphabet[rlbwt[idx].get_c()];

            if (movi_options->is_verbose())
                std::cerr << "\t up: " << up << " idx: " << idx << " c:" << c << "\n";

            // sanity check
            if (c == R[pos_on_r]) {
                // Observing a match after the reposition
                // The right match_len should be:
                // min(new_lcp, match_len + 1)
                // But we cannot compute lcp here
                offset = up ? get_n(idx) - 1 : 0;
                /* if (movi_options->is_verbose())
                    std::cerr << "\t idx: " << idx << " offset: " << offset << "\n"; */
            } else {
                std::cerr << "\t \t This should not happen!\n";
                std::cerr << "\t \t pos: " << pos_on_r << " r[pos]:" <<  R[pos_on_r] << " t[pointer]:" << c << "\n";
                std::cerr << "\t \t " << up << ", " << R[pos_on_r] << ", " << pos_on_r << "\n";
                std::cerr << "\t \t ";
                for (int k = 10; k > 0; --k)
                    std::cerr << alphabet[rlbwt[idx - k].get_c()] << "-";
                for (int k = 0; k < 10; k++)
                    std::cerr << alphabet[rlbwt[idx + k].get_c()] << "-";
                std::cerr << "\n";
                auto saved_idx = idx;

                movi_options->set_verbose(true);
#if USE_THRESHOLDS
                reposition_thresholds(saved_idx, offset, R[pos_on_r], scan_count);
#endif
                exit(0);
            }
        }
    
        sum_matching_lengths += match_len;
        mq.add_ml(match_len, movi_options->is_stdout());
        pos_on_r -= 1;

        // LF step
        ff_count = LF_move(offset, idx);
        ff_count_tot += ff_count;
        if (movi_options->is_logs()) {
            if (iteration_count % 200 == 0) {
                auto t2 = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
                mq.add_cost(elapsed);
            }
            mq.add_fastforward(ff_count);
            mq.add_scan(scan_count);
        }

        if (movi_options->is_multi_classify()) {
            if (match_len >= movi_options->get_min_match_len()) {
                /*uint64_t full_ind = run_offsets[idx] + offset;
                uint16_t cur_doc = doc_pats[full_ind];
                classify_cnts[cur_doc]++;
                if (classify_cnts[cur_doc] >= classify_cnts[best_doc]) {
                    best_doc = cur_doc;
                }*/ 

                uint64_t color_id;
#if COLOR_MODE == 1
                color_id = static_cast<uint64_t>(rlbwt[idx].color_id);
                // Skip doc sets that weren't saved (thrown away by compression).
                if (color_id >= unique_doc_sets.size()) continue;
#else
                if (movi_options->is_doc_sets_vector_of_vectors()) {
                    color_id = static_cast<uint64_t>(doc_set_inds[idx]);
                    // Skip doc sets that weren't saved (thrown away by compression).
                    if (color_id >= unique_doc_sets.size()) continue;
                } else {
                    color_id = doc_set_flat_inds[idx].get();
                    // Skip doc sets that weren't saved (thrown away by compression).
                    if (color_id >= flat_colors.size()) continue;
                }
#endif
                std::span<uint16_t> cur_set;
                if (movi_options->is_doc_sets_vector_of_vectors()) {
                    std::vector<uint16_t> &cur_set_vec = unique_doc_sets[color_id];
                    cur_set = std::span<uint16_t>(cur_set_vec.data(), cur_set_vec.size());
                } else {
                    uint32_t cur_set_size = flat_colors[color_id];
                    cur_set = std::span<uint16_t>(flat_colors.data() + color_id + 1,
                                                  flat_colors.data() + color_id + 1 + cur_set_size);
                }

                for (int doc : cur_set) {
                    if (!movi_options->is_pvalue_scoring()) {
                        classify_cnts[doc]++;
                        if (doc != best_doc) {
                            if (best_doc == std::numeric_limits<uint16_t>::max() || classify_cnts[doc] > classify_cnts[best_doc]) {
                                second_best_doc = best_doc;
                                best_doc = doc;
                            } else if (second_best_doc == std::numeric_limits<uint16_t>::max() || classify_cnts[doc] > classify_cnts[second_best_doc]) {
                                second_best_doc = doc;
                            }
                        }
                    } else {
                        // p value strategy
                        double val = match_len - (log_lens[doc] / log4);
                        if (val >= 0) {
                            doc_scores[doc] += std::min(val, 1.);
                            if (doc != best_doc) {
                                if (best_doc == std::numeric_limits<uint16_t>::max() || doc_scores[doc] > doc_scores[best_doc]) {
                                    second_best_doc = best_doc;
                                    best_doc = doc;
                                } else if (second_best_doc == std::numeric_limits<uint16_t>::max() || doc_scores[doc] > doc_scores[second_best_doc]) {
                                    second_best_doc = doc;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (movi_options->is_multi_classify()) {
        float PML_mean = static_cast<float>(sum_matching_lengths) / mq.query().length();
        if (PML_mean < UNCLASSIFIED_THRESHOLD || best_doc == std::numeric_limits<uint16_t>::max()) {
            // Not present
            out_file << "0,0\n";
        } else {
            if (second_best_doc == std::numeric_limits<uint16_t>::max()) {
                out_file << to_taxon_id[best_doc] << ",0";
            } else {
                float best_doc_cnt, second_best_doc_cnt, second_best_diff;
                if (movi_options->get_min_match_len() > 0) {
                    best_doc_cnt = classify_cnts[best_doc];
                    second_best_doc_cnt = classify_cnts[second_best_doc];
                    second_best_diff = (best_doc_cnt - second_best_doc_cnt);
                } else {
                    // p-value strategy
                    best_doc_cnt = doc_scores[best_doc];
                    second_best_doc_cnt = doc_scores[second_best_doc];
                    second_best_diff = (best_doc_cnt - second_best_doc_cnt);
                }

                if (second_best_diff < 0.05 * best_doc_cnt) {
                    out_file << to_taxon_id[best_doc] << "," << to_taxon_id[second_best_doc];
                } else {
                    out_file << to_taxon_id[best_doc] << ",0";
                }
            }
            out_file << "\n";
        }
    }

    return ff_count_tot;
}

#if USE_THRESHOLDS
bool MoveStructure::reposition_thresholds(uint64_t& idx, uint64_t offset, char r_char, uint64_t& scan_count) {
    // If offset is greather than or equal to the threshold, reposition down
    // otherwise, reposition up
    uint64_t saved_idx = idx;
    uint64_t alphabet_index = alphamap[static_cast<uint64_t>(r_char)];
    scan_count = 0;

    if (movi_options->is_verbose())
        std::cerr << "\t \t \t repositioning with thresholds ... \n";

    char rlbwt_char = alphabet[rlbwt[idx].get_c()];

    if (movi_options->is_verbose()) {
        std::cerr << "\t \t \t alphabet_index: " << alphabet_index << " r_char:" << r_char << " rlbwt_char:" << rlbwt_char << "\n";
        std::cerr << "\t \t \t idx:" << idx << "\n"
                  << "\t \t \t offset: " << offset << " threshold:" << get_thresholds(idx, alphamap_3[alphamap[rlbwt_char]][alphabet_index]) << "\n";
    }

    if (idx == end_bwt_idx) {
        if (movi_options->is_verbose()) std::cerr << "\t \t \t idx == end_bwt_idx"
                                << "\n\t \t \t idx: " << idx << " end_bwt_idx: " << end_bwt_idx << " alphabet_index: " << alphabet_index << "\n"
                                << "\t \t \t end_bwt_idx_thresholds[alphabet_index]: " << end_bwt_idx_thresholds[alphabet_index] << "\n";
        if (offset >= end_bwt_idx_thresholds[alphabet_index]) {
            if (movi_options->is_verbose())
                std::cerr << "\t \t \t Repositioning down with thresholds:\n";
#if USE_NEXT_POINTERS
            if (constant) {
                scan_count += 1;
                if (end_bwt_idx_next_down[alphabet_index] == std::numeric_limits<uint16_t>::max())
                    idx = r;
                else
                    idx = saved_idx + end_bwt_idx_next_down[alphabet_index];
            }
#endif

#if THRESHOLDS_WITHOUT_NEXTS
            idx = reposition_down(saved_idx, r_char, scan_count);
#endif
            if (r_char != alphabet[rlbwt[idx].get_c()])
                std::cerr << "1: " << r_char << " " << alphabet[rlbwt[idx].get_c()];
            return false;
        } else {
            if (movi_options->is_verbose())
                std::cerr << "\t \t \t Repositioning up with thresholds:\n";
#if USE_NEXT_POINTERS
            if (constant) {
                scan_count += 1;
                if (end_bwt_idx_next_up[alphabet_index] == std::numeric_limits<uint16_t>::max())
                    idx = r;
                else
                    idx = saved_idx - end_bwt_idx_next_up[alphabet_index];
            }
#endif

#if THRESHOLDS_WITHOUT_NEXTS
            idx = reposition_up(saved_idx, r_char, scan_count);
#endif
            if (r_char != alphabet[rlbwt[idx].get_c()])
                std::cerr << "2: " << r_char << " " << alphabet[rlbwt[idx].get_c()];
            return true;
        }
    }

    if (movi_options->is_verbose()) std::cerr << "\t \t \t rlbwt[idx].get_offset(): " << get_offset(idx)
                            << " get_thresholds(idx, alphabet_index): " << get_thresholds(idx, alphamap_3[alphamap[rlbwt_char]][alphabet_index]) 
                            << "\n\t \t \t idx:" << idx << "\n";

    alphabet_index = alphamap_3[alphamap[rlbwt_char]][alphabet_index];
    if (alphabet_index == 3)
        std::cerr << "error: alphamap_3 is not working in general - " 
                    << alphabet_index << "!\n";

    if (offset >= get_thresholds(idx, alphabet_index)) {
        if (movi_options->is_verbose())
            std::cerr << "\t \t \t Repositioning down with thresholds:\n";
#if USE_NEXT_POINTERS
        if (constant) {
            scan_count += 1;
            if (rlbwt[saved_idx].get_next_down(alphabet_index) == std::numeric_limits<uint16_t>::max())
                idx = r;
            else
                idx = saved_idx + rlbwt[saved_idx].get_next_down(alphabet_index);
        } else {
            std::cerr << "MODE is set to " << MODE <<", but the constant variable is false.\n";
            exit(0);
        }
#endif
        auto tmp = idx;
#if THRESHOLDS_WITHOUT_NEXTS
        idx = reposition_down(saved_idx, r_char, scan_count);
#endif
        if (r_char != alphabet[rlbwt[idx].get_c()]) {
            std::cerr << "3: " << r_char << " " << alphabet[rlbwt[idx].get_c()] << "\n";
            std::cerr << "idx: " << idx << " saved_idx: " << saved_idx << " tmp: " << tmp << "\n";
            std::cerr << "offset: " << offset << "\n";
            std::cerr << "get_thresholds(saved_idx, alphabet_index): " << get_thresholds(saved_idx, alphabet_index) << "\n";
        }
        return false;
    } else {
        if (movi_options->is_verbose())
            std::cerr << "\t \t \t Repositioning up with thresholds:\n";
#if USE_NEXT_POINTERS
        scan_count += 1;
        if (constant) {
            if (rlbwt[saved_idx].get_next_up(alphabet_index) == std::numeric_limits<uint16_t>::max())
                idx = r;
            else
                idx = saved_idx - rlbwt[saved_idx].get_next_up(alphabet_index);
        }
#endif

#if THRESHOLDS_WITHOUT_NEXTS
        idx = reposition_up(saved_idx, r_char, scan_count);
#endif
        if (r_char != alphabet[rlbwt[idx].get_c()]) {
            std::cerr << "idx: " << idx << " saved_idx: " << saved_idx << "\n";
            std::cerr << "4: " << r_char << " " << alphabet[rlbwt[idx].get_c()] << "\n";
            std::cerr << "offset: " << offset << "\n";
            std::cerr << "get_thresholds(saved_idx, alphabet_index): " << get_thresholds(saved_idx, alphabet_index) << "\n";
        }
        return true;
    }

    // TODO: default return?

    if (r_char != alphabet[rlbwt[idx].get_c()])
        std::cerr << "5: " << r_char << " " << alphabet[rlbwt[idx].get_c()];
    return false;
}
#endif

bool MoveStructure::reposition_randomly(uint64_t& idx, uint64_t& offset, char r_char, uint64_t& scan_count) {
    uint64_t saved_idx = idx;
    thread_local ThreadRandom random_generator;
    // uint16_t reposition_direction = random_generator.get_random() % 2;
    uint16_t reposition_direction =  offset * 2 < get_n(idx) ? 1 : 0;
    bool up = false;
    scan_count = 0;
    if (movi_options->is_verbose())
        std::cerr << "idx before repositioning: " << idx << "\n";


    // Selecting the right direction for trivial cases at the beginning and end of the move table
    if (idx == r - 1) {
        reposition_direction = 1;
    }
    if (idx == 0) {
        reposition_direction = 0;
    }

    if ( (reposition_direction == 1 && idx > 0) or idx == r - 1) {

        // repositioning up
        up = true;

        if (movi_options->is_verbose())
            std::cerr << "Repositioning up randomly:\n";

        idx = reposition_up(saved_idx, r_char, scan_count);
        if (movi_options->is_verbose())
            std::cerr << "idx after repositioning up: " << idx << "\n";

        if (idx >= r) {
            if (movi_options->is_verbose())
                std::cerr << "Up didn't work, try repositioning down:\n";

            // repositioning down
            up = false;
            idx = reposition_down(saved_idx, r_char, scan_count);
            if (movi_options->is_verbose())
                std::cerr << "idx after repositioning down: " << idx << "\n";
            if (idx == r) {
                // TODO
                std::cerr << "Neither up or down repositioning works.\n";
                std::cerr << "The character does not exist in the index.\n";
                exit(0);
            }
        }
    } else {

        // repositioning down
        up = false;

        if (movi_options->is_verbose())
            std::cerr << "Repositioning down randomly:\n";

        idx = reposition_down(saved_idx, r_char, scan_count);
        if (movi_options->is_verbose())
            std::cerr << "idx after repositioning down: " << idx << "\n";

        if (idx >= r) {
            if (movi_options->is_verbose())
                std::cerr << "Down didn't work, try repositioning up:\n";

            // repositioning up
            up = true;
            idx = reposition_up(saved_idx, r_char, scan_count);
            if (movi_options->is_verbose())
                std::cerr << "idx after repositioning up: " << idx << "\n";
            if (idx == r) {
                // TODO
                std::cerr << "Neither up or down repositioning works.\n";
                std::cerr << "The character does not exist in the index.\n";
                exit(0);
            }
        }
    }

    // sanity check
    char c = alphabet[rlbwt[idx].get_c()];
    if (c != r_char or idx == r) {
        if (movi_options->is_verbose()) {
            std::cerr << "c: " << c << "\n";
            std::cerr << "idx: " << idx << "\n";
        }
        // TODO
        std::cerr << "This should never happen.\n";
        exit(0);
    }

    return up;
}

bool MoveStructure::check_alphabet(char& c) {
    // for (char c_: alphabet) {
    //     if (c == c_)
    //         return true;
    // }
    // return false;
    if (movi_options->ignore_illegal_chars_status() > 0) {
        if (alphamap[static_cast<uint64_t>(c)] == alphamap.size()) {
            thread_local ThreadRandom random_generator;
            c = movi_options->ignore_illegal_chars_status() == 1 ?  'A' : alphabet[ random_generator.get_random() % alphabet.size() ];
            return true;
        }
    }
    return alphamap[static_cast<uint64_t>(c)] != alphamap.size();
}

void MoveStructure::serialize_doc_pats(std::string fname) {
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    std::cerr << "Writing doc pats to: " << fname << std::endl;

    fout.write(reinterpret_cast<char*>(&doc_pats[0]), length * sizeof(doc_pats[0]));
    fout.close();
}

void MoveStructure::deserialize_doc_pats(std::string fname) {
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    
    doc_pats.resize(length);
    fin.read(reinterpret_cast<char*>(&doc_pats[0]), length * sizeof(doc_pats[0]));
    fin.close();
    
    std::cerr << "Finished deserializing document patterns" << std::endl;
}

void MoveStructure::serialize_doc_sets(std::string fname) {
    std::ofstream fout(fname, std::ios::out | std::ios::binary);

    size_t unique_cnt = unique_doc_sets.size();
    std::cerr << "Number of unique document sets: " << unique_cnt << std::endl; 
    fout.write(reinterpret_cast<char*>(&unique_cnt), sizeof(unique_cnt));
    for (size_t i = 0; i < unique_doc_sets.size(); i++) {
        std::vector<uint16_t> &cur = unique_doc_sets[i];
        uint16_t doc_cnt = cur.size();
        fout.write(reinterpret_cast<char*>(&doc_cnt), sizeof(doc_cnt));
        fout.write(reinterpret_cast<char*>(&cur[0]), cur.size() * sizeof(cur[0]));
    }
    fout.write(reinterpret_cast<char*>(&doc_set_inds[0]), r * sizeof(doc_set_inds[0]));
    fout.close();
}

void MoveStructure::compute_color_ids_from_flat() {

    for (uint64_t color_offset = 0; color_offset < flat_colors.size(); ) {

        uint32_t color_id = color_offset_to_id.size();
        color_offset_to_id[color_offset] = color_id;

        uint32_t next_color =  color_offset + flat_colors[color_offset] + 1;
        color_offset = next_color;

    }
}

void MoveStructure::deserialize_doc_sets_flat() {
    std::string fname = movi_options->get_index_dir() + "/doc_sets_flat.bin";
    std::ifstream fin(fname, std::ios::in | std::ios::binary);

    uint64_t flat_colors_size = 0;
    fin.read(reinterpret_cast<char*>(&flat_colors_size), sizeof(uint64_t));
    std::cerr << "flat_colors_size: " << flat_colors_size << "\n";

    flat_colors.resize(flat_colors_size);
    fin.read(reinterpret_cast<char*>(&flat_colors[0]), flat_colors_size * sizeof(flat_colors[0]));
    std::cerr << "Read flat_colors\n";

    doc_set_flat_inds.resize(r);
    fin.read(reinterpret_cast<char*>(&doc_set_flat_inds[0]), r * sizeof(doc_set_flat_inds[0]));
    std::cerr << "Read doc_set_flat_inds\n";

    fin.close();
}

void MoveStructure::deserialize_doc_sets(std::string fname) {
    std::ifstream fin(fname, std::ios::in | std::ios::binary);

    size_t unique_cnt = 0;
    fin.read(reinterpret_cast<char*>(&unique_cnt), sizeof(unique_cnt));
    std::cerr << "Number of unique document sets: " << unique_cnt << std::endl; 
    unique_doc_sets.resize(unique_cnt);
    for (size_t i = 0; i < unique_cnt; i++) {
        std::vector<uint16_t> &cur = unique_doc_sets[i];
        uint16_t doc_cnt = 0;
        fin.read(reinterpret_cast<char*>(&doc_cnt), sizeof(doc_cnt));
        cur.resize(doc_cnt);
        fin.read(reinterpret_cast<char*>(&cur[0]), (size_t) doc_cnt * sizeof(cur[0]));
    }

#if COLOR_MODE == 0
    // The following are now stored in rlbwt_colored (colored move rows)
    doc_set_inds.resize(r);
    fin.read(reinterpret_cast<char*>(&doc_set_inds[0]), r * sizeof(doc_set_inds[0]));
#endif
    fin.close();

    /*uint64_t missing_cnt = 0;
    for (size_t i = 0; i < r; i++) {
        if (doc_set_inds[i] >= unique_doc_sets.size()) {
            missing_cnt++;
        }
    }
    std::cerr << "Fraction of runs without color: " << (double) missing_cnt / r << std::endl;*/
}

void MoveStructure::serialize() {
    mkdir(movi_options->get_index_dir().c_str(),0777);
    std::string fname = movi_options->get_index_dir() + "/index.movi";
    if (movi_options->is_color_move_rows()) {
        fname = movi_options->get_index_dir() + "/index_colored.movi";
    }
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    if (!fout) {
        throw std::runtime_error("Failed to open the index file at: " + movi_options->get_index_dir());
    }

    if (!movi_options->is_no_header()) {
        char index_type = static_cast<char>(MODE);
        fout.write(reinterpret_cast<char*>(&index_type), sizeof(index_type));
    }

    if (movi_options->is_verbose()) {
        std::cerr << "length: " << length << " r: " << r << " end_bwt_idx: " << end_bwt_idx << "\n";
    }
    fout.write(reinterpret_cast<char*>(&length), sizeof(length));
    fout.write(reinterpret_cast<char*>(&r), sizeof(r));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx), sizeof(end_bwt_idx));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_thresholds[0]), 4*sizeof(end_bwt_idx_thresholds[0]));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_next_down[0]), 4*sizeof(end_bwt_idx_next_down[0]));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_next_up[0]), 4*sizeof(end_bwt_idx_next_up[0]));

    uint64_t alphamap_size = alphamap.size();
    fout.write(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    fout.write(reinterpret_cast<char*>(&alphamap[0]), alphamap.size()*sizeof(alphamap[0]));

    uint64_t alphabet_size = alphabet.size();
    fout.write(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));    
    fout.write(reinterpret_cast<char*>(&alphabet[0]), alphabet.size()*sizeof(alphabet[0]));

    fout.write(reinterpret_cast<char*>(&splitting), sizeof(splitting));
    fout.write(reinterpret_cast<char*>(&constant), sizeof(constant));

    // This is not used any more as the onebit modes is deprecated
    fout.write(reinterpret_cast<char*>(&onebit), sizeof(onebit));

    if (movi_options->is_color_move_rows()) {
        fout.write(reinterpret_cast<char*>(&rlbwt_colored[0]), rlbwt_colored.size()*sizeof(rlbwt_colored[0]));
    } else {
        fout.write(reinterpret_cast<char*>(&rlbwt[0]), rlbwt.size()*sizeof(rlbwt[0]));
    }
#if TALLY_MODE
    fout.write(reinterpret_cast<char*>(&tally_checkpoints), sizeof(tally_checkpoints));
    uint64_t tally_ids_len = tally_ids[0].size();
    fout.write(reinterpret_cast<char*>(&tally_ids_len), sizeof(tally_ids_len)); 
    for (uint32_t i = 0; i < alphabet.size(); i++) {
        fout.write(reinterpret_cast<char*>(&tally_ids[i][0]), tally_ids[i].size()*sizeof(tally_ids[i][0]));
    }
#endif
    uint64_t n_overflow_size = n_overflow.size();
    fout.write(reinterpret_cast<char*>(&n_overflow_size), sizeof(n_overflow_size));
    fout.write(reinterpret_cast<char*>(&n_overflow[0]), n_overflow.size()*sizeof(uint64_t));
    uint64_t offset_overflow_size = offset_overflow.size();
    fout.write(reinterpret_cast<char*>(&offset_overflow_size), sizeof(offset_overflow_size));
    fout.write(reinterpret_cast<char*>(&offset_overflow[0]), offset_overflow.size()*sizeof(uint64_t));
    uint64_t thresholds_overflow_size = thresholds_overflow.size();
    fout.write(reinterpret_cast<char*>(&thresholds_overflow_size), sizeof(thresholds_overflow_size));
    for (uint32_t i = 0; i < thresholds_overflow_size; i++) {
        fout.write(reinterpret_cast<char*>(&thresholds_overflow[i][0]), (alphabet.size() - 1)*sizeof(thresholds_overflow[i][0]));
    }

    size_t orig_size = orig_string.size();
    fout.write(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
    fout.write(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));

    fout.write(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));

    uint64_t counts_size = counts.size();
    fout.write(reinterpret_cast<char*>(&counts_size), sizeof(counts_size));
    fout.write(reinterpret_cast<char*>(&counts[0]), counts.size()*sizeof(counts[0]));

    uint64_t last_runs_size = last_runs.size();
    fout.write(reinterpret_cast<char*>(&last_runs_size), sizeof(last_runs_size));
    fout.write(reinterpret_cast<char*>(&last_runs[0]), last_runs.size()*sizeof(last_runs[0]));
    fout.write(reinterpret_cast<char*>(&last_offsets[0]), last_offsets.size()*sizeof(last_offsets[0]));
    fout.write(reinterpret_cast<char*>(&first_runs[0]), first_runs.size()*sizeof(first_runs[0]));
    fout.write(reinterpret_cast<char*>(&first_offsets[0]), first_offsets.size()*sizeof(first_offsets[0]));

#if BLOCKED_MODE
    if (id_blocks.size() > 0) {
        uint64_t id_blocks_size = id_blocks[0].size();
        fout.write(reinterpret_cast<char*>(&id_blocks_size), sizeof(id_blocks_size));
        for (uint64_t i = 0; i < alphabet.size(); i++) {
            fout.write(reinterpret_cast<char*>(&id_blocks[i][0]), id_blocks_size*sizeof(id_blocks[i][0]));
        }
    } else {
        uint64_t id_blocks_size = 0;
        fout.write(reinterpret_cast<char*>(&id_blocks_size), sizeof(id_blocks_size));
    }
#endif
    fout.write(reinterpret_cast<char*>(&original_r), sizeof(original_r));

    fout.close();
}

void MoveStructure::deserialize() {
#if COLOR_MODE == 1
    std::string fname = movi_options->get_index_dir() + "/index_colored.movi";
#else
    std::string fname = movi_options->get_index_dir() + "/index.movi";
#endif
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    if (!fin) {
        // Attempt to read an index file built with the old index name
        fname = movi_options->get_index_dir() + "/movi_index.bin";
        fin.open(fname, std::ios::in | std::ios::binary);
        if (!fin) {
            throw std::runtime_error("Failed to open the index file at: " + movi_options->get_index_dir());
        }
    }
    std::cerr << "fname: " << fname << std::endl;
    fin.seekg(0, std::ios::beg); 

    if (!movi_options->is_no_header()) {
        char index_type;
        fin.read(reinterpret_cast<char*>(&index_type), sizeof(index_type));
        if (MODE != index_type) {
            throw std::runtime_error("MODE does not match the index_type in the header.");
        }
        std::cerr << "The " << program() << " index is being used.\n";
    }

    fin.read(reinterpret_cast<char*>(&length), sizeof(length));
    fin.read(reinterpret_cast<char*>(&r), sizeof(r));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx), sizeof(end_bwt_idx));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_thresholds[0]), 4*sizeof(end_bwt_idx_thresholds[0]));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_next_down[0]), 4*sizeof(end_bwt_idx_next_down[0]));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_next_up[0]), 4*sizeof(end_bwt_idx_next_up[0]));

    std::cerr << "n: " << length << "\nr: " << r << "\n" << "$: " << end_bwt_idx << "\n";

    uint64_t alphamap_size;
    fin.read(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    alphamap.resize(alphamap_size);
    fin.read(reinterpret_cast<char*>(&alphamap[0]), alphamap_size*sizeof(alphamap[0]));
    uint64_t alphabet_size;
    fin.read(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));
    alphabet.resize(alphabet_size);
    fin.read(reinterpret_cast<char*>(&alphabet[0]), alphabet_size*sizeof(alphabet[0]));
    std::cerr << "alphabet_size: " << alphabet_size << "\n";
    if (alphabet.size() > 4) {
        std::cerr << "Warning: There are more than 4 characters, the index expects only A, C, T and G in the reference.\n";
    }

    fin.read(reinterpret_cast<char*>(&splitting), sizeof(splitting));
    fin.read(reinterpret_cast<char*>(&constant), sizeof(constant));

    // This is not used any more as the onebit modes is deprecated
    fin.read(reinterpret_cast<char*>(&onebit), sizeof(onebit));

    rlbwt.resize(r);
    std::cerr << "sizeof(MoveRow): " << sizeof(MoveRow) << std::endl;
    fin.read(reinterpret_cast<char*>(&rlbwt[0]), r*sizeof(MoveRow));
#if TALLY_MODE
    fin.read(reinterpret_cast<char*>(&tally_checkpoints), sizeof(tally_checkpoints));
    std::cerr << "Tally mode with tally_checkpoints = " << tally_checkpoints << "\n";
    uint64_t tally_ids_len = 0;
    fin.read(reinterpret_cast<char*>(&tally_ids_len), sizeof(tally_ids_len));
    tally_ids.resize(alphabet.size());
    for (uint32_t i = 0; i < alphabet.size(); i++) {
        tally_ids[i].resize(tally_ids_len);
        fin.read(reinterpret_cast<char*>(&tally_ids[i][0]), tally_ids_len*sizeof(MoveTally));
    }
#endif

    uint64_t n_overflow_size;
    fin.read(reinterpret_cast<char*>(&n_overflow_size), sizeof(n_overflow_size));
    n_overflow.resize(n_overflow_size);
    fin.read(reinterpret_cast<char*>(&n_overflow[0]), n_overflow_size*sizeof(uint64_t));
    uint64_t offset_overflow_size;
    fin.read(reinterpret_cast<char*>(&offset_overflow_size), sizeof(offset_overflow_size));
    offset_overflow.resize(offset_overflow_size);
    fin.read(reinterpret_cast<char*>(&offset_overflow[0]), offset_overflow_size*sizeof(uint64_t));
    uint64_t thresholds_overflow_size;
    fin.read(reinterpret_cast<char*>(&thresholds_overflow_size), sizeof(thresholds_overflow_size));
    thresholds_overflow.resize(thresholds_overflow_size);
    for (uint32_t i = 0; i < thresholds_overflow_size; i++) {
        thresholds_overflow[i].resize(alphabet.size() - 1);
        fin.read(reinterpret_cast<char*>(&thresholds_overflow[i][0]), (alphabet.size() - 1)*sizeof(uint64_t));
    }

    std::cerr << "All the move rows are read.\n";

    size_t orig_size;
    fin.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));

    fin.read(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));
    reconstructed = false;

    fin.read(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));

    uint64_t counts_size = 0;
    fin.read(reinterpret_cast<char*>(&counts_size), sizeof(counts_size));
    counts.resize(counts_size);
    fin.read(reinterpret_cast<char*>(&counts[0]), counts_size*sizeof(uint64_t));

    uint64_t last_runs_size = 0;
    fin.read(reinterpret_cast<char*>(&last_runs_size), sizeof(last_runs_size));
    last_runs.resize(last_runs_size);
    fin.read(reinterpret_cast<char*>(&last_runs[0]), last_runs_size*sizeof(uint64_t));
    last_offsets.resize(last_runs_size);
    fin.read(reinterpret_cast<char*>(&last_offsets[0]), last_runs_size*sizeof(uint64_t));
    first_runs.resize(last_runs_size);
    fin.read(reinterpret_cast<char*>(&first_runs[0]), last_runs_size*sizeof(uint64_t));
    first_offsets.resize(last_runs_size);
    fin.read(reinterpret_cast<char*>(&first_offsets[0]), last_runs_size*sizeof(uint64_t));

#if BLOCKED_MODE
    uint64_t id_blocks_size = 0;
    fin.read(reinterpret_cast<char*>(&id_blocks_size), sizeof(id_blocks_size));
    std::cerr << "id_blocks_size: " << id_blocks_size << "\n";
    if (id_blocks_size > 0) {
        id_blocks.resize(alphabet.size());
        for (uint64_t i = 0; i < alphabet.size(); i++) {
            id_blocks[i].resize(id_blocks_size);
            fin.read(reinterpret_cast<char*>(&id_blocks[i][0]), id_blocks_size*sizeof(uint32_t));
        }
    }
#endif
    // To be able to load the indexes that haven't stored original_r
    if (fin.eof()) {
        std::cerr << "original_r was not stored in the index!\n";
    } else {
        fin.read(reinterpret_cast<char*>(&original_r), sizeof(original_r));
    }

    fin.close();
}

void MoveStructure::load_document_info() {
    // Read in document offsets.
    std::ifstream doc_offsets_file(movi_options->get_index_dir() + "/ref.fa.doc_offsets");
    if (!doc_offsets_file.good()) {
        std::cerr << "Error: doc_offsets file not found at \""
                  << movi_options->get_index_dir() << "/ref.fa.doc_offsets\"" << std::endl;
        exit(0);
    }
    uint64_t doc_offset;
    while ((doc_offsets_file >> doc_offset)) {
        doc_offsets.push_back(doc_offset);
    }
    doc_offsets_file.close();
    num_docs = doc_offsets.size();
    std::cerr << "num_docs: " << num_docs << std::endl;
    num_species = num_docs;

    // Read in document taxon id
    std::ifstream doc_ids_file(movi_options->get_index_dir() + "/ref.fa.doc_ids");
    if (doc_ids_file.good()) {
        uint32_t doc_id;
        while ((doc_ids_file >> doc_id)) {
            doc_ids.push_back(doc_id);
            taxon_id_compress[doc_id] = 0;
        }
        doc_ids_file.close();
    } else {
        // No doc ids information
        doc_ids.resize(num_docs);
        for (size_t i = 0; i < num_docs; i++) {
            doc_ids[i] = i + 1;
            taxon_id_compress[i + 1] = 0;
        }
    }

    // Compress taxon id to 0...(num_species - 1)
    num_species = 0;
    for (auto &[taxon_id, c_id] : taxon_id_compress) {
        to_taxon_id.push_back(taxon_id);
        c_id = num_species++;
    }
    for (size_t i = 0; i < doc_ids.size(); i++) {
        doc_ids[i] = taxon_id_compress[doc_ids[i]];
    }

    log_lens.resize(num_species);
    for (int i = 0; i < num_docs; i++) {
        uint64_t doc_len = doc_offsets[i] - (i == 0 ? 0 : doc_offsets[i - 1]);
        log_lens[doc_ids[i]] += doc_len;
    }
    for (int i = 0; i < num_species; i++) {
        log_lens[i] = log(log_lens[i]);
    }

    // Fill in powers of 2 array.
    pow2[0] = 1;
    for (size_t i = 1; i < ARR_SIZE; i++) {
        pow2[i] = (pow2[i - 1] << 1);
        if (pow2[i] >= MOD) {
            pow2[i] -= MOD;
        }
    }
}

void MoveStructure::initialize_classify_cnts() {
    classify_cnts.resize(num_species, 0);
    doc_scores.resize(num_species, 0);
}

void MoveStructure::verify_lfs() {
    uint64_t not_matched = 0;
    for (uint64_t i = 0; i < all_p.size(); i++) {
        std::uint64_t end_ = (i < all_p.size() - 1) ? all_p[i + 1] : length;
        for (uint64_t j = all_p[i]; j < end_; j++) {
            uint64_t offset_ = j - all_p[i];
            uint64_t idx_ = i;
            uint64_t lf = 0;
            if (i != end_bwt_idx) {
                uint64_t alphabet_index = rlbwt[idx_].get_c();
                lf = LF(j, alphabet_index);
            } else {
                std::cerr << "end_run = " << i << " len: " << rlbwt[i].get_n () << "\n";
            }
            LF_move(offset_, idx_);
            uint64_t lf_move = all_p[idx_] + offset_;
            if (lf != lf_move) {
                not_matched += 1;
                std::cerr << "j\t" << j << "\n";
                std::cerr << "idx\t" << i << "\n";
                std::cerr << "offset\t" << j - all_p[i] << "\n";
                std::cerr << "rlbwt[idx].get_id\t" << get_id(i) << "\n";
                std::cerr << "get_offset(i)\t" << get_offset(i) << "\n";
                for (uint64_t k = 0; k <= i; k++) {
                    std::cerr << rlbwt[k].get_n() << " ";
                }
                std::cerr << "\n\n";

                std::cerr << "lf\t" << lf << "\n";
                std::cerr << "lf_move\t" << lf_move << "\n";
                std::cerr << "idx_\t" << idx_ << "\n";
                std::cerr << "offset_\t" << offset_ << "\n";
                std::cerr << "all_p[idx_]\t" << all_p[idx_] << "\n";
                std::cerr << "\n\n\n";
            }
        }
    }
    if (not_matched == 0) {
        std::cerr << "All the LF_move operations are correct.\n";
    } else {
        original_r = 0;
        std::cerr << "There are " << not_matched << " LF_move operations that failed to match the true lf results.\n";
    }
}

void MoveStructure::analyze_rows() {
    for (int i = 0; i < first_runs.size(); i++) {
        std::cerr << i << "\t" << first_runs[i] << "\t" << last_runs[i] << "\n";
    }
    std::vector<uint64_t> counts_length(16,0);
    std::vector<uint64_t> counts_offset(16,0);
    std::vector<uint64_t> counts_threshold0(16,0);
    std::vector<uint64_t> counts_threshold1(16,0);
    std::vector<uint64_t> counts_threshold2(16,0);
    uint64_t counter = 0;
    uint64_t split_thresholds = 0;
    uint64_t end_row = 0;
    for (uint64_t i = 0; i < r; i++) {
        end_row += get_n(i);
        std::cout << end_row << "\t" << get_char(i) << "\n";
        counter += 1;
        if (i%100000 == 0) std::cerr << i << "\t" << split_thresholds << "\r";
        for (int j = 0; j < 16; j ++) {
            if (get_n(i) >= std::pow(2,j + 1)) {
                counts_length[j] += 1;
            }
            if (get_offset(i) >= std::pow(2,j + 1)) {
                counts_offset[j] += 1;
            }
#if SPLIT_THRESHOLDS_FALSE
            if (get_thresholds(i, 0) >= std::pow(2,j + 1)) {
                counts_threshold0[j] += 1;
            }
            if (get_thresholds(i, 1) >= std::pow(2,j + 1)) {
                counts_threshold1[j] += 1;
            }
            if (get_thresholds(i, 2) >= std::pow(2,j + 1)) {
                counts_threshold2[j] += 1;
            }
#endif
        }
    }
    std::cerr << "split_thresholds: " << split_thresholds << "\n";
    std::cerr << "counter: " << counter << "\n";
    std::cerr << "\ncounts_length:\n";
    for (int j=0; j < 16; j++) {
        std::cerr << j + 1 << "\t" << counts_length[j] << "\n";
    }
    std::cerr << "\ncounts_offset:\n";
    for (int j=0; j < 16; j++) {
        std::cerr << j + 1 << "\t" << counts_offset[j] << "\n";
    }
    std::cerr << "\ncounts_threshold0:\n";
    for (int j=0; j < 16; j++) {
        std::cerr << j << ": " << counts_threshold0[j] << "\n";
    }
    std::cerr << "\ncounts_threshold1:\n";
    for (int j=0; j < 16; j++) {
        std::cerr << j << ": " << counts_threshold1[j] << "\n";
    }
    std::cerr << "\ncounts_threshold2:\n";
    for (int j=0; j < 16; j++) {
        std::cerr << j << ": " << counts_threshold2[j] << "\n";
    }
}

void MoveStructure::print_stats() {

    std::cerr << "n:\t" << length << "\n";
    std::cerr << "r:\t" << r << "\n";
    std::cerr << "n/r:\t" << static_cast<double>(length)/r << "\n";
    std::cerr << "original_r:\t" << original_r << "\n";
    std::cerr << "end_bwt_idx:\t" << end_bwt_idx << "\n";

    for (int i = 0; i < first_runs.size(); i++) {
        std::cerr << alphabet[i] << "\t" << i << "\t" << first_runs[i] << "\t" << last_runs[i] << "\n";
    }

    if (original_r != 0) {
        std::cerr << "original_r: " << original_r << "\n";
        std::cerr << "n/original_r: " << static_cast<double>(length)/original_r << "\n";
    }
    std::cerr << "Size of the rlbwt table: " << sizeof(rlbwt[0]) * rlbwt.size() * (0.000000001) << "\n";
#if BLOCKED_MODE
    std::cerr << "Size of the block table: " << sizeof(id_blocks[0][0]) * id_blocks[0].size() * id_blocks.size() * (0.000000001) << "\n";
#endif
#if TALLY_MODE
    std::cerr << "Size of the tally table: " << sizeof(tally_ids[0][0]) * tally_ids[0].size() * tally_ids.size() * (0.000000001) << "\n";
#endif

    if (movi_options->is_output_ids()) {
        print_ids();
    }
}

void MoveStructure::write_doc_set_freqs(std::string fname) {
    // Get doc set counts.
    doc_set_cnts.resize(unique_doc_sets.size());
    for (size_t i = 0; i < r; i++) {
        doc_set_cnts[doc_set_inds[i]]++;
    }

    std::vector<std::pair<uint64_t, uint32_t>> freqs(unique_doc_sets.size());
    for (size_t i = 0; i < freqs.size(); i++) {
        freqs[i].first = doc_set_cnts[i];
        freqs[i].second = i;
    }
    sort(freqs.begin(), freqs.end(), std::greater<>());

    std::ofstream out(fname);
    for (size_t i = 0; i < 100000; i++) {
        out << freqs[i].first << " ";
        for (int doc : unique_doc_sets[freqs[i].second]) {
            out << doc << " ";
        }
        out << "\n";
    }
    out.close();
}

void MoveStructure::print_ids() {

    std::cerr << "\nThe ids are being written out..\n";

    std::string base_name = movi_options->get_index_dir() + "/ids";
    std::vector<std::unique_ptr<std::ofstream>> id_files;

    std::ofstream ids_all(base_name + ".all");

    for (int i = 0; i < alphabet.size(); i++) {
        std::string ext(1, alphabet[i]);
        id_files.push_back(std::make_unique<std::ofstream>(base_name + "." + ext));
    }

    for (uint64_t i = 0; i < r; i++) {

        if (i%100000 == 0) std::cerr << i << "\r";

        if (i != end_bwt_idx) {

            uint64_t id = get_id(i);
            uint64_t adjusted_id = id - first_runs[rlbwt[i].get_c() + 1];

            int char_index = static_cast<int>(rlbwt[i].get_c());
            *id_files[char_index] << adjusted_id << "\t" << i << "\n";

            ids_all << adjusted_id << "\n";
        }
    }

    ids_all.close();
    for (auto& id_file : id_files) {
        id_file->close();
    }

    std::cerr << "\nDone.\n";
}
