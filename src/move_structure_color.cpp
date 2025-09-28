#include "move_structure.hpp"

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

// Finds document sets for each run in the rlbwt.
void MoveStructure::build_doc_sets() {
    // Stores set of documents appearing in rows in each run.
    std::unordered_map<DocSet, uint32_t> unique; 
    doc_set_inds.resize(r);

    uint64_t unique_cnt = 0;
    uint16_t temp[MAX_RUN_LENGTH];
    for (uint64_t i = 0; i < r; i++) {
        if (i % 1000000 == 0) {
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
    if (!fin.good()) {
        throw std::runtime_error(ERROR_MSG("[build_tree_doc_sets] Failed to open the doc set similarities file: " + movi_options->get_index_dir() + "/doc_set_similarities.txt"));
    }

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
        // debug_out used to be a member of MoveStructure class
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

void MoveStructure::compute_run_lcs() {
    std::cout << "run,length,lcs\n";
    for (int i = 0; i < r; i++) {
        if (i % 1000000 == 0) std::cerr << i << "\r";
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

void MoveStructure::add_colors_to_rlbwt() {
#if MODE == 3 or MODE == 6
    rlbwt_colored.resize(rlbwt.size());
    for (uint64_t i = 0; i < rlbwt.size(); i++) {
        if (i % 1000000 == 0) {
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

void MoveStructure::compute_color_ids_from_flat() {

    std::cerr << "Start computing the color id table from offsets...\n";

    for (uint64_t color_offset = 0; color_offset < flat_colors.size(); ) {

        if (movi_options->is_report_color_ids()) {
            uint32_t color_id = color_offset_to_id.size();
            color_offset_to_id[color_offset] = color_id;
        }

        uint64_t next_color =  color_offset + flat_colors[color_offset] + 1;
        color_offset = next_color;

        if (num_colors % 1000000 == 0)
            std::cerr << color_offset << "/" << flat_colors.size()
                      << " (" << std::round(static_cast<double>(color_offset)/static_cast<double>(flat_colors.size()) * 100.0) << "%)\r";
        num_colors += 1;

    }

    std::cerr << "\n";
    std::cerr << "Number of colors: " << num_colors << "\n";

    if (movi_options->is_report_color_ids()) {
        std::cerr << "Color offset to color id table is created.\n";
    }

}

void MoveStructure::initialize_classify_cnts() {
    classify_cnts.resize(num_species, 0);
    doc_scores.resize(num_species, 0);
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