#include "move_structure.hpp"

std::ofstream MoveStructure::open_index_write() {
    mkdir(movi_options->get_index_dir().c_str(), 0777);
    std::string fname = movi_options->get_index_dir() + "/index.movi";
    if (movi_options->is_color_move_rows()) {
        fname = movi_options->get_index_dir() + "/index_colored.movi";
    }
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    if (!fout) {
        throw std::runtime_error(ERROR_MSG("Failed to open the index file at: " + movi_options->get_index_dir()));
    }
    return fout;
}

std::ifstream MoveStructure::open_index_read() {
    std::string fname;
#if COLOR_MODE == 1
    fname = movi_options->get_index_dir() + "/index_colored.movi";
#else
    fname = movi_options->get_index_dir() + "/index.movi";
#endif

    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    if (!fin) {
        // Attempt to read an index file built with the old index name
        fname = movi_options->get_index_dir() + "/movi_index.bin";
        fin.open(fname, std::ios::in | std::ios::binary);
        if (!fin) {
            throw std::runtime_error(ERROR_MSG("[deserialize index] Failed to open the index file at: " + movi_options->get_index_dir()));
        }
    }

    if (movi_options->is_verbose()) {
        std::cerr << "Index file name: " << fname << std::endl;
    }

    fin.seekg(0, std::ios::beg);
    return fin;
}

void MoveStructure::write_index_header(std::ofstream& fout) {
    if (!movi_options->is_no_header()) {
        if (movi_options->is_legacy_header()) {
            // Legacy header (version 1.x) - just a single byte for mode
            char index_type = static_cast<char>(MODE);
            fout.write(reinterpret_cast<char*>(&index_type), sizeof(index_type));
        } else {
            // Modern header (version 2.x) with magic number and metadata
            MoviHeader header;
            char index_type = static_cast<char>(MODE);
            header.init(index_type, length, r, original_r, end_bwt_idx);
            header.write(fout);
        }
    }

    if (movi_options->is_no_header() or movi_options->is_legacy_header()) {
        // Write basic index characteristics
        fout.write(reinterpret_cast<char*>(&length), sizeof(length));
        fout.write(reinterpret_cast<char*>(&r), sizeof(r));
        fout.write(reinterpret_cast<char*>(&end_bwt_idx), sizeof(end_bwt_idx));
    }

    if (movi_options->is_verbose()) {
        std::cerr << "length: " << length << " r: " << r << " original_r: " << original_r << " end_bwt_idx: " << end_bwt_idx << "\n";
    }
}

void MoveStructure::read_index_header(std::ifstream& fin) {
    if (!movi_options->is_no_header()) {
        if (movi_options->is_legacy_header()) {
            // Legacy header - just a single byte for mode
            char index_type;
            fin.read(reinterpret_cast<char*>(&index_type), sizeof(index_type));

            // Verify mode matches
            if (static_cast<char>(index_type) != static_cast<char>(MODE)) {
                throw std::runtime_error(ERROR_MSG("[deserialize index] Index mode mismatch. Expected mode " +
                    std::to_string(MODE) + " but index has mode " + std::to_string(static_cast<uint8_t>(index_type))));
            }
        } else {
            // Modern header (version 2.x)
            MoviHeader header;
            fin.read(reinterpret_cast<char*>(&header), sizeof(MoviHeader));

            // Verify magic number
            if (header.magic != MOVI_MAGIC) {
                throw std::runtime_error(ERROR_MSG("[deserialize index] Invalid magic number in header. Not a valid Movi 2.0 index file."));
            }

            // Verify mode matches
            if (header.type != static_cast<uint8_t>(MODE)) {
                throw std::runtime_error(ERROR_MSG("[deserialize index] Index mode mismatch. Expected mode " +
                    std::to_string(MODE) + " but index has mode " + std::to_string(header.type)));
            }

            length = header.length;
            r = header.r;
            original_r = header.original_r;
            end_bwt_idx = header.end_bwt_idx;
        }
        std::cerr << GENERAL_MSG("The " + program() + " index is being used.");
    }

    // Read basic index characteristics for older versions
    if (movi_options->is_no_header() or movi_options->is_legacy_header()) {
        fin.read(reinterpret_cast<char*>(&length), sizeof(length));
        fin.read(reinterpret_cast<char*>(&r), sizeof(r));
        fin.read(reinterpret_cast<char*>(&end_bwt_idx), sizeof(end_bwt_idx));
    }

    std::cerr << INFO_MSG("Basic index characteristics:\n\tn: " + std::to_string(length) + "\n\tr: " + std::to_string(r) + "\n\t$: " + std::to_string(end_bwt_idx));
}

void MoveStructure::write_basic_index_data(std::ofstream& fout) {
    // Write the thresholds of the end_bwt_ids
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_thresholds[0]), 4*sizeof(end_bwt_idx_thresholds[0]));

    // Write the pointers of the end_bwt_ids for repositioning in constant mode
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_next_down[0]), 4*sizeof(end_bwt_idx_next_down[0]));
    fout.write(reinterpret_cast<char*>(&end_bwt_idx_next_up[0]), 4*sizeof(end_bwt_idx_next_up[0]));

    // Write alphamap
    uint64_t alphamap_size = alphamap.size();
    fout.write(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    fout.write(reinterpret_cast<char*>(&alphamap[0]), alphamap.size()*sizeof(alphamap[0]));

    // Write list of the symbols in the alphabet
    uint64_t alphabet_size = alphabet.size();
    fout.write(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));    
    fout.write(reinterpret_cast<char*>(&alphabet[0]), alphabet.size()*sizeof(alphabet[0]));

    // The fallowing flag are not used any more, the header information is enough
    fout.write(reinterpret_cast<char*>(&nt_splitting), sizeof(nt_splitting));
    fout.write(reinterpret_cast<char*>(&constant), sizeof(constant));

    // Only stored in old Movi versions
    if (movi_options->is_no_header() or movi_options->is_legacy_header()) {
        // This is not used any more as the onebit modes is deprecated
        bool onebit = false;
        fout.write(reinterpret_cast<char*>(&onebit), sizeof(onebit));
    }
}

void MoveStructure::read_basic_index_data(std::ifstream& fin) {
    // Read the thresholds of the end_bwt_ids
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_thresholds[0]), 4*sizeof(end_bwt_idx_thresholds[0]));

    // Read the pointers of the end_bwt_ids for repositioning in constant mode
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_next_down[0]), 4*sizeof(end_bwt_idx_next_down[0]));
    fin.read(reinterpret_cast<char*>(&end_bwt_idx_next_up[0]), 4*sizeof(end_bwt_idx_next_up[0]));

    // Read alphamap
    uint64_t alphamap_size;
    fin.read(reinterpret_cast<char*>(&alphamap_size), sizeof(alphamap_size));
    alphamap.resize(alphamap_size);
    fin.read(reinterpret_cast<char*>(&alphamap[0]), alphamap_size*sizeof(alphamap[0]));

    // Read the list of the symbols in the alphabet
    uint64_t alphabet_size;
    fin.read(reinterpret_cast<char*>(&alphabet_size), sizeof(alphabet_size));
    if (movi_options->is_verbose()) {
        std::cerr << "alphabet_size: " << alphabet_size << "\n";
    }
    alphabet.resize(alphabet_size);
    fin.read(reinterpret_cast<char*>(&alphabet[0]), alphabet_size*sizeof(alphabet[0]));
    if (alphabet.size() > 4) {
        std::cerr << "Warning: There are more than 4 characters, the index expects only A, C, T and G in the reference.\n";
    }

    fin.read(reinterpret_cast<char*>(&nt_splitting), sizeof(nt_splitting));
    fin.read(reinterpret_cast<char*>(&constant), sizeof(constant));

    // Only stored in old Movi versions
    if (movi_options->is_no_header() or movi_options->is_legacy_header()) {
        // This is not used any more as the onebit modes is deprecated
        bool onebit;
        fin.read(reinterpret_cast<char*>(&onebit), sizeof(onebit));
    }
}

void MoveStructure::write_overflow_tables(std::ofstream& fout) {
    // Write the overflow table (Used in Movi large, constant, and split)
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

    // The following is not related to overflow tables
    // Only stored in old Movi versions
    if (movi_options->is_no_header() or movi_options->is_legacy_header()) {
        // Write the length of the original string used for building the index
        size_t orig_size = 0;
        fout.write(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));

        // A flag to indicate that the original string was reconstructed
        bool reconstructed = false;
        fout.write(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));

        // The eof_row is not used anymore, it is the same as end_bwt_idx
        uint64_t eof_row = 0;
        fout.write(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));
    }
}

void MoveStructure::read_overflow_tables(std::ifstream& fin) {
    // Read the overflow table (Used in Movi large, constant, and split)
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

    // The following is not related to overflow tables
    // Only stored in old Movi versions
    if (movi_options->is_no_header() or movi_options->is_legacy_header()) {
        // Read the length of the original string used for building the index
        size_t orig_size;
        fin.read(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));

        // A flag to indicate that the original string was reconstructed
        bool reconstructed = false;
        fin.read(reinterpret_cast<char*>(&reconstructed), sizeof(reconstructed));

        // The eof_row is not used anymore, it is the same as end_bwt_idx
        uint64_t eof_row = 0;
        fin.read(reinterpret_cast<char*>(&eof_row), sizeof(eof_row));
    }
}

void MoveStructure::write_counts_data(std::ofstream& fout) {
    uint64_t counts_size = counts.size();
    fout.write(reinterpret_cast<char*>(&counts_size), sizeof(counts_size));
    fout.write(reinterpret_cast<char*>(&counts[0]), counts.size()*sizeof(counts[0]));

    uint64_t last_runs_size = last_runs.size();
    fout.write(reinterpret_cast<char*>(&last_runs_size), sizeof(last_runs_size));
    fout.write(reinterpret_cast<char*>(&last_runs[0]), last_runs.size()*sizeof(last_runs[0]));
    fout.write(reinterpret_cast<char*>(&last_offsets[0]), last_offsets.size()*sizeof(last_offsets[0]));
    fout.write(reinterpret_cast<char*>(&first_runs[0]), first_runs.size()*sizeof(first_runs[0]));
    fout.write(reinterpret_cast<char*>(&first_offsets[0]), first_offsets.size()*sizeof(first_offsets[0]));
}

void MoveStructure::read_counts_data(std::ifstream& fin) {
    // Read the counts of each character in the alphabet
    uint64_t counts_size = 0;
    fin.read(reinterpret_cast<char*>(&counts_size), sizeof(counts_size));
    counts.resize(counts_size);
    fin.read(reinterpret_cast<char*>(&counts[0]), counts_size*sizeof(uint64_t));

    // The backward search intervals for every sequence of length one (one for each character in the alphabet)
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
}

#if BLOCKED_MODES
void MoveStructure::write_id_blocks(std::ofstream& fout) {
    // Write the id blocks table for blocked indexes
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
    fout.write(reinterpret_cast<char*>(&block_size), sizeof(block_size));
}

void MoveStructure::read_id_blocks(std::ifstream& fin) {
    // Read the id blocks table for blocked indexes
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
    if (!fin.eof() and movi_options->is_adjusted_block()) {
        fin.read(reinterpret_cast<char*>(&block_size), sizeof(block_size));
    }
}
#endif

#if TALLY_MODES
void MoveStructure::write_tally_table(std::ofstream& fout) {
    // Write the tally table for tally indexes
    fout.write(reinterpret_cast<char*>(&tally_checkpoints), sizeof(tally_checkpoints));
    uint64_t tally_ids_len = tally_ids[0].size();
    fout.write(reinterpret_cast<char*>(&tally_ids_len), sizeof(tally_ids_len)); 
    for (uint32_t i = 0; i < alphabet.size(); i++) {
        fout.write(reinterpret_cast<char*>(&tally_ids[i][0]), tally_ids[i].size()*sizeof(tally_ids[i][0]));
    }
}

void MoveStructure::read_tally_table(std::ifstream& fin) {
    // Read the tally table for tally indexes
    fin.read(reinterpret_cast<char*>(&tally_checkpoints), sizeof(tally_checkpoints));
    std::cerr << "Tally mode with tally_checkpoints = " << tally_checkpoints << "\n";
    uint64_t tally_ids_len = 0;
    fin.read(reinterpret_cast<char*>(&tally_ids_len), sizeof(tally_ids_len));
    tally_ids.resize(alphabet.size());
    for (uint32_t i = 0; i < alphabet.size(); i++) {
        tally_ids[i].resize(tally_ids_len);
        fin.read(reinterpret_cast<char*>(&tally_ids[i][0]), tally_ids_len*sizeof(MoveTally));
    }
}
#endif

void MoveStructure::write_main_table(std::ofstream& fout) {
    // Write the main Movi table
    if (movi_options->is_color_move_rows()) {
        fout.write(reinterpret_cast<char*>(&rlbwt_colored[0]), rlbwt_colored.size()*sizeof(rlbwt_colored[0]));
    } else {
        fout.write(reinterpret_cast<char*>(&rlbwt[0]), rlbwt.size()*sizeof(rlbwt[0]));
    }
}

void MoveStructure::read_main_table(std::ifstream& fin, std::streamoff rlbwt_offset) {
    size_t rlbwt_size = r * sizeof(MoveRow);

    if (movi_options->is_mmap()) {
        // The rlbwt table is stored as a file called rlbwt.movi in the index directory
        std::string rlbwt_fname = movi_options->get_index_dir() + "/rlbwt.movi";
        int fd = open(rlbwt_fname.c_str(), O_RDONLY);
        if (fd == -1) {
            throw std::runtime_error(ERROR_MSG("Failed to open index file: " + std::string(strerror(errno))));
        }

        // Get the file size
        struct stat sb;
        if (fstat(fd, &sb) == -1) {
            close(fd);
            throw std::runtime_error(ERROR_MSG("fstat failed: " + std::string(strerror(errno))));
        }

        void* rlbwt_mem = mmap(nullptr, rlbwt_size, PROT_READ, MAP_SHARED, fd, 0);
        if (rlbwt_mem == MAP_FAILED) {
            close(fd);
            throw std::runtime_error(ERROR_MSG("Failed to mmap rlbwt: " + std::string(strerror(errno))));
        }

        rlbwt_view = std::span<MoveRow>(static_cast<MoveRow*>(rlbwt_mem), r);

        // Continue reading the rest of the file
        fin.seekg(rlbwt_offset + rlbwt_size, std::ios::beg);
    } else {
        rlbwt.resize(r);
        if (movi_options->is_verbose()) {
            std::cerr << "sizeof(MoveRow): " << sizeof(MoveRow) << std::endl;
        }
        fin.read(reinterpret_cast<char*>(&rlbwt[0]), r*sizeof(MoveRow));
    }
    std::cerr << INFO_MSG("All the move rows are read.");
}

void MoveStructure::serialize() {
    std::ofstream fout = open_index_write();

    write_index_header(fout);

    write_basic_index_data(fout);

    write_main_table(fout);

#if TALLY_MODES
    write_tally_table(fout);
#endif

    write_overflow_tables(fout);

    write_counts_data(fout);

#if BLOCKED_MODES
    write_id_blocks(fout);
#endif

    // Only stored in old Movi versions
    if (movi_options->is_no_header() or movi_options->is_legacy_header()) {
        // Write the original_r (number of rows before splitting)
        fout.write(reinterpret_cast<char*>(&original_r), sizeof(original_r));
    }

    fout.close();
}

void MoveStructure::deserialize() {
    std::ifstream fin = open_index_read();

    read_index_header(fin);

    read_basic_index_data(fin);

    // Read the main Movi table
    std::streamoff rlbwt_offset = fin.tellg();
    read_main_table(fin, rlbwt_offset);

#if TALLY_MODES
    read_tally_table(fin);
#endif

    read_overflow_tables(fin);

    read_counts_data(fin);

#if BLOCKED_MODES
    read_id_blocks(fin);
#endif

    if (movi_options->is_no_header() or movi_options->is_legacy_header()) {
        // To be able to load the indexes that haven't stored original_r
        if (fin.eof()) {
            std::cerr << "original_r was not stored in the index!\n";
        } else {
            // Read the original_r (number of rows before splitting)
            fin.read(reinterpret_cast<char*>(&original_r), sizeof(original_r));
        }
    }

    fin.close();
}

void MoveStructure::flat_and_serialize_colors_vectors() {
    std::unordered_map<uint32_t, uint64_t> flat_index;
    for (uint64_t i = 0; i < unique_doc_sets.size(); i++) {
        if (i % 1000000 == 0) {
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
        if (i % 1000000 == 0) {
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

void MoveStructure::serialize_doc_pats(std::string fname) {
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    std::cerr << "Writing doc pats to: " << fname << std::endl;

    fout.write(reinterpret_cast<char*>(&doc_pats[0]), length * sizeof(doc_pats[0]));
    fout.close();
}

void MoveStructure::deserialize_doc_pats(std::string fname) {
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    if (!fin.good()) {
        throw std::runtime_error(ERROR_MSG("[deserialize doc pats] Failed to open document patterns file at " + fname));
    }

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

void MoveStructure::deserialize_doc_sets_flat() {
    std::string fname = movi_options->get_index_dir() + "/doc_sets_flat.bin";
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    if (!fin.good()) {
        throw std::runtime_error(ERROR_MSG("[deserialize doc sets flat] Failed to open document sets flat file at " + fname));
    }

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
    if (!fin.good()) {
        throw std::runtime_error(ERROR_MSG("[deserialize doc sets] Failed to open document sets file at " + fname));
    }

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

void MoveStructure::load_document_info() {
    // Read in document offsets.
    std::string doc_offsets_fname = movi_options->get_index_dir() + "/ref.fa.doc_offsets";
    std::ifstream doc_offsets_file(doc_offsets_fname);
    if (!doc_offsets_file.good()) {
        throw std::runtime_error(ERROR_MSG("[load document info] doc_offsets file not found at \"" + doc_offsets_fname + "\""));
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
    std::string doc_ids_fname = movi_options->get_index_dir() + "/ref.fa.doc_ids";
    std::ifstream doc_ids_file(doc_ids_fname);
    if (doc_ids_file.good()) {
        uint32_t doc_id;
        while ((doc_ids_file >> doc_id)) {
            doc_ids.push_back(doc_id);
            taxon_id_compress[doc_id] = 0;
        }
        doc_ids_file.close();
    } else {
        // No doc ids information
        std::cerr << GENERAL_MSG("No doc ids information, setting doc ids to 1...");
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

void MoveStructure::serialize_sampled_SA() {
    std::string fname = movi_options->get_index_dir() + "/ssa.movi";
    std::ofstream fout(fname, std::ios::out | std::ios::binary);
    uint64_t SA_sample_rate = movi_options->get_SA_sample_rate();
    fout.write(reinterpret_cast<char*>(&SA_sample_rate), sizeof(SA_sample_rate));
    uint64_t sampled_SA_entries_size = sampled_SA_entries.size();
    fout.write(reinterpret_cast<char*>(&sampled_SA_entries_size), sizeof(sampled_SA_entries_size));
    fout.write(reinterpret_cast<char*>(&sampled_SA_entries[0]), sampled_SA_entries_size*sizeof(sampled_SA_entries[0]));
    uint64_t all_p_size = all_p.size();
    fout.write(reinterpret_cast<char*>(&all_p_size), sizeof(all_p_size));
    fout.write(reinterpret_cast<char*>(&all_p[0]), all_p_size*sizeof(all_p[0]));
    fout.close();
}

void MoveStructure::deserialize_sampled_SA() {
    std::string fname = movi_options->get_index_dir() + "/ssa.movi";
    std::ifstream fin(fname, std::ios::in | std::ios::binary);
    if (!fin.good()) {
        throw std::runtime_error(ERROR_MSG("[deserialize sampled SA] Failed to open sampled SA entries file at " + fname +
                                           "\nBuild the sampled SA by running the build-SA command."));
    }

    uint64_t SA_sample_rate = 0;0;
    fin.read(reinterpret_cast<char*>(&SA_sample_rate), sizeof(SA_sample_rate));
    movi_options->set_SA_sample_rate(SA_sample_rate);
    uint64_t sampled_SA_entries_size = 0;
    fin.read(reinterpret_cast<char*>(&sampled_SA_entries_size), sizeof(sampled_SA_entries_size));
    sampled_SA_entries.resize(sampled_SA_entries_size);
    fin.read(reinterpret_cast<char*>(&sampled_SA_entries[0]), sampled_SA_entries_size*sizeof(sampled_SA_entries[0]));
    uint64_t all_p_size = 0;
    fin.read(reinterpret_cast<char*>(&all_p_size), sizeof(all_p_size));
    all_p.resize(all_p_size);
    fin.read(reinterpret_cast<char*>(&all_p[0]), all_p_size*sizeof(all_p[0]));
    fin.close();
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
    for (size_t i = 0; i < 1000000; i++) {
        out << freqs[i].first << " ";
        for (int doc : unique_doc_sets[freqs[i].second]) {
            out << doc << " ";
        }
        out << "\n";
    }
    out.close();
}

void MoveStructure::write_ftab() {
    size_t ftab_k = movi_options->get_ftab_k();

    uint64_t ftab_size = std::pow(4, ftab_k);
    if (ftab_size != ftab.size()) {
        throw std::runtime_error(ERROR_MSG("[write ftab] The size of the ftab is not correct: " +
                                 std::to_string(ftab_size) + " != " + std::to_string(ftab.size())));
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
            if (!fin.good()) {
                throw std::runtime_error(ERROR_MSG("[read ftab] Failed to open the ftab file at " + fname));
            }

            fin.read(reinterpret_cast<char*>(&ftab_k), sizeof(ftab_k));

            uint64_t ftab_size = 0;
            fin.read(reinterpret_cast<char*>(&ftab_size), sizeof(ftab_size));
            if (ftab_size != std::pow(4, ftab_k)) {
                throw std::runtime_error(ERROR_MSG("[read ftab] The size of the ftab is not correct: " +
                                 std::to_string(ftab_size) + " != " + std::to_string(std::pow(4, ftab_k))));
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
        if (!fin.good()) {
            throw std::runtime_error(ERROR_MSG("[read ftab] Failed to open the ftab file at " + fname));
        }
        fin.read(reinterpret_cast<char*>(&ftab_k), sizeof(ftab_k));
        uint64_t ftab_size = 0;
        fin.read(reinterpret_cast<char*>(&ftab_size), sizeof(ftab_size));
        std::cerr << GENERAL_MSG("Reading the ftab.\nftab_k: " + std::to_string(ftab_k) + "\t" + "ftab_size: " + std::to_string(ftab_size) + "\n");
        if (ftab_size != std::pow(4, ftab_k)) {
            throw std::runtime_error(ERROR_MSG("[read ftab] The size of the ftab is not correct: " +
                                 std::to_string(ftab_size) + " != " + std::to_string(std::pow(4, ftab_k))));
        }

        ftab.resize(ftab_size);
        fin.read(reinterpret_cast<char*>(&ftab[0]), ftab_size*sizeof(ftab[0]));
        fin.close();
    }
}