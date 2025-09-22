#include "utils.hpp"

uint32_t alphamap_3[4][4] = {{3, 0, 1, 2},
                             {0, 3, 1, 2},
                             {0, 1, 3, 2},
                             {0, 1, 2, 3}};

std::string program() {
#if MODE == 0
    return "large"; // old name: "default"
#endif
#if MODE == 1
    return "constant";
#endif
#if MODE == 4 // Like the default constant mode, but without the pointers to the neighbors with the other characters
    return "split";
#endif
#if MODE == 3
    return "regular"; // old name: compact
#endif
#if MODE == 6
    return "regular-thresholds"; // old name: compact-thresholds
#endif
#if MODE == 2
    return "blocked";
#endif
#if MODE == 8
    return "blocked-thresholds";
#endif
#if MODE == 5
    return "tally";
#endif
#if MODE == 7
    return "tally-thresholds";
#endif
    std::cerr << "The mode is not defined.";
    exit(0);
}

std::string query_type(MoviOptions& movi_options) {
    if (movi_options.is_pml()) {
        if (movi_options.is_random_repositioning()) {
            return "rpml";
        } else {
            return "pml";
        }
    } else if (movi_options.is_zml()) {
        return "zml";
    } else {
        // TODO: handle other query types
        throw std::runtime_error("Invalid query type");
    }
}

kseq_t* open_kseq(gzFile& fp, std::string file_address) {
    std::cerr << "file_address: " << file_address << "\n";
    kseq_t *seq;
    // Solution for handling the stdin input: https://biowize.wordpress.com/2013/03/05/using-kseq-h-with-stdin/
    if (file_address == "-") {
        FILE *instream = stdin;
        fp = gzdopen(fileno(instream), "r"); // STEP 2: open the file handler
    } else {
        fp = gzopen(file_address.c_str(), "r"); // STEP 2: open the file handler
    }
    seq = kseq_init(fp); // STEP 3: initialize seq
    return seq;
}

void close_kseq(kseq_t *seq, gzFile& fp) {
    kseq_destroy(seq); // STEP 5: destroy seq
    // std::cerr << "kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    // std::cerr << "fp file closed!\n";
}

char complement(char c) {
    // # is the separator, complement(#) = #
    char c_comp = c == '#' ? '#' : (c == 'A' ? 'T' : ( c == 'C' ? 'G' : (c == 'G' ? 'C' : 'A')));
    return c_comp;
}

std::string reverse_complement(std::string& fw) {
    std::string rc = "";
    rc.resize(fw.size());
    for (int i = fw.size() - 1;  i >= 0; i--) {
        rc[fw.size() - i - 1] = complement(fw[i]);
    }
    return rc;
}

std::string reverse_complement_from_pos(MoveQuery& mq_fw, int32_t pos_on_r, uint64_t match_len) {
    std::string rc = "";
    rc.resize(match_len);
    for (int i = pos_on_r + match_len - 1;  i >= pos_on_r; i--) {
        rc[pos_on_r + match_len - 1 - i] = complement(mq_fw.query()[i]);
    }
    return rc;
}

std::string number_to_kmer(size_t j, size_t m, std::vector<unsigned char>& alphabet, std::vector<uint64_t>& alphamap) {
    std::string kmer = "";
    for (int i = m - 2; i >= 0; i -= 2) {
        size_t pair = (j >> i) & 0b11;
        // std::cerr << i << " " << pair << " ";
        kmer += alphabet[pair + alphamap['A']];
    }
    return kmer;
}

uint64_t kmer_to_number(size_t k, std::string& r, int32_t pos, std::vector<uint64_t>& alphamap, bool rc) {
    if (r.length() < k) {
        std::cerr << "The k does not match the kmer length!\n";
        exit(0);
    }

    uint64_t res = 0;
    for (size_t i = 0; i < k; i++) {
        if (alphamap[static_cast<uint64_t>(r[pos + i])] == alphamap.size()) {
            return std::numeric_limits<uint64_t>::max();
        }
        if (rc) {
            uint64_t char_code = alphamap[complement(r[pos + i])] - alphamap['A'];
            res = (char_code << ((i)*2)) | res;
        } else {
            uint64_t char_code = alphamap[r[pos + i]] - alphamap['A'];
            res = (char_code << ((k-i-1)*2)) | res;
        }
    }
    return res;
}

uint8_t F_char(std::vector<uint64_t>& first_runs, uint64_t run) {
    for (uint8_t i = first_runs.size() - 1; i > 0; i--) {
        if (run >= first_runs[i]) {
            return i - 1;
        }
    }
    std::cerr << "Undefined behaviour.\n";
    exit(0);
}

void read_thresholds(std::string tmp_filename, std::vector<uint64_t>& thresholds) {
    // Read the 5 Bytes thresholds to uin64_t variables

    struct stat filestat;
    FILE *fd;

    if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
        std::cerr <<("open() file " + tmp_filename + " failed");

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        std::cerr <<("stat() file " + tmp_filename + " failed");

    if (filestat.st_size % THRBYTES != 0)
        std::cerr <<("invilid file " + tmp_filename);

    uint64_t length_thr = filestat.st_size / THRBYTES;
    uint64_t threshold = 0;

    thresholds.resize(length_thr);

    // when the thresholds overflow 5 Bytes, a method is implemented to recover the true values
    uint64_t MAX_5_BYTES = 1ULL << 40;
    uint64_t step_count = 0;

    for (uint64_t i = 0; i < length_thr; ++i) {
        if (i % 1000000 == 0) {
            std::cerr << "read thresholds:\t" << i << "\r";
        }

        threshold = 0;
        if ((fread(&threshold, THRBYTES, 1, fd)) != 1)
            std::cerr <<("fread() file " + tmp_filename + " failed");

        // Detect a sudden drop in thresholds value (overflow)
        if (i > 0 and threshold != 0                                           /* Ignore first iteration and the row with sentinel character */
            and threshold < (thresholds[i-1] - step_count * MAX_5_BYTES)/10    /* Unusual drop in the threshold values */
            and (thresholds[i-1] - step_count * MAX_5_BYTES) > MAX_5_BYTES/10) /* last threshold was not too small */ {

            std::cerr << "\n" << threshold << "\t" << thresholds[i-1] << "\t" << (thresholds[i-1] - step_count * MAX_5_BYTES) << "\n";
            step_count += 1;
        }

        thresholds[i] = threshold + step_count * MAX_5_BYTES;
    }

    std::cerr << "Finished reading " << length_thr << " thresholds.\n";
}

template <typename T>
void output_lens(const std::vector<T>& matching_lengths, std::ofstream& mls_file) {

    uint64_t matching_lengths_size = matching_lengths.size();
    mls_file.write(reinterpret_cast<const char*>(&matching_lengths_size), sizeof(matching_lengths_size));

    mls_file.write(reinterpret_cast<const char*>(&matching_lengths[0]), matching_lengths_size * sizeof(T));

}

void output_matching_lengths(bool to_stdout, std::ofstream& mls_file, std::string read_id, MoveQuery& mq, bool color, bool no_output) {

    if (to_stdout) {

        std::cout << ">" << read_id << "\n";
        std::reverse(mq.get_matching_lengths_string().begin(), mq.get_matching_lengths_string().end());
        std::cout << mq.get_matching_lengths_string();
        std::cout << "\n";

    } else {

        uint16_t st_length = read_id.length();
        // uint16_t st_length = seq->name.m;

        mls_file.write(reinterpret_cast<char*>(&st_length), sizeof(st_length));
        // mls_file.write(reinterpret_cast<char*>(&seq->name.s[0]), st_length);

        mls_file.write(reinterpret_cast<char*>(&read_id[0]), st_length);

        if (color) {

            std::vector<uint64_t>& matching_lengths = mq.get_matching_colors();
            output_lens(matching_lengths, mls_file);

        } else {

            std::vector<uint32_t>& matching_lengths = mq.get_matching_lengths();
            output_lens(matching_lengths, mls_file);

        }
    }
}

void output_sa_entries(std::ofstream& sa_entries_file, std::string read_id, MoveQuery& mq, bool no_output) {
    if (!no_output) {
        sa_entries_file << ">" << read_id << "\n";
        for (auto& sa_entry : mq.get_sa_entries()) {
            sa_entries_file << sa_entry << " ";
        }
        sa_entries_file << "\n";
    }
}

void output_counts(bool to_stdout, std::ofstream& count_file, std::string read_id, size_t query_length, int32_t pos_on_r, uint64_t match_count, bool no_output) {
    if (!no_output) {
        if (to_stdout) {
            // std::cout << seq->name.s << "\t";
            std::cout << read_id << "\t";
            std::cout << query_length - pos_on_r << "/" << query_length << "\t" << match_count << "\n";
        } else {
            // count_file << seq->name.s << "\t";
            count_file << read_id << "\t";
            count_file << query_length - pos_on_r << "/" << query_length << "\t" << match_count << "\n";
        }
    }
}

void output_kmers(bool to_stdout, std::ofstream& kmer_file, std::string read_id, size_t all_kmer_count, MoveQuery& mq, bool no_output) {
    if (!no_output) {
        if (to_stdout) {
            std::cout << read_id << "\t";
            std::cout << mq.found_kmer_count << "/" << all_kmer_count << "\t" << mq.get_matching_lengths_string() << "\n";
        } else {
            kmer_file << read_id << "\t";
            kmer_file << mq.found_kmer_count << "/" << all_kmer_count << "\t" << mq.get_matching_lengths_string() << "\n";
        }
    }
}

void output_logs(std::ofstream& costs_file, std::ofstream& scans_file, std::ofstream& fastforwards_file, std::string read_id, MoveQuery& mq, bool no_output) {
    if (!no_output) {
        // costs_file << ">" << seq->name.s << "\n";
        // scans_file << ">" << seq->name.s << "\n";
        costs_file << ">" << read_id << "\n";
        scans_file << ">" << read_id << "\n";
        fastforwards_file << ">" << read_id << "\n";
        for (auto& cost : mq.get_costs()) {
            costs_file << cost.count() << " ";
        }
        for (auto& scan: mq.get_scans()) {
            scans_file << scan << " ";
        }
        for (auto& fast_forward : mq.get_fastforwards()) {
            fastforwards_file << fast_forward << " ";
        }
        costs_file << "\n";
        scans_file << "\n";
        fastforwards_file << "\n";
    }
}

void output_read(std::string& read_id, std::string& read, bool no_output) {
    if (!no_output) {
        std::cout << ">" << read_id << "\n";
        std::cout << read << "\n";
    }
}

void open_output_files(MoviOptions& movi_options, OutputFiles& output_files) {

    std::string index_type = program();

    // Handle multi_classify out_file (used in ReadProcessor)
    if (movi_options.is_multi_classify()) {
        output_files.out_file = std::ofstream(movi_options.get_out_file());
    }

    bool should_open_files = ((!movi_options.is_stdout() || movi_options.is_classify()) &&
                             !movi_options.is_filter() &&
                             !movi_options.is_no_output());

    if (should_open_files) {

        // Handle color files
        if (movi_options.is_report_colors()) {
            std::string colors_file_name = movi_options.get_read_file()  + "." + index_type + ".colors.bin";
            output_files.colors_file = std::ofstream(colors_file_name, std::ios::out | std::ios::binary);
        } else if (movi_options.is_report_color_ids()) {
            std::string colors_file_name = movi_options.get_read_file()  + "." + index_type + ".color_ids.bin";
            output_files.colors_file = std::ofstream(colors_file_name, std::ios::out | std::ios::binary);
        }

        // Determine output file name prefix based on context
        std::string out_file_name_prefix = movi_options.get_out_file() != "" ? movi_options.get_out_file() :
                                                                               movi_options.get_read_file() + "." + index_type;
        if (!movi_options.is_kmer()) {
            out_file_name_prefix += "." + query_type(movi_options);
        }

        // Handle PML/ZML, count query files and kmer query files
        if (movi_options.is_pml() || movi_options.is_zml()) {
            std::string mls_file_name = out_file_name_prefix + ".bin";
            output_files.mls_file = std::ofstream(mls_file_name, std::ios::out | std::ios::binary);
        } else if (movi_options.is_count()) {
            std::string matches_file_name = out_file_name_prefix + ".matches";
            output_files.matches_file = std::ofstream(matches_file_name);
        } else if (movi_options.is_kmer()) {
            std::string kmer_file_name = out_file_name_prefix + ".kmers." + std::to_string(movi_options.get_k());
            output_files.kmer_file = std::ofstream(kmer_file_name);
        }

        // Handle SA entries file
        if (movi_options.is_get_sa_entries()) {
            std::string sa_entries_file_name = out_file_name_prefix + ".sa_entries";
            output_files.sa_entries_file = std::ofstream(sa_entries_file_name);
        }

        // Handle log files
        if (movi_options.is_logs()) {
            output_files.costs_file = std::ofstream(out_file_name_prefix + ".costs");
            output_files.scans_file = std::ofstream(out_file_name_prefix + ".scans");
            output_files.fastforwards_file = std::ofstream(out_file_name_prefix + ".fastforwards");
        }
    }
}

void close_output_files(MoviOptions& movi_options, OutputFiles& output_files) {
    if (!movi_options.is_no_output() && !movi_options.is_filter()) {
        if (movi_options.is_pml() || movi_options.is_zml()) {
            if (!movi_options.is_stdout()) {
                output_files.mls_file.close();
            }
            if (movi_options.is_get_sa_entries()) {
                output_files.sa_entries_file.close();
            }
            std::cerr << "The output file for the matching lengths closed.\n";
        } else if (movi_options.is_count()) {
            if (!movi_options.is_stdout()) {
                output_files.matches_file.close();
            }
            std::cerr << "The count file is closed.\n";
        }

        if (movi_options.is_logs()) {
            output_files.costs_file.close();
            output_files.scans_file.close();
            output_files.fastforwards_file.close();
        }
    }
}

// Borrowed from spumoni written by Omar Ahmed: https://github.com/oma219/spumoni/tree/main
std::string parse_null_reads(const char* ref_file, const char* output_path) {
    /* Parses out null reads in the case that we don't use a file-list */

    // Variables for generating null reads ...
    std::srand(time(0));
    char grabbed_seq[NULL_READ_CHUNK+1];
    grabbed_seq[NULL_READ_CHUNK] = '\0';

    size_t curr_total_null_reads = 0;
    std::ofstream output_null_fd (output_path, std::ofstream::out);

    // Variables for parsing FASTA ...
    gzFile fp = gzopen(ref_file, "r");
    kseq_t* seq = kseq_init(fp);

    // Go through FASTA file, and extract reads until done
    bool go_for_extraction = (curr_total_null_reads < NULL_READ_BOUND);
    while (kseq_read(seq)>=0 && go_for_extraction) {
        size_t reads_to_grab = (curr_total_null_reads >= NUM_NULL_READS) ? 25 : 100; // downsample if done

        for (size_t i = 0; i < reads_to_grab && go_for_extraction && (seq->seq.l > NULL_READ_CHUNK); i++) {
            size_t random_index = rand() % (seq->seq.l-NULL_READ_CHUNK);
            std::strncpy(grabbed_seq, (seq->seq.s+random_index), NULL_READ_CHUNK);

            // Make sure we don't extract reads of Ns
            if (std::string(grabbed_seq).find("N") == std::string::npos) {
                output_null_fd << ">read_" << curr_total_null_reads << "\n";
                output_null_fd << grabbed_seq << "\n";
                curr_total_null_reads++;
                go_for_extraction = (curr_total_null_reads < NULL_READ_BOUND);
            }
        }

        // Special case - if sequence is less than or equal to 150 bp
        if (seq->seq.l <= NULL_READ_CHUNK) {
            output_null_fd << ">read_" << curr_total_null_reads << "\n";
            output_null_fd << seq->seq.s << "\n";
            curr_total_null_reads++;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    output_null_fd.close();
    return output_path;
}