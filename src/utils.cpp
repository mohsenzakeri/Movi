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
        return "pml";
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
    std::cerr << "kseq destroyed!\n";
    gzclose(fp); // STEP 6: close the file handler
    std::cerr << "fp file closed!\n";
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

void read_thresholds(std::string tmp_filename, sdsl::int_vector<>& thresholds) {
    int log_n = 100;

    struct stat filestat;
    FILE *fd;

    if ((fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
        std::cerr <<("open() file " + tmp_filename + " failed");

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        std::cerr <<("stat() file " + tmp_filename + " failed");

    if (filestat.st_size % THRBYTES != 0)
        std::cerr <<("invilid file " + tmp_filename);

    size_t length_thr = filestat.st_size / THRBYTES;
    size_t threshold = 0;

    thresholds = sdsl::int_vector<>(length_thr, 0, log_n);

    size_t i = 0;
    for (i = 0; i < length_thr; ++i) {
        if (i % 100000 == 0) {
            std::cerr << "read thresholds:\t" << i << "\r";
        }
        size_t threshold = 0;
        if ((fread(&threshold, THRBYTES, 1, fd)) != 1)
            std::cerr <<("fread() file " + tmp_filename + " failed");
        thresholds[i] = threshold;
    }
    std::cerr << "Finished reading " << i << " thresholds.\n";
}


void output_matching_lengths(bool stdout, std::ofstream& mls_file, std::string read_id, std::vector<uint16_t>& matching_lengths) {

    uint64_t matching_lengths_size = matching_lengths.size();

    if (stdout) {
        // std::cout << ">" << seq->name.s << " \n";
        printf(">%s \n", read_id);
        for (int64_t i = matching_lengths_size - 1; i >= 0; i--) {
            printf("%d ", matching_lengths[i]);
        }
        printf("\n");
    } else {
        // uint16_t st_length = seq->name.m;
        uint16_t st_length = read_id.length();
        mls_file.write(reinterpret_cast<char*>(&st_length), sizeof(st_length));
        // mls_file.write(reinterpret_cast<char*>(&seq->name.s[0]), st_length);
        mls_file.write(reinterpret_cast<char*>(&read_id[0]), st_length);
        mls_file.write(reinterpret_cast<char*>(&matching_lengths_size), sizeof(matching_lengths_size));
        mls_file.write(reinterpret_cast<char*>(&matching_lengths[0]), matching_lengths_size * sizeof(matching_lengths[0]));
    }
}

void output_counts(bool stdout, std::ofstream& count_file, std::string read_id, size_t query_length, int32_t pos_on_r, uint64_t match_count) {
    if (stdout) {
        // std::cout << seq->name.s << "\t";
        std::cout << read_id << "\t";
        std::cout << query_length - pos_on_r << "/" << query_length << "\t" << match_count << "\n";
    } else {
        // count_file << seq->name.s << "\t";
        count_file << read_id << "\t";
        count_file << query_length - pos_on_r << "/" << query_length << "\t" << match_count << "\n";
    }
}

void output_logs(std::ofstream& costs_file, std::ofstream& scans_file, std::ofstream& fastforwards_file, std::string read_id, MoveQuery& mq) {
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