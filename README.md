# Movi 2 ![GitHub](https://img.shields.io/github/license/mohsenzakeri/movi?color=green)

Movi is a full text index for pangenomes. It takes advantage of the move data structure (by Nishimoto and Tabei). The high locality of the reference and latency hiding in Movi results in low and predictabile query latencies. These properties make Movi ideal for applications like Nanopore Adaptive sampling which requires real-time classification of the reads.

## Getting Started

Required dependences: `sdsl`, `zlib`, `cmake`, and `gcc` with C++20 support.

To build from source:
```
git clone https://github.com/mohsenzakeri/Movi
cd Movi
mkdir build && cd build
cmake ..
make -j4
```

## Usage
```
 movi <action> [OPTION...]

 actions:       build, inspect, query, view

 main options:
      --type arg  The type of the index: regular-thresholds (default), regular, blocked-thresholds, blocked,
                  sampled-thresholds, sampled
  -h, --help      Print help
      --help-all  Print help with all options (including advanced)
      --verbose   Enable verbose mode

 build options:
  -i, --index arg       Index directory [REQUIRED]
  -f, --fasta arg       Reference file [REQUIRED unless -l is passed]
  -l, --list arg        List of fasta files, only works with 'movi' binary [REQUIRED unless -f is passed]
      --separators      Use separators in between fasta entries
      --checkpoint arg  Create checkpoint for id field of sampled and sampled-thresholds indexes every n move rows
      --verify          Verify if all the LF-move operations are correct

 inspect options:
  -i, --index arg           Index directory [REQUIRED]
      --output-ids          Output the adjusted ids of all the runs to ids.* files, one file per character
      --flat-color-vectors  Flat and serialize the colors vectors

 query options:
  -i, --index arg     Index directory [REQUIRED]
  -r, --read arg      fasta/fastq Read file for query [REQUIRED]
  -o, --out-file arg  Output file prefix if computing PMLs for classification
  -t, --threads arg   Number of threads for query
      --pml           Compute the pseudo-matching lengths (default)
      --zml           Compute the Ziv-Merhav cross parsing length)
      --count         Compute the count queries
      --classify      Enable binary classification of the reads
      --filter        Filter the reads based on the matching lengths, output the filtered reads to stdout
  -v, --invert        Output the not found reads during filtering
      --stdout        Write the output to stdout, writes the matching lengths by default, or the report of classification
                      if --classify is passed

 view options:
      --bpf arg    The base profile format (BPF) file to view [REQUIRED]
      --small-bpf  Read the file with PMLs stored as uint16_t (default: uint32_t).
      --large-bpf  Read the file with PMLs stored as uint64_t (default: uint32_t)
```

## Build the Movi index

```
./movi build --fasta <fasta file> --index <index directory> --type <index type>
```

`<fasta file>` is a file including the reference genomes to be indexed.

`<index directory>` is the directory where you want the Movi index to be located.

Note: The flag `--type <index type>` determines the strategy to build the main table of the Movi index. If no value is passed, the `regular-thresholds` index is built.
Possible index types: `large` `constant` `split` `regular` `regular-thresholds` `blocked` `blocked-thresholds` `sampled` `sampled-thresholds`

If the thresholds are note available, please pass `--type regular` to build an index without the thresholds data structure.

## Compute Pseudo Matching Lengths (PML)

To compute PMLs using the movi index, please run the following command on the fastq or fasta file of the reads:
```
./movi query --index <index directory> --read <reads file> -o <output prefix>
```

`<reads file>` is the address of the fasta or fastq file containing the reads.

After the query command finishes, the output will be stored at `<output prefix>.<query type>.bpf`. Since this file is in the binary format, to view the PMLs please run the following command:
```
./movi view --bpf <bpf file> | less
```
`<bpf file>` is the file generated in the query step.

The output of the last command shows each read's name following by pseudo matching lengths computed for it. A pseudo matching length is outputed for every base of the read.

Alternatively, you can directly output the PMLs to a text file using the `--stdout` flag:
```
./movi query --index <index directory> --read <reads file> --stdout > <output file>
```
## Binary classification of the reads

You can apply different classification schemes on the pseudo matching lengths produced by Movi to determine whether the read is found in the index or not.
Movi implements the classification scheme adopted from the [SPUMONI](https://github.com/oma219/spumoni/tree/main) software written by Omar Y. Ahmed.
For binary classification you should pass the `--classify` flag to the query step:
```
./movi query --index <index directory> --read <reads file> --classify
```
Then, a file called `<reads file>.<index type>.<query type>.report` will be generated.

## Citing Movi
If you use  Movi in yoour research project, please cite:

>[Zakeri, Mohsen, Brown, Nathaniel K., Ahmed, Omar Y., Gagie, Travis, and Langmead, Ben. "Movi: a fast and cache-efficient full-text pangenome index". iScience (2024)](https://www.cell.com/iscience/fulltext/S2589-0042(24)02691-9)

>[Nishimoto, Takaaki, and Yasuo Tabei. "Optimal-time queries on BWT-runs compressed indexes." arXiv preprint arXiv:2006.05104 (2020)](https://arxiv.org/abs/2006.05104).


# Multi-class classification with Movi Color

Movi Color augments the Movi index with information about the origin of sequences, enabling multi-class and taxonomic classification.

>[Tan, Steven, Sina Majidian, Ben Langmead, and Mohsen Zakeri. "Movi Color: fast and accurate long-read classification with the move structure." bioRxiv (2025)](https://www.biorxiv.org/content/10.1101/2025.05.22.655637v1.abstract)
### Add colors to an existing Movi index

**Note: The following step is not required if the index was built using the `--color` flag.**

```
./movi color --index <index directory>
```
This step requires the file `ref.fa.doc_offsets` to be available in the index directory. Every line in the `ref.fa.doc_offsets` should specify the offset in the concatenated text corresponding to the rigth most end of a document + 1.
It will avoid rebuilding the index from the scratch and only adds the color information to an index.

Alternatively, you can pass the --color flag to the build command to construct the color table during index building.

### Perform Multi-class classification
```
./movi query --index <index directory> --read <reads file> --multi-classify --out-file <output file>
```

The `<output file>` will be in the csv format, each line will include three columns: `read name`, `primary document id`, `secondary document id`. 

For more information about the document IDs and multi-class classification, please see [this page](https://github.com/mohsenzakeri/Movi/wiki/Multi%E2%80%90class-classification-with-Movi-Color).
