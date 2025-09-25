# Movi ![GitHub](https://img.shields.io/github/license/mohsenzakeri/movi?color=green)

Movi is a full text index for pangenomes. It takes advantage of the move data structure (by Nishimoto and Tabei). The high locality of the reference and latency hiding in Movi results in low and predictabile query latencies. These properties make Movi ideal for applications like Nanopore Adaptive sampling which requires real-time classification of the reads.

>[Zakeri, Mohsen, Brown, Nathaniel K., Ahmed, Omar Y., Gagie, Travis, and Langmead, Ben. "Movi: a fast and cache-efficient full-text pangenome index". iScience (2024)](https://www.cell.com/iscience/fulltext/S2589-0042(24)02691-9)

>[Nishimoto, Takaaki, and Yasuo Tabei. "Optimal-time queries on BWT-runs compressed indexes." arXiv preprint arXiv:2006.05104 (2020)](https://arxiv.org/abs/2006.05104).

## Install Movi and its dependencies from source

Required dependences: `sdsl`, `zlib`, `cmake`, and `gcc` with C++17 support.

To build from source:
```
git clone https://github.com/mohsenzakeri/Movi
cd Movi
git checkout movi-color
mkdir build-movi-color
cd build-movi-color
cmake ..
make
```

Building the Movi index on a fasta file requires preprocessing the fasta using the [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds) software. 
So, Movi clones and installs this repository during build. Furthermore, for the `split` and `constant` indexes, Movi uses the [r-permute](https://github.com/drnatebrown/r-permute) implemented by Nate Brown, to create the splitting bit vector using the strategy introduced by Nishimoto and Tabei. These preprocessing steps are now included in the Movi build step.


## Build the Movi index

```
./movi build --fasta <fasta file> --index <index directory> --type <index type> --color
```

`<fasta file>` is a file including the reference genomes to be indexed.

`<index directory>` is the directory where you want the Movi index to be located.

The flag `--color` specifies building the Movi Color idnex.

Note: The flag `--type <index type>` determines the strategy to build the main table of the Movi index. If no value is passed, the `regular-thresholds` index is built.
Possible index types: `large` `constant` `split` `regular` `regular-thresholds` `blocked` `blocked-thresholds` `tally` `tally-thresholds`

If the thresholds are note available, please pass `--type regular` to build an index without the thresholds data structure.

## Compute Pseudo Matching Lengths (PML)

To compute PMLs using the movi index, please run the following command on the fastq or fasta file of the reads:
```
./movi query --index <index directory> --read <reads file>
```

`<reads file>` is the address of the fasta or fastq file containing the reads.

After the query command finishes, a file in base profile format (bpf) with the same name as the reads file and the extension `<index type>.<query_type>.bpf` is generated in the directory that also includes the reads file.
Since this file is in the binary format, to view the PMLs please run the following command:
```
./movi view --bpf <bpf file> | less
```
`<bpf file>` is the file generated in the query step.

The output of the last command shows each read's name following by pseudo matching lengths computed for it. A pseudo matching length is outputed for every base of the read. This is the same as the output produced by [SPUMONI](https://github.com/oma219/spumoni).

Alternatively, you can directly output the PMLs to a text file using the `--stdout` flag:
```
./movi query --index <index directory> --read <reads file> --stdout > <output file>
```
## Binary classification of the reads

You can apply different classification schemes on the pseudo matching lengths produced by Movi to determine whether the read is found in the index or not.
However, Movi implements the classification scheme adopted from the [SPUMONI](https://github.com/oma219/spumoni/tree/main) software written by Omar Y. Ahmed. More details about the classification strategy can be found in the [SPUMONI 2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02958-1).

To use the in-built classification scheme in Movi, you should first build the null database by running:
```
./movi null --index <index directory> --gen-reads
```

Then you should pass the `--classify` flag to the query step:
```
./movi query --index <index directory> --read <reads file> --classify
```
Then, a file called `<reads file>.<index type>.<query type>.report` will be generated with 5 columns:
```
read id: the read name or id
status: whether the read is found or not present
avg max-value: average of the max PML found in each window
above thr: how many windows have a mx PML above the threshold -- the threshold is determined from a null distribution
below thr: how many windows have a mx PML below the threshold
```

## Multi-class classification with Movi Color

Movi Color augments the Movi index with information about the origin of sequences, enabling multi-class and taxonomic classification.

>[Tan, Steven, Sina Majidian, Ben Langmead, and Mohsen Zakeri. "Movi Color: fast and accurate long-read classification with the move structure." bioRxiv (2025)](https://www.biorxiv.org/content/10.1101/2025.05.22.655637v1.abstract)
### Add colors to an existing Movi index

**Note: The following step is not required if the index was built using the `--color` flag.**

```
./movi color --index <index directory>
```
This step requires the file `ref.fa.doc_offsets` to be available in the index directory. Every line in the `ref.fa.doc_offsets` should specify the offset in the concatenated text corresponding to the rigth most end of a document + 1.
It will avoid rebuilding the index from the scratch and only adds the color information to an index.

### Perform Multi-class classification
```
./movi query --index <index directory> --read <reads file> --multi-classify --out-file <output file>
```

The `<output file>` will be in the csv format, each line will include three columns: `read name`, `primary document id`, `secondary document id`. 

For more information about the document IDs and multi-class classification, please see [this page](https://github.com/mohsenzakeri/Movi/wiki/Multi%E2%80%90class-classification-with-Movi-Color).
