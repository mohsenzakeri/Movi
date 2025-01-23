# Movi ![GitHub](https://img.shields.io/github/license/mohsenzakeri/movi?color=green)

Movi is a full text index for pangenomes. It takes advantage of the move data structure (by Nishimoto and Tabei). Movi works with genomic file formats. The high locality of the reference and latency hiding in Movi results in low and predictabile query latencies. These properties make Movi ideal for real time applications like Nanopore Adaptive sampling.

>[Zakeri, Mohsen, Brown, Nathaniel K., Ahmed, Omar Y., Gagie, Travis, and Langmead, Ben. "Movi: a fast and cache-efficient full-text pangenome index". iScience (2024)](https://www.cell.com/iscience/fulltext/S2589-0042(24)02691-9)

>Nishimoto, Takaaki, and Yasuo Tabei. "Optimal-time queries on BWT-runs compressed indexes." arXiv preprint arXiv:2006.05104 (2020).

## Install Movi and its dependencies from source


Required dependences: `sdsl`, `zlib`, `cmake`, and `gcc`

To build from source:
```
git clone https://github.com/mohsenzakeri/Movi
cd Movi
mkdir build
cd build
cmake ..
make
```

Building the Movi index on a fasta file requires preprocessing the fasta using the [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds) software. 
So, Movi clones and installs this repository during build. Furthermore, for the `split` and `constant` indexes, Movi uses the [r-permute](https://github.com/drnatebrown/r-permute) implemented by Nate Brown, to
create the splitting bit vector using the strategy introduced by Nishimoto and Tabei.


## Build the Movi index

```
./movi build --fasta <fasta file> --index <index directory> --type <index type>
```

`<fasta file>` is a file including the reference genomes to be indexed.

`<index directory>` is the directory where you want the Movi index to be located.

`<index type` determines the type of Movi index to be built. If no value is passed, the `regular-thresholds` index is built.
Possible index types: `large` `constant` `split` `regular` `regular-thresholds` `blocked` `blocked-thresholds` `tally` `tally-thresholds`


The index will be located at `<index directory>/movi.index`

## Compute Pseudo Matching Lengths (PML) using Movi

To compute PMLs using the movi index, please run the following command on the fastq or fasta file of the reads:
```
./movi query --index <index directory> --read <reads file>
```

`<reads file>` is the address of the fasta or fastq file containing the reads.

After the query command finishes, a file with the same name as the reads file and the extension `<index type>.mls.bin` is generated in the directory that also includes the reads file.
Since this file is in the binary format, to view the PMLs please run the following command:
```
./movi view --mls-file <mls file> | less
```
`<mls file>` is the file generated in the query step.

The output of the last command shows each read's name following by pseudo matching lengths computed for it. A pseudo matching length is outputed for every base of the read. This is the same as the output produced by SPUMONI.

Alternatively, you can directly output the PMLs to a text file using the `--stdout` flag:
```
./movi query --index <index directory> --read <reads file> --stdout > <output file>
```

## Binary classification of the reads

You can apply different classification schemes on the pseudo matching lengths produced by Movi to determine whether the read is found in the index or not.
You can also apply the in-built classification scheme by Movi by passing the flag `--classify` to the query step:
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
The classification scheme is adopted from the [SPUMONI](https://github.com/oma219/spumoni/tree/main) software written by Omar Y. Ahmed.

In order for the classification to work, the null distribution should be also generated in the index building step which is now the default behavior. However, if you want to use  an old index, you can add the null distribution to the index by running:
```
./movi null --index <index directory>
```

### Movi is now published in iScience, you can read more about it here:
> [https://www.cell.com/iscience/fulltext/S2589-0042(24)02691-9](https://www.cell.com/iscience/fulltext/S2589-0042(24)02691-9)

