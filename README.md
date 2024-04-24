# Movi ![GitHub](https://img.shields.io/github/license/mohsenzakeri/movi?color=green)

Movi is a full text index for pangenomes. It takes advantage of the move data structure (by Nishimoto and Tabei). Movi works with genomic file formats. The high locality of the reference and latency hiding in Movi results in low and predictabile query latencies. These properties make Movi ideal for real time applications like Nanopore Adaptive sampling.

>[Zakeri, Mohsen, Brown, Nathaniel K., Ahmed, Omar Y., Gagie, Travis, and Langmead, Ben. "Movi: a fast and cache-efficient full-text pangenome index". bioRxiv preprint (2023)](https://www.biorxiv.org/content/10.1101/2023.11.04.565615v2)

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
Please download and install [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds) before proceeding to the next step.

After installing `pfp-thresholds`, please edit the `preprocess_ref.sh` script in the main directory of Movi to include the path to the `pfp-thresholds` binary:
```
pfp=<PATH TO PFP_THRESHOLDS BINARY>
```

Please also make sure the `movi_default`, `movi_constant` and `prepare_ref` variables include the correct paths to the corresponding binaries according to where the binaries (in the `build` directory after the installation) are located on your system:
```
movi_default=<PATH TO MOVI-DEFAULT BINARY>
movi_constant=<PATH TO MOVI-CONSTANT BINARY>
prepare_ref=<PATH TO PREPARE_REF BINARY>
```

## Build the Movi index

### default index
After including the correct paths in the script, to build the `default` Movi index, you may run the build-script with the following command:
```
bash preprocess_ref.sh default <fasta list file> <index directory>
```
`<fasta list file>` is a file which contains the address of all the reference fasta files to be indexed with Movi. Each line in the list file should be the address of a separate fasta file.

`<index directory>` is the directory where you want the Movi index to be located.

The index will be located at `<index directory>/movi_index.bin`

### constant index
Building the constant index requires further preprocessing the reference file using the [r-permute](https://github.com/drnatebrown/r-permute) tool. 
Please download and install [r-permute](https://github.com/drnatebrown/r-permute) before proceeding to the next step.

After installing `r-permute`, please edit the `preprocess_ref.sh` script in the main directory of Movi to include the paths to the `build_constructor` and `run_construct` binaries:
```
bconstructor=<PATH TO BUILD_CONSTRUCTOR BINARY>
rconstructor=<PATH TO RUN_CONSTRUCTOR BINARY>
```
`build_constructor` and `run_constructor` binaries are found at `r-permute/build/test/src/` after building the `r-permute` tool.

**To build the `constant` Movi index, please run the following:**
```
bash preprocess_ref.sh constant <fasta list file> <index directory>
```
The constant index will be located at `<index directory>/constant_index/movi_index.bin`

## Compute Pseudo Matching Lengths (PML) using Movi

To compute PMLs using the **default** movi index, please run the following command on the fastq or fasta file of the reads:
```
movi-default query --pml --index <default index directory> --read <reads file>
```
or the following command for using the **constant** movi index:
```
movi-constant query --pml --index <constant index directory> --read <reads file>
```

`<reads file>` is the address of the fasta or fastq file containing the reads.

After the query command finishes, a file with the same name as the reads file and the extension `mpml.bin` is generated in the directory that also includes the reads file.
Since this file is in the binary format, to view the PMLs please run the following command:
```
movi-default view --pml-file <mpml file> | less
```
`<mpml file>` is the file generated in the query step.

The output of the last command shows each read's name following by pseudo matching lengths computed for it. A pseudo matching length is outputed for every base of the read. This is the same as the output produced by SPUMONI.

### You can read more about Movi here:
> [https://www.biorxiv.org/content/10.1101/2023.11.04.565615v2](https://www.biorxiv.org/content/10.1101/2023.11.04.565615v2)

