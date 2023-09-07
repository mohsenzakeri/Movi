pfp="~/pfp-thresholds/build/pfp_thresholds"

t="/usr/bin/time -v -o"
fasta="$1"
index_dir="$2"
mkdir ${index_dir}

cmd="$t $index_dir/build.pfp.time $pfp -f $fasta"
echo $cmd
eval $cmd
# cmd="$t $index_dir/build.movi.time $marlin build reg $fasta $index_dir"
# echo $cmd
# eval $cmd

