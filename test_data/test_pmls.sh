if [ ! -z "$3" ]; then
    expected=$3
    echo "A"
else
    expected="../test_data/expected_pmls"
    echo "B"
fi
index_dir=$1
reads="../test_data/reads.fasta"
movi_dir="../build"
$movi_dir/movi-default query --pml --index $index_dir --read $reads
$movi_dir/movi-default view --pml-file $reads.default.mpml.bin > test_out
ls -l ../test_data/*.bin
diff $expected test_out | wc -l


index_dir=$2
if [ ! -z "$index_dir" ]; then
    $movi_dir/movi-constant query --pml --index $index_dir --read $reads
    $movi_dir/movi-default view --pml-file $reads.constant.mpml.bin > test_out
    ls -l ../test_data/*.bin
    diff ../test_data/expected_pmls test_out | wc -l
fi
