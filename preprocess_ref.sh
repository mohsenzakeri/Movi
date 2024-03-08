pfp="<PATH TO PFP_THRESHOLDS BINARY>"

prepare_ref="<PATH TO PREPARE_REF BINARY>"
movi_default="<PATH TO MOVI-DEFAULT BINARY>"
movi_constant="<PATH TO MOVI-CONSTANT BINARY>"

bconstructor="<PATH TO BUILD_CONSTRUCTOR BINARY>"
rconstructor="<PATH TO RUN_CONSTRUCTOR BINARY>"

t="/usr/bin/time -o "

fasta_list="$2"
index_dir="$3"
clean_fasta="$index_dir/ref.fa"

mkdir -p ${index_dir}

# check if index is already built
if [ "$1" == "default" ]; then
  if test -f "$index_dir/movi_index.bin"
  then
    echo "The index is already built at $index_dir/movi_index.bin"
    exit
  fi
fi

# check if the constant index is already built
if [ "$1" == "constant" ]; then
  if test -f "$index_dir/constant_index/movi_index.bin"
  then
    echo "The constant index is already built at $index_dir/constant_index/movi_index.bin"
    exit
  fi
fi


# check if the clean fasta is already built, otherwise, build
if test -f "$index_dir/ref.fa"
then
  echo "The clean_fasta is already made, skipping.."
else
  cmd="$t $index_dir/build.prepare_ref.time $prepare_ref $fasta_list $clean_fasta list"
  echo $cmd
  eval $cmd
fi

# check if the bwt and thresholds files are already built, otherwise, build
if test -f "$index_dir/ref.fa.thr_pos"
then
  if test -f "$index_dir/ref.fa.bwt" 
  then
    echo "The pfp_threshold step is already done, skipping.."
  else
    cmd="$t $index_dir/build.pfp.time $pfp -f $clean_fasta"
    echo $cmd
    eval $cmd
  fi
else
  cmd="$t $index_dir/build.pfp.time $pfp -f $clean_fasta"
  echo $cmd
  eval $cmd
fi


if [ "$1" == "constant" ]; then
  cmd="$t $index_dir/build_rlbwt.time $movi_constant rlbwt --bwt-file $clean_fasta.bwt"
  echo $cmd
  eval $cmd
  cmd="$t $index_dir/build_constructor.time $bconstructor $clean_fasta"
  echo $cmd
  eval $cmd
  cmd="$t $index_dir/run_constructor.time $rconstructor $clean_fasta -d 5"
  echo $cmd
  eval $cmd
  cmd="$t $index_dir/build.movi.time $movi_constant build --fasta $clean_fasta --index $index_dir/constant_index"
  echo $cmd
  eval $cmd

  echo "Removing extra files.."
  rm $index_dir/ref.fa*
  echo "done"
fi



if [ "$1" == "default" ]; then
  cmd="$t $index_dir/build.movi.time $movi_default build --fasta $clean_fasta --index $index_dir"
  echo $cmd
  eval $cmd

  echo "Removing extra files.."
  rm $index_dir/ref.fa*
  echo "done"
fi
