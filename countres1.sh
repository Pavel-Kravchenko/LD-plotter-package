echo "/Bash script run successfully/"
home=`pwd`
mask=0.25
in=*.fasta
for x in $in
do
 echo "/Clustalo run/"
 clustalw -INFILE=$x -OUTFILE="$x.aligned.fasta" -TYPE=DNA -OUTPUT=FASTA
 echo "/Done/"
done

a=*.fasta.aligned.fasta
for x in $a
do
 echo $x
 rm -R ${x%.fasta.aligned.fasta}
 mkdir ${x%.fasta.aligned.fasta}
 pwd
 cd ${x%.fasta.aligned.fasta}
 pwd
 python "$home/list_maker.py" ../$x ${x%.fasta.aligned.fasta} $mask
 for o in `cat list_file.txt`
 do
  name=out_2_slices_file_*_$o*.txt
  python "$home/LD_calc.py" $name $o
  Rscript "$home/R_plot.R" --no-save --no-restore --args `pwd` $mask $o
 done
 mv ../${x%.fasta.aligned.fasta}.dnd ./
 mv ../$x ./
 cd ..
 echo "Done"
done
echo "/Bash script end successfully/"



