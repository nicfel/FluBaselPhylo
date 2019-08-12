for i in `seq 1 740`;
do
    bsub -J mat -W "4:00" matlab -nodisplay -nojvm -singleCompThread -r "getNucleotideDifference($i)"
done
