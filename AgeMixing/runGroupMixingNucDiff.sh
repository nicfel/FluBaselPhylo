for a in `seq 1 1`;
do
    for b in `seq 1 5`;
    do
        for i in `seq 1 100`;
        do
            bsub -J mat -W "4:00" matlab -nodisplay -nojvm -singleCompThread -r "getGroupMixingNucDiff($i,$b)"
        done
    done
done
