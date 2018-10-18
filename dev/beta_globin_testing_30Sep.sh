#!/bin/bash

# set up directory structure for outputs
outdir="/Users/STAT2/Dropbox/_Venters_Lab_Resources/3_Rotation_Students/4_Cody/superE"

# testing with 10 reps, fixed WT Bglobin
#Rscript test_superE.R inputs/beta_globin_10_fix.csv $outdir/Bglobin_30Sep18/10rep_fix/test1/ -i 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-10, 10)"

#Rscript test_superE.R inputs/beta_globin_10_fix.csv $outdir/Bglobin_30Sep18/10rep_fix/test2/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-50, 50)"

#Rscript test_superE.R inputs/beta_globin_10_fix.csv $outdir/Bglobin_30Sep18/10rep_fix/test3/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-100, 100)"

#Rscript test_superE.R inputs/beta_globin_10_fix.csv $outdir/Bglobin_30Sep18/10rep_fix/test4/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-150, 150)"

#Rscript test_superE.R inputs/beta_globin_10_fix.csv $outdir/Bglobin_30Sep18/10rep_fix/test5/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-500, 500)"

#Rscript test_superE.R inputs/beta_globin_10_fix.csv $outdir/Bglobin_30Sep18/10rep_fix/test6/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-1000, 1000)"


# testing with 100 reps, fixed WT Bglobin
#Rscript test_superE.R inputs/beta_globin_100_fix.csv $outdir/Bglobin_30Sep18/100rep_fix/test1/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-10, 10)"

#Rscript test_superE.R inputs/beta_globin_100_fix.csv $outdir/Bglobin_30Sep18/100rep_fix/test2/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-50, 50)"

#Rscript test_superE.R inputs/beta_globin_100_fix.csv $outdir/Bglobin_30Sep18/100rep_fix/test3/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-100, 100)"

#Rscript test_superE.R inputs/beta_globin_100_fix.csv $outdir/Bglobin_30Sep18/100rep_fix/test4/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-150, 150)"

Rscript test_superE.R inputs/beta_globin_100_fix.csv $outdir/Bglobin_30Sep18/100rep_fix/test5/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-500, 500)"

Rscript test_superE.R inputs/beta_globin_100_fix.csv $outdir/Bglobin_30Sep18/100rep_fix/test6/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-1000, 1000)"


# testing with 10 reps, normalized WT Bglobin
Rscript test_superE.R inputs/beta_globin_10_norm.csv $outdir/Bglobin_30Sep18/10rep_norm/test1/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-10, 10)"

Rscript test_superE.R inputs/beta_globin_10_norm.csv $outdir/Bglobin_30Sep18/10rep_norm/test2/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-50, 50)"

Rscript test_superE.R inputs/beta_globin_10_norm.csv $outdir/Bglobin_30Sep18/10rep_norm/test3/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-100, 100)"

Rscript test_superE.R inputs/beta_globin_10_norm.csv $outdir/Bglobin_30Sep18/10rep_norm/test4/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-150, 150)"

Rscript test_superE.R inputs/beta_globin_10_norm.csv $outdir/Bglobin_30Sep18/10rep_norm/test5/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-500, 500)"

Rscript test_superE.R inputs/beta_globin_10_norm.csv $outdir/Bglobin_30Sep18/10rep_norm/test6/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-1000, 1000)"


# testing with 100 reps, normalized WT Bglobin
Rscript test_superE.R inputs/beta_globin_100_norm.csv $outdir/Bglobin_30Sep18/100rep_norm/test1/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-10, 10)"

Rscript test_superE.R inputs/beta_globin_100_norm.csv $outdir/Bglobin_30Sep18/100rep_norm/test2/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-50, 50)"

Rscript test_superE.R inputs/beta_globin_100_norm.csv $outdir/Bglobin_30Sep18/100rep_norm/test3/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-100, 100)"

Rscript test_superE.R inputs/beta_globin_100_norm.csv $outdir/Bglobin_30Sep18/100rep_norm/test4/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-150, 150)"

Rscript test_superE.R inputs/beta_globin_100_norm.csv $outdir/Bglobin_30Sep18/100rep_norm/test5/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-500, 500)"

Rscript test_superE.R inputs/beta_globin_100_norm.csv $outdir/Bglobin_30Sep18/100rep_norm/test6/ -i 100 500 1000 2000 5000 10000 20000 -f ~E1+E2+E3+E4+E5+E6 -ab "c(-1000, 1000)"


