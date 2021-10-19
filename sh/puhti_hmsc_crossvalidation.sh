#!/bin/bash -l
#SBATCH --account=project_2001175
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=longrun
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000

t1=`date +%s`

# Load r-env-singularity
module load r-env-singularity/4.0.2

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
    sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2001175" >> ~/.Renviron

# Run the R script
srun singularity_wrapper exec Rscript --no-save puhti_hmsc_crossvalidation.R $1 

t2=`date +%s`
tdiff=`echo 'scale=3;('$t2'-'$t1')/3600' | bc`
echo 'FILE '${1}' TOTAL RUNTIME '${tdiff}' HOURS'