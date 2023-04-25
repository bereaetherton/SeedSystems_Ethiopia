#!/bin/sh
#SBATCH --account=epi                # what account you are with
#SBATCH --qos=epi-b                   # use which account 
#SBATCH --job-name=EthINA           # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=betherton@ufl.edu           # Where to send mail	
#SBATCH --nodes=1                   # Use one node
#SBATCH --nodes=1                     # Use one node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Use 1 core
#SBATCH --mem=1GB                   # Memory limit
#SBATCH --time=90:00:00               # Time limit hrs:min:sec
#SBATCH --output=EthOut.out   # Standard output and error log

pwd;hostname;date

module load ufrc
module load R/3.6

Rscript --vanilla NewEthScript.R
date