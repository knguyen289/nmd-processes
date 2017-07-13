#!/bin/bash
#SBATCH --mail-user=nguyenk@ufl.edu
#SBATCH --job-name=nmdmod
#SBATCH --mail-type=FAIL
#SBATCH --output=/ufrc/ewang/nguyenk/nmd_reg/outlog_mod/nmd_reg_%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 ## threads
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --qos=ewang-b

module load python

python get_mods.py