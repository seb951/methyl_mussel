#!/bin/bash

#SBATCH --time=0-00:09
#SBATCH --mem=1G
#SBATCH --account=def-bruneaua
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1


cd $SCRATCH/sophie_breton/Bouvet_Venustaconcha/00p_pipeline/

rm ../*0.methXtractor.temp
#rm ../*1.methXtractor.temp
#rm ../*2.methXtractor.temp
#rm ../*3.methXtractor.temp
#rm ../*4.methXtractor.temp
#rm ../*5.methXtractor.temp
#rm ../*6.methXtractor.temp
#rm ../*8.methXtractor.temp
#rm ../*9.methXtractor.temp


