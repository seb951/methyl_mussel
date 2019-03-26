#!/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/r/3.5.0/bin/Rscript

setwd("/scratch/renaut/sophie_breton/Bouvet_Venustaconcha")

#list all sequences
system("ls 001_fastq/*R1.fastq.gz >sequences")
sequences = read.table("sequences", header = F, stringsAsFactors = F)
sequences_2 = sequences
sequences_2[,1] = (gsub("R1","R2",sequences[,1]))

for(i in 1:nrow(sequences))
{
  #gzip extract if neccessary
  if(file.exists(gsub(".gz","",sequences[i,1])) == F)
    {gzip_1 = paste("gzip -cd ",sequences[i,1]," >",gsub(".gz","",sequences[i,1]),sep = "")
    gzip_2 = paste("gzip -cd ",sequences_2[i,1]," >",gsub(".gz","",sequences_2[i,1]),sep = "")
    system(gzip_1)
    system(gzip_2)
    }
  
 #file names
  input_forward.fq = gsub(".gz"," ",sequences[i,1])
  input_reverse.fq = gsub(".gz"," ",sequences_2[i,1])
  output_forward_paired.fq = gsub(".gz","_paired.fq ",sequences[i,1])
  output_forward_unpaired.fq = gsub(".gz","_unpaired.fq ",sequences[i,1])
  output_reverse_paired.fq = gsub(".gz","_paired.fq ",sequences_2[i,1])
  output_reverse_unpaired.fq = gsub(".gz","_unpaired.fq ",sequences_2[i,1])
  seq = strsplit(sequences[i,1],"/")[[1]][2]
  dedup = gsub(".gz","_paired_bismark_bt2_pe.bam",seq)
  methyl_extract = gsub(".gz","_paired_bismark_bt2_pe.deduplicated.bam",seq)
  log = gsub(".gz",".trimmo.log",sequences[i,1])
  
 #trimmomatic
  trimmo = paste("java -jar /home/renaut/bin/trimmomatic-0.36.jar PE -threads 1 -phred33 ",input_forward.fq, input_reverse.fq, output_forward_paired.fq, output_forward_unpaired.fq, output_reverse_paired.fq,output_reverse_unpaired.fq,"ILLUMINACLIP:/home/renaut/trinityrnaseq-Trinity-v2.8.4/trinity-plugins/Trimmomatic-0.36/adapters/all_primers.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:59 2>",log,sep = "")
  
  system(trimmo)
  
  #prep genome
  if(file.exists("00r_reference/Bisulfite_Genome") == F) {
    genome_prep = c("/home/renaut/bismark_v0.20.1/bismark_genome_preparation  --path_to_bowtie /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/bowtie2/2.3.4.1/bin/ --verbose 00r_reference")
    system(genome_prep)
    }

  #Running bismark (Typical alignment example, tolerating one non-bisulfite mismatch per read)
  bismark = paste("/home/renaut/bismark_v0.20.1/bismark --samtools_path /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/samtools/1.9/bin/samtools --path_to_bowtie /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/bowtie2/2.3.4.1/bin/ -n 1 --genome 00r_reference -1 ",output_forward_paired.fq," -2 ",output_reverse_paired.fq,sep = "")
  
  system(bismark)


 # Running deduplicate_bismark
 dedup = paste("/home/renaut/bismark_v0.20.1/deduplicate_bismark -p --samtools_path /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/samtools/1.9/bin/samtools --bam ",dedup,sep = "")
 system(dedup)

 # Running bismark_methylation_extractor (need to  --scaffolds due to highly broken genome draft)
 bismark_methylation_extractor = paste("/home/renaut/bismark_v0.20.1/bismark_methylation_extractor --scaffolds --samtools_path /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/samtools/1.9/bin/samtools --gzip --bedGraph ",methyl_extract,sep = "")
 system(bismark_methylation_extractor)

 #clean up clutter
 #system("rm *methXtractor.temp")

 #cluter in sequencing directory
 system("rm 001_fastq/*unpaired.fq 001_fastq/*trimmo.log 001_fastq/*fastq")  
 
 #move into a specific bismark folder
 mv = paste("mv *",strsplit(seq,"_R1")[[1]][1],"*  003_results/.",sep = "")
 system(mv)
 
 # report
 setwd("003_results/")
 report = paste("/home/renaut/bismark_v0.20.1/bismark2summary",sep = "")
 system(report)
 system("mv bismark_summary_report* ../.")
}

