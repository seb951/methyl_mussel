
setwd("/Users/jerry/Dropbox/CSBQ/sophie_breton/methylation_2019")

library(dplyr)
library(ggplot2)


#Length of sequence
len = read.table("data/assembly_k41_B24G_H4_kc3-10.fa1000_nocontam.len",stringsAsFactors = F)
len_df = data.frame(name = as.numeric(gsub(">","",len[seq(1,nrow(len),by = 2),1])), len = as.numeric(len[seq(2,nrow(len),by = 2),1],stringsAsFactors = F))
#len = len[seq(2,nrow(len),by = 2),1]
#len_o = log(len[order(len)],10)

#coverage file
cov = read.table("data/Vell_8b_R1.fastq_paired_bismark_bt2_pe.deduplicated.bismark.cov")
colnames(cov) = c("contig","start","stop","percent","methyl","non_methyl")

cov[,1] = as.numeric(cov[,1]);cov[,2] = as.numeric(cov[,2]);cov[,3] = as.numeric(cov[,3])

#add contig size and reorder according to it.
cov$contig_size = len_df[match(cov[,1],len_df[,1]),2]
cov = cov[order(cov$contig_size),]


#sliding window of 1000 sites, by sampling only 0.1% of the data to make object smaller
cov$sliding_mean = NA
for(i in seq(1000,nrow(cov)-1000,by = 1000)) #A sequence every 1,000 site...
  {
    cov$sliding_mean[i] = mean(cov[(i-1000):(i+1000),4],na.rm = T) #mean
    
    if(i %% 500000 == 0) print(paste("Done ",i," of: ",nrow(cov),", The time is: ",Sys.time(),sep = ""))
    }


###
#methyl sites summary graphs...
###

system("wget https://raw.githubusercontent.com/seb951/methyl_mussel/master/bismark_summary_report.txt")
  
bismark_summary_report = read.table("bismark_summary_report.txt",header = T, stringsAsFactors = F,sep = "\t")

bismark_summary_report$CPG_mfract =  bismark_summary_report[,10] / rowSums(bismark_summary_report[,10:11]) * 100
bismark_summary_report$CPG_fract =  bismark_summary_report[,11] / rowSums(bismark_summary_report[,10:11]) * 100

bismark_summary_report$CHG_mfract =  bismark_summary_report[,12] / rowSums(bismark_summary_report[,12:13]) * 100
bismark_summary_report$CHG_fract =  bismark_summary_report[,13] / rowSums(bismark_summary_report[,12:13]) * 100

bismark_summary_report$CHH_mfract =  bismark_summary_report[,14] / rowSums(bismark_summary_report[,14:15]) * 100
bismark_summary_report$CHH_fract =  bismark_summary_report[,15] / rowSums(bismark_summary_report[,14:15]) * 100

#data in data.frame format
bismark_df = data.frame(sites = unlist(bismark_summary_report[,c(10:15)]))

bismark_df$fraction = unlist(bismark_summary_report[,c(16:21)])

bismark_df$type = rep(c(rep("methyl",6),rep("unmethyl",6)),3)

bismark_df$position = c(rep("CPG",12),rep("CHG",12),rep("CHH",12))

bismark_df$individuals = rep(sapply(strsplit(bismark_summary_report[,1],split = "_R1"),'[',1),6)

bismark_df$sites_log = log(bismark_df$sites,2)

p=ggplot() + labs(title = "Methylation", x ="Individuals", y = "Number of sites") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
    geom_bar(aes(y = sites, x = individuals,fill = interaction(type,position)),
    data = bismark_df,stat="identity")

dev.new(width = 8,height = 4,units = "cm",noRStudioGD = TRUE)
p + facet_grid(rows=vars(position),scales = "free") + scale_fill_manual(labels = rep(c("methylated","unmethylated"),3),name = "Sites",values=c("#0000BB55", "#0000BBFF","#BB837855","#BB8378FF","#00640055","#006400FF"))
dev.print(device=pdf,"figures/methylation.pdf", onefile=FALSE)
dev.off()

dev.new()
png("figures/methylation.png",width=8, res =400,height=4,units= 'in')
p + facet_grid(rows=vars(position),scales = "free") + scale_fill_manual(labels = rep(c("methylated","unmethylated"),3),name = "Sites", values=c("#0000BB55", "#0000BBFF","#BB837855","#BB8378FF","#00640055","#006400FF"))
dev.off()

###
###Does length of sequence influence the methylation?
###
#summarise by contig
by_cov = cov %>% group_by(contig) %>% summarise_all(list(sum,mean))

#adding a total percentage
#by_cov$percent = by_cov[,5] / rowSums(by_cov[,5:6]) * 100

#data.frame (name,  pas rapport...,mean contig size,methyl / non -methyl sub)
by_cov = as.data.frame(by_cov[,c(1,10,14,5,6)])

#now calculate a mean methylation per size class (1,300000 )
methyl_size_class = data.frame(size_class = seq(1000,320000,by =1000), mean_methyl = 0, total_methyl_per_contig= 0, total_sites_per_contig = 0,total_methyl= 0, total_sites = 0)

for(m in 2:nrow(methyl_size_class))
{
  temp = by_cov[by_cov[,3]>methyl_size_class[m-1,1] & by_cov[,3]<methyl_size_class[m,1],4:5]
  if(nrow(temp)>0) methyl_size_class[m,2] = sum(temp[,1])/sum(temp[,1:2])
  if(nrow(temp)>0) methyl_size_class[m,3] = sum(temp[,1]) / nrow(temp)
  if(nrow(temp)>0) methyl_size_class[m,4] = sum(temp[,1:2]) / nrow(temp)
  if(nrow(temp)>0) methyl_size_class[m,5] = sum(temp[,1])
  if(nrow(temp)>0) methyl_size_class[m,6] = sum(temp[,2])
  }


#
#total sites
png("figures/methylation_contig_size.png",width=8, res =400,height=5,units= 'in')
#dev.new(width=8, height=5,units = "cm",noRStudioGD = TRUE)
par(mar=c(6,6,4,6),mgp = c(3.5,2,1))
plot(methyl_size_class[methyl_size_class[,2]!=0,1],log(methyl_size_class[methyl_size_class[,2]!=0,6],2),
     type = "l",pch = 20,lwd = 6,col = "#000000FF",ylim = c(0,18),yaxt= "n",xaxt= "n",
     xlab = "contig size (X 1,000 bp)",ylab = "Number of CPG sites",main = "Relationship between contig size and methylation")

#fraction of sites methylated per size class
points(methyl_size_class[methyl_size_class[,2]!=0,1],methyl_size_class[methyl_size_class[,2]!=0,2]*log(250000,2),col = "darkred",pch = 20)

points(methyl_size_class[methyl_size_class[,2]!=0,1],log(methyl_size_class[methyl_size_class[,2]!=0,5],2),type = "l",pch = 20,lwd = 6,col = "#00000088")
mean = sum(cov[,5])/sum(cov[,5:6])
at = c(0,mean*log(250000,2),log(250000,2)/4,log(250000,2)/2,log(250000,2)/4*3,log(250000,2))
label = c(1,10,100,1000,10000,50000,250000)
label2 = c(1,10,100,"1k","10k","50k","250k")
label3 = c(0,signif(mean*100,4),25,50,75,100)
axis(2,at = log(label,2), label = label2,lwd =4)
axis(4,at = at,label = label3,lwd =4,col = "darkred",col.axis="darkred",col.lab = "darkred")
mtext("fraction of sites methylated",4,line =3.5,col = "darkred")
text(x=360000,y=2,label = signif(mean*100,3), xpd=T,col = "darkred")


axis(1,at = c(0,50,100,150,200,250,300)*1000, label = c(0,50,100,150,200,250,300),lwd  =4)


legend(x=30000,y=18,legend = c("Unmethylated","Methylated"),cex =0.8,fill = c("#000000FF","#00000088"),bg = "#FFFFFF88")

abline(h = mean*log(250000,2),col = "#8B0000AA",lwd = 2,lty = 2)  
  
#dev.print(device=pdf,"figures/methylation_contig_size.pdf", onefile=FALSE)
dev.off()






###
###sandbox
###
by_cov_sum = mutate(by_cov_10k,sum = rle(cov[,1])$lengths)

by_cov_sum_10 = by_cov_sum[by_cov_sum$sum>1,]

by_cov_sum_10$pos = sample(c(1:nrow(by_cov_sum_10)))

plot(by_cov_sum_10$pos, by_cov_sum_10$V4,pch = 20,lwd = 0.1)

by_cov$fraction = by_cov$methyl/rowSums(by_cov[,5:6])* 100

by_cov_sub = as.data.frame(by_cov[seq(1,193000,by = 100),])

plot(log(by_cov_sub$contig_size,10),by_cov_sub$fraction)



  

#cov_sd_plot = cov #keep all data...
cov_sd_plot = cov[!is.na(cov$sliding_mean),] #keep only 0.1% of data...

rownames(cov_sd_plot) = 1:nrow(cov_sd_plot)


#plot the methylation percentage.
dev.new()
plot(log(cov_sd_plot$contig_size,10),cov_sd_plot[,4],pch = 20, col = "#99999950",
     main = "Methylation",
     xlab = "Contig size (log10)",  
     ylab = "Methylation fraction per site")
#axis(side = 4,at = c(1,10,100),label = c(1,10,100),col = "blue")

#add the trend line (mean per 1,000 sites)
points(log(cov_sd_plot$contig_size,10),cov_sd_plot$sliding_mean,type= "l",col = "darkblue",lwd = 1)




