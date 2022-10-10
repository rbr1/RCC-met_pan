# #### Run on Oxford compute set up
# # qlogin -pe shmem 12 -P immune-rep.prjc -q short.qc
module purge
module load HDF5/1.10.5-gompi-2019a
module load umap-learn/0.3.10-foss-2019a-Python-3.7.2
module load Seurat/3.1.2-foss-2019a-R-3.6.0
module load Harmony/1.0.0-foss-2019a-R-3.6.0
R

################################ 
library("Seurat")
library('harmony') 
library(ggplot2)
library(pryr)
library(future) 

options(future.globals.maxSize = 30000 * 1024^2)
plan(multiprocess) 

MLtiplet_file = "/users/immune-rep/mfj169/INSTALLED_PROGRAMS/R_MODULES/MLtiplet_2.0/R/R_functions.R"
listFunctions <- function(filename) {
  temp.env <- new.env()
  sys.source(filename, envir = temp.env)
  functions <- lsf.str(envir=temp.env)
  rm(temp.env)
  return(functions)
}
listFunctions(MLtiplet_file)
source(MLtiplet_file)

a = 1
if(a==1){
	file="/users/immune-rep/mfj169/SINGLE_CELL_RNA_SEQ/10X_PIPELINE/Samples_10X_RCC_final.txt"
	p <- as.matrix(read.csv(file, head=T, sep="\t"))
	p=p[which(p[,"To_use_in_RCCpanmet"]=="Yes"),]
	sample_id = as.character(p[,"Sample_Name"])
	sample_output_id = as.character(p[,"Sample_Name"])
	GEX = as.character(p[,"Location_of_GEX"])
	BCR.location = as.character(p[,"Location_of_BCR"])
	TCR.location = as.character(p[,"Location_of_TCR"])
	TCR.location[grep('GU0', sample_id)] = NA ## we are not using these - these samples do not include BCR information
	Overall_sample_group = as.character(p[,"Patient"])
	Site = as.character(p[,"Sample_type"])
	Disease = as.character(p[,"Disease"])
	batch = "RCCmetP_RCC"
	out_dir = "/well/immune-rep/shared/10X_GENOMICS/RCC_WORKING_DATA/"
	out_dir_raw = "/well/immune-rep/shared/10X_GENOMICS/RCC_WORKING_DATA/"
	PLOTS = "FINAL/"
}

################# check files exist
TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
for(c in BCRs){
	if(is.na(BCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = gsub("_annotations.csv",".fasta",BCR.location[c] )
		csv_file = BCR.location[c]
		command = concat(c("cp ", csv_file, " ", out_dir_raw, "filtered_contig_annotations_BCR_", sample,".csv"))
		if(file.exists(fasta_file)==F){print (fasta_file)}
		if(file.exists(fasta_file)==F){print (csv_file)}
}}

for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = gsub("_annotations.csv",".fasta",TCR.location[c] )
		csv_file = TCR.location[c]
		if(file.exists(fasta_file)==F){print (fasta_file)}
		if(file.exists(fasta_file)==F){print (csv_file)}
	}}

TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
incorrect = NULL
for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = gsub("_annotations.csv",".fasta",TCR.location[c] )
		csv_file = TCR.location[c]
		p <- as.matrix(read.csv(csv_file, head=T, sep=","))
		t = table(p[,"chain"])
		if("TRA" %in% names(t)==F){
			print (concat(c(sample, c)))
			print (t)
			incorrect = c(incorrect, c)
		}
		}}
################# filter BCRs and TCRs
VDJAnnotator_10X.py.location = "/well/immune-rep/shared/CODE/VDJANNOTATOR/VDJAnnotator_10X.py"
BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))

# 1: copies raw annotations to central directory (which makes it easier to check that things are run properly, this is an optinoal step)
# 2: copies raw fasta sequences to central directory (which makes it easier to check that things are run properly, this is an optinoal step)
# 3: annotates sequences
for(c in BCRs){
  if(is.na(BCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = gsub("_annotations.csv",".fasta",BCR.location[c] )
		csv_file = BCR.location[c]
		command = concat(c("cp ", csv_file, " ", out_dir_raw, "filtered_contig_annotations_BCR_", sample,".csv"))
		system(command)
		command = concat(c("cp ", fasta_file, " ", out_dir_raw, "filtered_contig_BCR_", sample,".fasta"))
		system(command)
		command = concat(c("/usr/bin/python ", VDJAnnotator_10X.py.location," ", out_dir_raw," ", sample ," ", csv_file," ", fasta_file," IG HOMO_SAPIENS"))
		system(command)
		print(concat(c(c, " ", sample)))
		}
	}

#See guidelines for outputs

TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
TCRs = setdiff(TCRs, grep("GU0", TCR.location))
for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = gsub("_annotations.csv",".fasta",TCR.location[c] )
		csv_file = TCR.location[c]
		command = concat(c("cp ", csv_file, " ", out_dir_raw, "filtered_contig_annotations_TCR_", sample,".csv"))
		system(command)
		command = concat(c("cp ", fasta_file, " ", out_dir_raw, "filtered_contig_TCR_", sample,".fasta"))
		system(command)
		command = concat(c("/usr/bin/python  ", VDJAnnotator_10X.py.location," ", out_dir_raw," ", sample ," ", csv_file ," ", fasta_file," TCR HOMO_SAPIENS"))
		system(command)
		print(concat(c(c, " ", sample)))}
}
#See guidelines for outputs

################# batch for IMGT
library(ape)
outfile = concat(c(out_dir_raw, "All_filtered_contig_BCR_", batch,"_"))
index = 1
total_seqs = NULL
max_n_sequences_per_file = 1000000-10
for(c in BCRs){
	if(is.na(BCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "filtered_contig_BCR_", sample,".fasta"))
		seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
		names(seqs) = paste0(names(seqs),"||", sample, sep = "")
		if( length(seqs) + length(total_seqs)> max_n_sequences_per_file){
			write.dna(total_seqs, concat(c(outfile,index,".fasta")), format = "fasta", append = FALSE)
			index = index+1
			total_seqs = NULL
		}
		total_seqs = c(total_seqs, seqs)
		print(c)
		print (length(total_seqs))
}}
write.dna(total_seqs, concat(c(outfile,index,".fasta")), format = "fasta", append = FALSE,colsep = "")
print(concat(c("scp -p mfj169@rescomp2.well.ox.ac.uk:", outfile,"* ./ " )))


outfile = concat(c(out_dir_raw, "All_filtered_contig_TCR_", batch,"_"))
index = 1
total_seqs = NULL
max_n_sequences_per_file = 1000000-10
for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "filtered_contig_TCR_", sample,".fasta"))
		seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
		names(seqs) = paste0(names(seqs),"||", sample, sep = "")
		if( length(seqs) + length(total_seqs)> max_n_sequences_per_file){
			write.dna(total_seqs, concat(c(outfile,index,".fasta")), format = "fasta", append = FALSE)
			index = index+1
			total_seqs = NULL
		}
		total_seqs = c(total_seqs, seqs)
		print(c)
		print (length(total_seqs))
}}
write.dna(total_seqs, concat(c(outfile,index,".fasta")), format = "fasta", append = FALSE,colsep = "")
print(concat(c("scp -p mfj169@rescomp2.well.ox.ac.uk:", outfile,"* ./ " )))

############## upload to IMGT *************** and return outputs to directory
### unpackage IMGT
dir_IMGT_BCR = "/well/immune-rep/shared/10X_GENOMICS/RCC_WORKING_DATA/FINAL/IMGT_BCR/"
dir_IMGT_TCR = "/well/immune-rep/shared/10X_GENOMICS/RCC_WORKING_DATA/FINAL/IMGT_TCR/"
TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
TCRs = setdiff(TCRs, grep("GU0", TCR.location))
############
type = "BCR"
p <- as.matrix(read.csv(concat(c(dir_IMGT_BCR,"2_IMGT-gapped-nt-sequences.txt")), head=T, sep="\t"))
p1 = p[sort(c(which(p[,"V.D.J.REGION"]!=''),which(p[,"V.J.REGION"]!=''))),]
id = as.character(p1[,"Sequence.ID"])
seq = toupper (as.character(p1[,"V.D.J.REGION"]))
seq = gsub(".","",seq,fixed = T)
seq1 = toupper(as.character(p1[,"V.J.REGION"]))
w = which(seq=='')
seq[w] = seq1[w]
w = which(seq=='')

orig_sample = strsplit(id, "||", fixed = T)
w = setdiff(c(1:length(id)),  grep("||", id, fixed = T))
orig_id = orig_sample
for(i in c(1:length(orig_sample))){
	orig_id[i] = orig_id[[i]][1]
	orig_sample[i] = orig_sample[[i]][2]}
orig_sample = unlist(orig_sample)
orig_sample1 = strsplit(id, "__", fixed = T)
# for(i in c(1:length(w))){
	# orig_sample[w[i]] = orig_sample1[[w[i]]][2]}
orig_sample = unlist(orig_sample)
w = which(is.na(orig_sample))
head(id[w])

orig_sample1 = strsplit(id, "__", fixed = T)
for(i in w){
	orig_sample[i] = orig_sample1[[i]][2]}
unique(orig_sample)
orig_sample [which(orig_sample=="01_Kidney_met_panc1_b")] = "01_Kidney_met_panc1_biopsy"
orig_sample [which(orig_sample=="02_Kidney_met_panc1_b")] = "02_Kidney_met_panc1_biopsy"
orig_sample [which(orig_sample=="04_Kidney_met_panc1_b")] = "04_Kidney_met_panc1_blood"

orig_id = unlist(orig_id)
orig_id = strsplit(orig_id, "__", fixed = T)
for(i in c(1:length(orig_id))){
	orig_id[i] = orig_id[[i]][1]}
orig_id = unlist(orig_id)

orig_sample1 = strsplit(orig_sample,"_CD4")
for(i in c(1:length(orig_sample1))){
	orig_sample1[i] = orig_sample1[[i]][1]}
orig_sample1 = unlist(orig_sample1)

sample_id1 = strsplit(sample_id,"_CD4")
for(i in c(1:length(sample_id1))){
	sample_id1[i] = sample_id1[[i]][1]}
sample_id1 = unlist(sample_id1)

### split by sample
used = NULL
library(ape)
for(c in BCRs){
	if(is.na(BCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "filtered_contig_", type,"_", sample,".fasta"))
		seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
		csv_file = BCR.location[c]
		p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
		seq_names = p2[,"contig_id"]
		w = which(orig_sample1== sample_id1[c])
		w1 = which(orig_id %in% seq_names)
		seq_inds = intersect(w,w1)
		total_seqs = seq[seq_inds]
		names(total_seqs) = orig_id[seq_inds]
		outfile = concat(c(out_dir_raw, "IMGT_filtered_contig_", type,"_", sample,".fasta"))
		write.dna(total_seqs, outfile, format = "fasta", append = FALSE, nbcol = -1,colw = 10000000)
		used = c(used,seq_inds)
		print (c)
}}
length(used)
table(table(used))
length(orig_id)

##############################
type = "TCR"
p <- as.matrix(read.csv(concat(c(dir_IMGT_TCR,"2_IMGT-gapped-nt-sequences.txt")), head=T, sep="\t"))
p1 = p[sort(c(which(p[,"V.D.J.REGION"]!=''),which(p[,"V.J.REGION"]!=''))),]
id = as.character(p1[,"Sequence.ID"])
seq = toupper (as.character(p1[,"V.D.J.REGION"]))
seq = gsub(".","",seq,fixed = T)
seq1 = toupper(as.character(p1[,"V.J.REGION"]))
w = which(seq=='')
seq[w] = seq1[w]
w = which(seq=='')

orig_sample = strsplit(id, "||", fixed = T)
w = setdiff(c(1:length(id)),  grep("||", id, fixed = T))
orig_id = orig_sample
for(i in c(1:length(orig_sample))){
	orig_id[i] = orig_id[[i]][1]
	orig_sample[i] = orig_sample[[i]][2]}
orig_sample = unlist(orig_sample)
orig_sample1 = strsplit(id, "__", fixed = T)
# for(i in c(1:length(w))){
	# orig_sample[w[i]] = orig_sample1[[w[i]]][2]}
orig_sample = unlist(orig_sample)
w = which(is.na(orig_sample))
head(id[w])

orig_sample1 = strsplit(id, "__", fixed = T)
for(i in w){
	orig_sample[i] = orig_sample1[[i]][2]}
unique(orig_sample)
orig_sample [grep("01_Kidney_met", orig_sample)] = "01_Kidney_met_panc1_biopsy"
orig_sample [grep("02_Kidney_met", orig_sample)] = "02_Kidney_met_panc1_biopsy"
orig_sample [grep("04_Kidney_met", orig_sample)] = "04_Kidney_met_panc1_blood"

orig_id = unlist(orig_id)
orig_id = strsplit(orig_id, "__", fixed = T)
for(i in c(1:length(orig_id))){
	orig_id[i] = orig_id[[i]][1]}
orig_id = unlist(orig_id)

orig_sample1 = strsplit(orig_sample,"_CD4")
for(i in c(1:length(orig_sample1))){
	orig_sample1[i] = orig_sample1[[i]][1]}
orig_sample1 = unlist(orig_sample1)

sample_id1 = strsplit(sample_id,"_CD4")
for(i in c(1:length(sample_id1))){
	sample_id1[i] = sample_id1[[i]][1]}
sample_id1 = unlist(sample_id1)


### split by sample
used = NULL
library(ape)
for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "filtered_contig_", type,"_", sample,".fasta"))
		seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
		csv_file = TCR.location[c]
		p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
		seq_names = p2[,"contig_id"]
		w = which(orig_sample1== sample_id1[c])
		w1 = which(orig_id %in% seq_names)
		seq_inds = intersect(w,w1)
		total_seqs = seq[seq_inds]
		names(total_seqs) = orig_id[seq_inds]
		outfile = concat(c(out_dir_raw, "IMGT_filtered_contig_", type,"_", sample,".fasta"))
		write.dna(total_seqs, outfile, format = "fasta", append = FALSE, nbcol = -1,colw = 10000000)
		used = c(used,seq_inds)
		print (c)
}}
length(used)
table(table(used))
length(orig_id)

######################### split the annotation files
############
type = "BCR"
p2 <- as.matrix(read.csv(concat(c(dir_IMGT_BCR,"1_Summary.txt")), head=T, sep="\t"))
p3 = p2[which(p2[,"V.DOMAIN.Functionality"]!="No results"),]
id = as.character(p3[,"Sequence.ID"])

v_mm = p3[,"V.REGION.identity.nt"]
j_mm = p3[,"J.REGION.identity.nt"]
v_mm =strsplit(gsub(" nt","",v_mm),"/",fixed = T)
j_mm =strsplit(gsub(" nt","", j_mm),"/",fixed = T)
for(i in c(1:length(v_mm))){
	v_mm[i] = as.numeric(v_mm[[i]][2])-as.numeric(v_mm[[i]][1])
	j_mm[i] = as.numeric(j_mm[[i]][2])-as.numeric(j_mm[[i]][1])}
v_mm = unlist(v_mm)
names(v_mm) = id
j_mm = unlist(j_mm)
names(j_mm) = id

p <- as.matrix(read.csv(concat(c(dir_IMGT_BCR,"6_Junction.txt")), head=T, sep="\t"))
p1 = p[which(p[,"V.DOMAIN.Functionality"]!="No results"),]
library(stringr)
p1 = p1[which(str_length(p1[,"JUNCTION..AA."])>=3),]
id = as.character(p1[,"Sequence.ID"])
junction_aa = p1[,"JUNCTION..AA."]
junction_nn = p1[,"JUNCTION"]
v = p1[,"V.GENE.and.allele"]
j = p1[,"J.GENE.and.allele"]
v =strsplit(v," ",fixed = T)
j =strsplit(j," ",fixed = T)
for(i in c(1:length(v))){
	v[i] = v[[i]][2]
	j[i] = j[[i]][2]}
v = unlist(v)
j = unlist(j)
v_mms = v_mm[id]
j_mms = j_mm[id]
x1 = cbind(junction_aa, junction_nn,v,j, v_mms, j_mms)

summary(str_length(junction_aa))
rownames(x1) = id

orig_sample = strsplit(id, "||", fixed = T)
w = setdiff(c(1:length(id)),  grep("||", id, fixed = T))
orig_id = orig_sample
for(i in c(1:length(orig_sample))){
	orig_id[i] = orig_id[[i]][1]
	orig_sample[i] = orig_sample[[i]][2]}
orig_sample = unlist(orig_sample)
orig_sample1 = strsplit(id, "__", fixed = T)
# for(i in c(1:length(w))){
	# orig_sample[w[i]] = orig_sample1[[w[i]]][2]}
orig_sample = unlist(orig_sample)
w = which(is.na(orig_sample))
head(id[w])

orig_sample1 = strsplit(id, "__", fixed = T)
for(i in w){
	orig_sample[i] = orig_sample1[[i]][2]}
unique(orig_sample)
orig_sample [grep("01_Kidney_met", orig_sample)] = "01_Kidney_met_panc1_biopsy"
orig_sample [grep("02_Kidney_met", orig_sample)] = "02_Kidney_met_panc1_biopsy"
orig_sample [grep("04_Kidney_met", orig_sample)] = "04_Kidney_met_panc1_blood"

orig_id = unlist(orig_id)
orig_id = strsplit(orig_id, "__", fixed = T)
for(i in c(1:length(orig_id))){
	orig_id[i] = orig_id[[i]][1]}
orig_id = unlist(orig_id)

orig_sample1 = strsplit(orig_sample,"_CD4")
for(i in c(1:length(orig_sample1))){
	orig_sample1[i] = orig_sample1[[i]][1]}
orig_sample1 = unlist(orig_sample1)

sample_id1 = strsplit(sample_id,"_CD4")
for(i in c(1:length(sample_id1))){
	sample_id1[i] = sample_id1[[i]][1]}
sample_id1 = unlist(sample_id1)

### split by sample
used = NULL
library(ape)
for(c in BCRs){
	if(is.na(BCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "filtered_contig_", type,"_", sample,".fasta"))
		seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
		csv_file = BCR.location[c]
		p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
		seq_names = p2[,"contig_id"]
		w = which(orig_sample1== sample_id1[c])
		w1 = which(orig_id %in% seq_names)
		seq_inds = intersect(w,w1)
		x = x1[seq_inds,]
		name = orig_id[seq_inds]
		outfile = concat(c(out_dir_raw, "IMGT_filtered_annotation_", type,"_", sample,".txt"))
		write.table(x, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
		print (c)
}}

########## TCR
type = "TCR"
p2 <- as.matrix(read.csv(concat(c(dir_IMGT_TCR,"1_Summary.txt")), head=T, sep="\t"))
p3 = p2[which(p2[,"V.DOMAIN.Functionality"]!="No results"),]
id = as.character(p3[,"Sequence.ID"])

v_mm = p3[,"V.REGION.identity.nt"]
j_mm = p3[,"J.REGION.identity.nt"]
v_mm =strsplit(gsub(" nt","",v_mm),"/",fixed = T)
j_mm =strsplit(gsub(" nt","", j_mm),"/",fixed = T)
for(i in c(1:length(v_mm))){
	v_mm[i] = as.numeric(v_mm[[i]][2])-as.numeric(v_mm[[i]][1])
	j_mm[i] = as.numeric(j_mm[[i]][2])-as.numeric(j_mm[[i]][1])}
v_mm = unlist(v_mm)
names(v_mm) = id
j_mm = unlist(j_mm)
names(j_mm) = id

p <- as.matrix(read.csv(concat(c(dir_IMGT_TCR,"6_Junction.txt")), head=T, sep="\t"))
p1 = p[which(p[,"V.DOMAIN.Functionality"]!="No results"),]
p1 = p1[which(str_length(p1[,"JUNCTION..AA."])>=3),]
id = as.character(p1[,"Sequence.ID"])
junction_aa = p1[,"JUNCTION..AA."]
junction_nn = p1[,"JUNCTION"]
v = p1[,"V.GENE.and.allele"]
j = p1[,"J.GENE.and.allele"]
v =strsplit(v," ",fixed = T)
j =strsplit(j," ",fixed = T)
for(i in c(1:length(v))){
	v[i] = v[[i]][2]
	j[i] = j[[i]][2]}
v = unlist(v)
j = unlist(j)
v_mms = v_mm[id]
j_mms = j_mm[id]
x1 = cbind(junction_aa, junction_nn,v,j, v_mms, j_mms)
rownames(x1) = id


orig_sample = strsplit(id, "||", fixed = T)
w = setdiff(c(1:length(id)),  grep("||", id, fixed = T))
orig_id = orig_sample
for(i in c(1:length(orig_sample))){
	orig_id[i] = orig_id[[i]][1]
	orig_sample[i] = orig_sample[[i]][2]}
orig_sample = unlist(orig_sample)
orig_sample1 = strsplit(id, "__", fixed = T)
# for(i in c(1:length(w))){
	# orig_sample[w[i]] = orig_sample1[[w[i]]][2]}
orig_sample = unlist(orig_sample)
w = which(is.na(orig_sample))
head(id[w])

orig_sample1 = strsplit(id, "__", fixed = T)
for(i in w){
	orig_sample[i] = orig_sample1[[i]][2]}
unique(orig_sample)
orig_sample [grep("01_Kidney_met", orig_sample)] = "01_Kidney_met_panc1_biopsy"
orig_sample [grep("02_Kidney_met", orig_sample)] = "02_Kidney_met_panc1_biopsy"
orig_sample [grep("04_Kidney_met", orig_sample)] = "04_Kidney_met_panc1_blood"

orig_id = unlist(orig_id)
orig_id = strsplit(orig_id, "__", fixed = T)
for(i in c(1:length(orig_id))){
	orig_id[i] = orig_id[[i]][1]}
orig_id = unlist(orig_id)

orig_sample1 = strsplit(orig_sample,"_CD4")
for(i in c(1:length(orig_sample1))){
	orig_sample1[i] = orig_sample1[[i]][1]}
orig_sample1 = unlist(orig_sample1)

sample_id1 = strsplit(sample_id,"_CD4")
for(i in c(1:length(sample_id1))){
	sample_id1[i] = sample_id1[[i]][1]}
sample_id1 = unlist(sample_id1)

### split by sample
used = NULL
library(ape)
for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "filtered_contig_", type,"_", sample,".fasta"))
		seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
		csv_file = TCR.location[c]
		p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
		seq_names = p2[,"contig_id"]
		w = which(orig_sample1== sample_id1[c])
		w1 = which(orig_id %in% seq_names)
		seq_inds = intersect(w,w1)
		x = x1[seq_inds,]
		name = orig_id[seq_inds]
		outfile = concat(c(out_dir_raw, "IMGT_filtered_annotation_", type,"_", sample,".txt"))
		write.table(x, file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
		used = c(used, seq_inds)
		print (c)
}}

length(used)
table(table(used))
length(orig_id)
################### split amino acid sequence files
############
type = "BCR"
p <- as.matrix(read.csv(concat(c(dir_IMGT_BCR,"5_AA-sequences.txt")), head=T, sep="\t"))
p1 = p[sort(c(which(p[,"V.D.J.REGION"]!=''),which(p[,"V.J.REGION"]!=''))),]
id = as.character(p1[,"Sequence.ID"])
seq = toupper (as.character(p1[,"V.D.J.REGION"]))
seq = gsub(".","",seq,fixed = T)
seq1 = toupper(as.character(p1[,"V.J.REGION"]))
w = which(seq=='')
seq[w] = seq1[w]
w = which(seq=='')

orig_sample = strsplit(id, "||", fixed = T)
w = setdiff(c(1:length(id)),  grep("||", id, fixed = T))
orig_id = orig_sample
for(i in c(1:length(orig_sample))){
	orig_id[i] = orig_id[[i]][1]
	orig_sample[i] = orig_sample[[i]][2]}
orig_sample = unlist(orig_sample)
orig_sample1 = strsplit(id, "__", fixed = T)
# for(i in c(1:length(w))){
	# orig_sample[w[i]] = orig_sample1[[w[i]]][2]}
orig_sample = unlist(orig_sample)
w = which(is.na(orig_sample))
head(id[w])

orig_sample1 = strsplit(id, "__", fixed = T)
for(i in w){
	orig_sample[i] = orig_sample1[[i]][2]}
unique(orig_sample)
orig_sample [grep("01_Kidney_met", orig_sample)] = "01_Kidney_met_panc1_biopsy"
orig_sample [grep("02_Kidney_met", orig_sample)] = "02_Kidney_met_panc1_biopsy"
orig_sample [grep("04_Kidney_met", orig_sample)] = "04_Kidney_met_panc1_blood"

orig_id = unlist(orig_id)
orig_id = strsplit(orig_id, "__", fixed = T)
for(i in c(1:length(orig_id))){
	orig_id[i] = orig_id[[i]][1]}
orig_id = unlist(orig_id)

orig_sample1 = strsplit(orig_sample,"_CD4")
for(i in c(1:length(orig_sample1))){
	orig_sample1[i] = orig_sample1[[i]][1]}
orig_sample1 = unlist(orig_sample1)

sample_id1 = strsplit(sample_id,"_CD4")
for(i in c(1:length(sample_id1))){
	sample_id1[i] = sample_id1[[i]][1]}
sample_id1 = unlist(sample_id1)

### split by sample
used = NULL
library(ape)
for(c in BCRs){
	if(is.na(BCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "filtered_contig_", type,"_", sample,".fasta"))
		seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
		csv_file = BCR.location[c]
		p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
		seq_names = p2[,"contig_id"]
		w = which(orig_sample1== sample_id1[c])
		w1 = which(orig_id %in% seq_names)
		seq_inds = intersect(w,w1)
		total_seqs = seq[seq_inds]
		names(total_seqs) = orig_id[seq_inds]
		outfile = concat(c(out_dir_raw, "IMGT_filtered_amino_acids_", type,"_", sample,".fasta"))
		write.dna(total_seqs, outfile, format = "fasta", append = FALSE, nbcol = -1,colw = 10000000)
		used = c(used,seq_inds)
		print (c)
}}
length(used)
table(table(used))
length(orig_id)

##############################
type = "TCR"
p <- as.matrix(read.csv(concat(c(dir_IMGT_TCR,"5_AA-sequences.txt")), head=T, sep="\t"))
p1 = p[sort(c(which(p[,"V.D.J.REGION"]!=''),which(p[,"V.J.REGION"]!=''))),]
id = as.character(p1[,"Sequence.ID"])
seq = toupper (as.character(p1[,"V.D.J.REGION"]))
seq = gsub(".","",seq,fixed = T)
seq1 = toupper(as.character(p1[,"V.J.REGION"]))
w = which(seq=='')
seq[w] = seq1[w]
w = which(seq=='')

orig_sample = strsplit(id, "||", fixed = T)
w = setdiff(c(1:length(id)),  grep("||", id, fixed = T))
orig_id = orig_sample
for(i in c(1:length(orig_sample))){
	orig_id[i] = orig_id[[i]][1]
	orig_sample[i] = orig_sample[[i]][2]}
orig_sample = unlist(orig_sample)
orig_sample1 = strsplit(id, "__", fixed = T)
# for(i in c(1:length(w))){
	# orig_sample[w[i]] = orig_sample1[[w[i]]][2]}
orig_sample = unlist(orig_sample)
w = which(is.na(orig_sample))
head(id[w])

orig_sample1 = strsplit(id, "__", fixed = T)
for(i in w){
	orig_sample[i] = orig_sample1[[i]][2]}
unique(orig_sample)
orig_sample [grep("01_Kidney_met", orig_sample)] = "01_Kidney_met_panc1_biopsy"
orig_sample [grep("02_Kidney_met", orig_sample)] = "02_Kidney_met_panc1_biopsy"
orig_sample [grep("04_Kidney_met", orig_sample)] = "04_Kidney_met_panc1_blood"

orig_id = unlist(orig_id)
orig_id = strsplit(orig_id, "__", fixed = T)
for(i in c(1:length(orig_id))){
	orig_id[i] = orig_id[[i]][1]}
orig_id = unlist(orig_id)

orig_sample1 = strsplit(orig_sample,"_CD4")
for(i in c(1:length(orig_sample1))){
	orig_sample1[i] = orig_sample1[[i]][1]}
orig_sample1 = unlist(orig_sample1)

sample_id1 = strsplit(sample_id,"_CD4")
for(i in c(1:length(sample_id1))){
	sample_id1[i] = sample_id1[[i]][1]}
sample_id1 = unlist(sample_id1)


### split by sample
used = NULL
library(ape)
for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "filtered_contig_", type,"_", sample,".fasta"))
		seqs <- read.dna(fasta_file, format = "fasta", as.character = T)
		csv_file = TCR.location[c]
		p2 <- as.matrix(read.csv(csv_file, head=T, sep=","))
		seq_names = p2[,"contig_id"]
		w = which(orig_sample1== sample_id1[c])
		w1 = which(orig_id %in% seq_names)
		seq_inds = intersect(w,w1)
		total_seqs = seq[seq_inds]
		names(total_seqs) = orig_id[seq_inds]
		outfile = concat(c(out_dir_raw, "IMGT_filtered_amino_acids_", type,"_", sample,".fasta"))
		write.dna(total_seqs, outfile, format = "fasta", append = FALSE, nbcol = -1,colw = 10000000)
		used = c(used,seq_inds)
		print (c)
}}
length(used)
table(table(used))
length(orig_id)

################# filter BCRs and TCRs on local annotation
VDJAnnotator_10X.py.location = "/well/immune-rep/shared/CODE/VDJANNOTATOR/VDJAnnotator_10X_IMGT.py"
BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
type = "BCR"

'''for(c in BCRs){
	if(is.na(BCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "IMGT_filtered_contig_", type,"_", sample,".fasta"))
		csv_file = BCR.location[c]
		command = concat(c("cp ", csv_file, " ", out_dir_raw, "filtered_contig_annotations_BCR_", sample,".csv"))
		system(command)
		command = concat(c("cp ", fasta_file, " ", out_dir_raw, "filtered_contig_BCR_", sample,".fasta"))
		system(command)
		}}'''
		
for(c in BCRs){
	if(is.na(BCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "IMGT_filtered_contig_", type,"_", sample,".fasta"))
		csv_file = BCR.location[c]
		command = concat(c("/usr/bin/python ", VDJAnnotator_10X.py.location," ", out_dir_raw," ", sample ," ", csv_file," ", fasta_file," IG HOMO_SAPIENS"))
		print(command)
		system(command)
		#print(concat(c(c, " ", sample)))}
	}
}

type = "TCR"
TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
'''for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = gsub("_annotations.csv",".fasta",TCR.location[c] )
		csv_file = TCR.location[c]
		command = concat(c("cp ", csv_file, " ", out_dir_raw, "filtered_contig_annotations_TCR_", sample,".csv"))
		system(command)
		command = concat(c("cp ", fasta_file, " ", out_dir_raw, "filtered_contig_TCR_", sample,".fasta"))
		system(command)
		print(concat(c(c, " ", sample)))}
}'''


for(c in TCRs){
	if(is.na(TCR.location[c])==F){
		sample = sample_id[c]
		fasta_file = concat(c(out_dir_raw, "IMGT_filtered_contig_", type,"_", sample,".fasta"))
		csv_file = TCR.location[c]
		command = concat(c("/usr/bin/python ", VDJAnnotator_10X.py.location," ", out_dir_raw," ", sample ," ", csv_file," ", fasta_file," TCR HOMO_SAPIENS"))
		print(command)
		system(command)
		#print(concat(c(c, " ", sample)))}
	}
}

################# cluster TCRs/BCRs per group
Network_generation.py_location = "/well/immune-rep/shared/CODE/NETWORK_GENERATION/Network_generation_Single_cell_grouped_2.0_IMGT.py"
Overall_sample_groups = unique(Overall_sample_group)
Overall_sample_groups = setdiff(Overall_sample_groups, Overall_sample_groups[grep("GU0", Overall_sample_groups)])
### BCR
BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
run = NULL
for(o in c(1:length(Overall_sample_groups))){
	samples_group = sample_id [intersect(BCRs ,which(Overall_sample_group ==Overall_sample_groups[o]))]
	if(length(samples_group)>0){
		run = c(run, Overall_sample_groups[o])
		#check which have BCRs/TCRs
		sample_fastas = paste0(concat(c(out_dir_raw,"Trimmed_sequences_")), samples_group,"_IG_IGH.fasta" )
		info = file.info(sample_fastas)
		w = which(info$size != 0)
		sample_fastas = sample_fastas[w]
		grouped_sample_id = Overall_sample_groups[o]
		samples_ids = paste0(samples_group[w], collapse = ",")
		sample_fastas = paste0(concat(c(out_dir_raw,"Trimmed_sequences_")), samples_group[w],"_IG_IGH.fasta" )
		sample_fastas = paste0(sample_fastas, collapse = ",")
		cell_info_files = paste0(concat(c(out_dir_raw,"Cell_annotation_")), samples_group[w],"_IG.txt" )
		amino_acid_fastas = paste0(concat(c(out_dir_raw, "IMGT_filtered_amino_acids_BCR_")), samples_group[w],".fasta" , collapse = ",")

		cell_info_files = paste0(cell_info_files, collapse = ",")
		command = concat(c("/usr/bin/python ", Network_generation.py_location," ", out_dir_raw," ", grouped_sample_id ,"_IGH ", sample_fastas ," ", samples_ids ," ", cell_info_files," ", amino_acid_fastas))
		print(command)
		# system(command)
		run = c(run, command)
		command = concat(c("/usr/bin/python ", Network_generation.py_location," ", out_dir_raw," ", grouped_sample_id ,"_IGL ", gsub("_IG_IGH.fasta","_IG_IGL.fasta", sample_fastas) ," ", samples_ids ," ", cell_info_files," ", amino_acid_fastas))
		print(command)
		# system(command)
		run = c(run, command)
		print(concat(c(o, " ", grouped_sample_id)))
	}
}
### TCR
TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
# run = NULL
for(o in c(1:length(Overall_sample_groups))){
	samples_group = sample_id [intersect(TCRs ,which(Overall_sample_group ==Overall_sample_groups[o]))]
	if(length(samples_group)>0){
		run = c(run, Overall_sample_groups[o])
		#check which have BCRs/TCRs
		sample_fastas = paste0(concat(c(out_dir_raw,"Trimmed_sequences_")), samples_group,"_TCR_TRB.fasta" )
		info = file.info(sample_fastas)
		w = which(info$size != 0)
		sample_fastas = sample_fastas[w]
		grouped_sample_id = Overall_sample_groups[o]
		samples_ids = paste0(samples_group[w], collapse = ",")
		sample_fastas = paste0(concat(c(out_dir_raw,"Trimmed_sequences_")), samples_group[w],"_TCR_TRB.fasta" )
		sample_fastas = paste0(sample_fastas, collapse = ",")
		cell_info_files = paste0(concat(c(out_dir_raw,"Cell_annotation_")), samples_group[w],"_TCR.txt" )
		
		amino_acid_fastas = paste0(concat(c(out_dir_raw, "IMGT_filtered_amino_acids_TCR_")), samples_group[w],".fasta" , collapse = ",")

		cell_info_files = paste0(cell_info_files, collapse = ",")
		command = concat(c("python ", Network_generation.py_location," ", out_dir_raw," ", grouped_sample_id ,"_TRB ", sample_fastas ," ", samples_ids ," ", cell_info_files," ", amino_acid_fastas))
		print(command)
		run = c(run, command)
		command = concat(c("python ", Network_generation.py_location," ", out_dir_raw," ", grouped_sample_id ,"_TRA ", gsub("_TCR_TRB.fasta","_TCR_TRA.fasta", sample_fastas) ," ", samples_ids ," ", cell_info_files," ", amino_acid_fastas))
		print(command)
		run = c(run, command)
		print(concat(c(o, " ", grouped_sample_id)))
	}
}	

outfile = concat(c(out_dir_raw, "VDJ_commands.txt"))
write.table(cbind(run), file = outfile, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,col.names = F, qmethod = c("escape", "double"),fileEncoding = "")
print(concat(c("Run commands from file
cat ", outfile)))
################# Get VDJ annotation information object

BCRs = unique(sort(intersect(which(is.na(BCR.location)==F), which(BCR.location!=''))))
TCRs = unique(sort(intersect(which(is.na(TCR.location)==F), which(TCR.location!=''))))
TCRs = setdiff(TCRs, grep("GU0", TCR.location))

labels.vector =c() ### get all cell IDs
for(c in c(1:length(sample_id))){
	if(c %in% BCRs){
		sample = sample_id[c]
		cell_info_BCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_IG.txt" )
		p1 <- as.matrix(read.csv(cell_info_BCR, head=TRUE, sep="\t"))
		cells = gsub("-1",concat(c("||", sample)), p1[,"X.cell"])
		labels.vector = c(labels.vector, cells)}
	if(c %in% TCRs){
		sample = sample_id[c]
		cell_info_TCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_TCR.txt" )
		p1 <- as.matrix(read.csv(cell_info_TCR, head=TRUE, sep="\t"))
		cells = gsub("-1",concat(c("||", sample)), p1[,"X.cell"])
		labels.vector = c(labels.vector, cells)
	}
	print (c)
}
labels.vector = sort(unique(labels.vector))

headers = c( "X.cell","contig1","chain1","constant_region1","n_umis1","V_gene_10X1",
"J_gene_10X1","cdr3_aa1","cdr3_nn1","V_gene1","J_gene1","V_mm1" , "J_mm1","mixed_contig_chain1","mixed_contig_n_umis1","contig2", "chain2", "constant_region2",
"n_umis2","V_gene2","J_gene2", "cdr3_aa2","cdr3_nn2","V_gene2.1", "J_gene2.1","V_mm2","J_mm2","mixed_contig_chain2","mixed_contig_n_umis2")
m_VDJ_BCR = matrix(data = "-", nrow = length(labels.vector), ncol = length(headers), dimnames = c(list(labels.vector), list(headers)))

m_VDJ_TCR = matrix(data = "-", nrow = length(labels.vector), ncol = length(headers), dimnames = c(list(labels.vector), list(headers)))

for(c in c(1:length(BCRs))){
	sample = sample_id[BCRs[c]]
	cell_info_TCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_TCR.txt" )
	cell_info_BCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_IG.txt" )
	p1 <- as.matrix(read.csv(cell_info_BCR, head=TRUE, sep="\t"))
	cells = gsub("-1",concat(c("||", sample)), p1[,"X.cell"])
	w = which(cells %in% labels.vector)
	m_VDJ_BCR[cells[w], ] = p1[w,headers]
	print(concat(c(c," ",sample_id[c])))
	print(table(m_VDJ_BCR[,"chain1"]))
}

for(c in c(1:length(TCRs))){
	sample = sample_id[TCRs[c]]
	cell_info_TCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_TCR.txt" )
	cell_info_BCR = paste0(concat(c(out_dir_raw,"Cell_annotation_")), sample,"_IG.txt" )
	p1 <- as.matrix(read.csv(cell_info_TCR, head=TRUE, sep="\t"))
	cells = gsub("-1",concat(c("||", sample)), p1[,"X.cell"])
	w = which(cells %in% labels.vector)
	m_VDJ_TCR[cells[w], ] = p1[w,headers]
	print(concat(c(c," ",sample_id[c])))
	print(table(m_VDJ_TCR[,"chain1"]))
}

clone1 = rep("-", length(labels.vector))
clone2 = rep("-", length(labels.vector))
m_VDJ_BCR = cbind(m_VDJ_BCR, clone1, clone2)
m_VDJ_TCR = cbind(m_VDJ_TCR, clone1, clone2)

run = NULL
for(o in c(1:length(Overall_sample_groups))){
	cluster_file1 = paste0(concat(c(out_dir_raw,"Cluster_identities_")), Overall_sample_groups[o],"_TRA.txt" )
	cluster_file2 = paste0(concat(c(out_dir_raw,"Cluster_identities_")), Overall_sample_groups[o],"_TRB.txt" )
cluster_files = c(cluster_file1, cluster_file2)
	types = c("TRA","TRB")
	info = file.info(cluster_files)
	w = which(info$size != 0)
	if(length(w)>0){
		for(c in c(1:length(cluster_files))){
			p1 <- as.matrix(read.csv(cluster_files[c], head=F, sep="\t"))
			p1=p1[-1,]
			cells = p1[,3]#gsub("-1",concat(c("||", sample_output_id[c])), p1[,3])
			clone = p1[,2]
			clone = paste (as.numeric(clone), concat(c("||",Overall_sample_groups[o])))
			clone = gsub(" ", "", clone)
			cells = strsplit(cells,"||", fixed = T)
			cell=NULL
			sample_source = NULL
			for(i in c(1:length(cells))){
				cell = c(cell, cells[[i]][1])
				sample_source = c(sample_source, cells[[i]][2])
			}
			m = sample_id [match(sample_source, sample_id)]
			cells = apply(cbind(cell, m), 1, paste, collapse="||")
			names(clone) = cells
			w = which(cells %in% labels.vector==T)
			cells = cells[w]
			if(types[c] %in% c("TRA", "IGH")){
				ids_clones = names(which(m_VDJ_TCR[cells, "chain1"]== types[c]))
				m_VDJ_TCR[ids_clones,"clone1"] = clone[ids_clones]
			}else{
				chain = m_VDJ_TCR[cells, "chain2"]
				chain[which(chain=="IGK")] = "IGL"
				ids_clones = names(which(chain == types[c]))
				m_VDJ_TCR[ids_clones,"clone2"] = clone[ids_clones]
			}
		}
	}
	print(length(which(m_VDJ_TCR[,"clone1"]!='-')))
}
	
for(o in c(1:length(Overall_sample_groups))){
	cluster_file3 = paste0(concat(c(out_dir_raw,"Cluster_identities_")), Overall_sample_groups[o],"_IGH.txt" )
	cluster_file4 = paste0(concat(c(out_dir_raw,"Cluster_identities_")), Overall_sample_groups[o],"_IGL.txt" )
	cluster_files = c( cluster_file3, cluster_file4)
	types = c("IGH","IGL")
	info = file.info(cluster_files)
	w = which(info$size != 0)
	if(length(w)>0){
		for(c in c(1:length(cluster_files))){
			p1 <- as.matrix(read.csv(cluster_files[c], head=F, sep="\t"))
			p1=p1[-1,]
			cells = p1[,3]
			clone = p1[,2]
			clone = paste (as.numeric(clone), concat(c("||",Overall_sample_groups[o])))
			clone = gsub(" ", "", clone)
			cells = strsplit(cells,"||", fixed = T)
			cell=NULL
			sample_source = NULL
			for(i in c(1:length(cells))){
				cell = c(cell, cells[[i]][1])
				sample_source = c(sample_source, cells[[i]][2])
			}
			m = sample_id [match(sample_source, sample_id)]
			cells = apply(cbind(cell, m), 1, paste, collapse="||")
			names(clone) = cells
			w = which(cells %in% labels.vector==T)
			cells = cells[w]
			if(types[c] %in% c("TRA", "IGH")){
				ids_clones = names(which(m_VDJ_BCR[cells, "chain1"]== types[c]))
				m_VDJ_BCR[ids_clones,"clone1"] = clone[ids_clones]
			}else{
				chain = m_VDJ_BCR[cells, "chain2"]
				chain[which(chain=="IGK")] = "IGL"
				ids_clones = names(which(chain == types[c]))
				m_VDJ_BCR[ids_clones,"clone2"] = clone[ids_clones]
			}
		}
	}
	print(length(which(m_VDJ_BCR[,"clone1"]!='-')))
}


m_VDJ_BCR[which(m_VDJ_BCR[,"constant_region1"]=="None"),"constant_region1"] = '-'
m_VDJ_BCR[which(m_VDJ_BCR[,"constant_region2"]=="None"),"constant_region2"] = '-'
m_VDJ_TCR[which(m_VDJ_TCR[,"constant_region1"]=="None"),"constant_region1"] = '-'
m_VDJ_TCR[which(m_VDJ_TCR[,"constant_region2"]=="None"),"constant_region2"] = '-'

VDJ_object = c(list(m_VDJ_BCR),list(m_VDJ_TCR))
names(VDJ_object) = c("BCR","TCR")

saveRDS(file=concat(c(out_dir_raw,"/VDJ_information_", batch,".VDJ")), VDJ_object)

## output as matrix file
out_file_table = concat(c(out_dir_raw,"/VDJ_information_BCR_", batch,".txt"))
write.table(m_VDJ_BCR, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
out_file_table = concat(c(out_dir_raw,"/VDJ_information_TCR_", batch,".txt"))
write.table(m_VDJ_TCR, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,col.names = T, qmethod = c("escape", "double"),fileEncoding = "")
print(concat(c("scp -p mfj169@rescomp2.well.ox.ac.uk:", out_file_table," ./ " )))



