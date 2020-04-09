# https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
# https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/

### reference of mouse for fasta file of CDS of "Mus_musculus.GRCm38.cds.all.fa.gz"
## as well as the gft file "Mus_musculus.GRCm38.98.gtf",we collect them via ENSEMBL ftp:
## https://useast.ensembl.org/info/data/ftp/index.html

### download only gft obtained by Salmon

#load library 
require("biomaRt")
library("tximport")
library("readr")
library("GenomicFeatures")

folder_salmon_lujambio_original <- "/LujambioLab/data/salmon/quant_all_original/"
folder_salmon_lujambio <- "/LujambioLab/data/salmon/quant_all/"

### rewrite the quant.sf but need to remove the version of transcript
files_quant<-list.files(path = folder_salmon_lujambio_original, pattern = "_quant.sf", all.files = TRUE,
                        full.names = FALSE, recursive = TRUE)
for(i in 1:length(files_quant)){
  quant <- read.table(file=paste0(folder_salmon_lujambio_original,files_quant[i]),sep="\t",header=TRUE)
  quant$Name <-gsub("\\.[0-9]*","",quant$Name)
  write.table(quant,file=paste0(folder_salmon_lujambio,files_quant[i]),sep="\t",row.names = FALSE,
              quote=FALSE)
}

folder_salmon_ref="/LujambioLab/data/salmon/"
folder_salmon_lujambio_final <- "/LujambioLab/data/salmon/quant_final/"

###Connection to BioMart ENSEMBL
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", 
               path="/biomart/martservice" ,dataset="mmusculus_gene_ensembl")

gffFile <- paste0(folder_salmon_ref,"Mus_musculus.GRCm38.98.gtf")

txdb <- makeTxDbFromGFF(file=gffFile, format="gtf", dataSource="ENSEMBL", organism="Mus musculus")

#list of chromosome
seqlevels(txdb)
#colnomes of data
columns(txdb)
#list of transcrip
GR <- transcripts(txdb)

#reading this tx2gene data.frame can be accomplished from a TxDb object
#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

all_samples <- c()
file_salmon<-list.files(path=folder_salmon_lujambio,pattern = "_quant.sf")
name_file_salmon <- gsub("_quant.sf","",file_salmon)
files <- file.path(folder_salmon_lujambio, file_salmon)
names(files) <- name_file_salmon

#### Creation table with Mouse and Human genes for Reads count
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

expression.Read_count <-as.data.frame(txi$counts)

write.table(expression.Read_count,
            file=paste0(folder_salmon_lujambio_final,"all_samples_mouse_genes_read_count.tsv"),
            sep="\t",quote=FALSE,row.names = FALSE)

#### give the gene as human name
musGenes_ensemb <- rownames(expression.Read_count)

# Basic function to convert mouse to human gene names
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("mgi_symbol","ensembl_gene_id"), filters = "ensembl_gene_id",
                 values = musGenes_ensemb , mart = mouse, 
                   attributesL = c("hgnc_symbol","ensembl_gene_id","external_gene_name"), 
                 martL = human, uniqueRows=T)

colnames(genesV2)[2]<-"Gene.stable.ID_mouse"
colnames(genesV2)[4]<-"Gene.stable.ID_human"

expression.Read_count$Gene.stable.ID_mouse <- rownames(expression.Read_count)
dim(expression.Read_count)
dim(genesV2)
expression.Read_count_human <- merge(genesV2,expression.Read_count,by="Gene.stable.ID_mouse")
dim(expression.Read_count_human)

write.table(expression.Read_count_human,
            file=paste0(folder_salmon_lujambio_final,"all_samples_mouse_human_genes_read_count.tsv"),
            sep="\t",quote=FALSE,row.names = FALSE)

table_mouse<-table(expression.Read_count_human$Gene.stable.ID_mouse)
repetition_mouse_Gene<-table_mouse[table_mouse>1]
write.table(repetition_mouse_Gene,
            file=paste0(folder_salmon_lujambio_final,"repeat_mouseGeneName_vs_humanGeneName.tsv"),
            sep="\t",quote=FALSE,row.names = FALSE)

table_human<-table(expression.Read_count_human$Gene.stable.ID_human)
repetition_human_Gene<-table_human[table_human>1]
write.table(repetition_human_Gene,
            file=paste0(folder_salmon_lujambio_final,"repeat_humanGeneName_vs_mouseGeneName.tsv"),
            sep="\t",quote=FALSE,row.names = FALSE)

### remove line where sum for all samples == 0 per gene
sum_reads <- rowSums(expression.Read_count_human[6:76])
no_reads_row <- which(sum_reads ==0)
expression.Read_count_human_clean <- expression.Read_count_human[-no_reads_row,]
write.table(expression.Read_count_human_clean,
            file=paste0(folder_salmon_lujambio_final,"all_samples_mouse_human_genes_read_count_noZerorow.tsv"),
            sep="\t",quote=FALSE,row.names = FALSE)



###### Creation table with Mouse and Human genes for TPM
txi_scaledTPM <- tximport(files, type="salmon", tx2gene=tx2gene,countsFromAbundance="scaledTPM")

expression.scaledTPM <-as.data.frame(txi_scaledTPM$counts)

write.table(expression.scaledTPM,
            file=paste0(folder_salmon_lujambio_final,"all_samples_mouse_genes_scaledTPM.tsv"),
            sep="\t",quote=FALSE,row.names = FALSE)

#### give the gene as human name
musGenes_ensemb_scaledTPM <- rownames(expression.scaledTPM)

# Basic function to convert mouse to human gene names
genesV2_scaledTPM = getLDS(attributes = c("mgi_symbol","ensembl_gene_id"), filters = "ensembl_gene_id",
                 values = musGenes_ensemb_scaledTPM , mart = mouse, 
                 attributesL = c("hgnc_symbol","ensembl_gene_id","external_gene_name"), 
                 martL = human, uniqueRows=T)

colnames(genesV2_scaledTPM)[2]<-"Gene.stable.ID_mouse"
colnames(genesV2_scaledTPM)[4]<-"Gene.stable.ID_human"

expression.scaledTPM$Gene.stable.ID_mouse <- rownames(expression.scaledTPM)
dim(expression.scaledTPM)
dim(genesV2_scaledTPM)
expression.scaledTPM_human <- merge(genesV2_scaledTPM,expression.scaledTPM,
                                    by="Gene.stable.ID_mouse")
dim(expression.scaledTPM_human)

write.table(expression.scaledTPM_human,
            file=paste0(folder_salmon_lujambio_final,"all_samples_mouse_human_genes_scaledTPM.tsv"),
            sep="\t",quote=FALSE,row.names = FALSE)

### remove line where sum for all samples == 0 per gene
sum_reads_scaledTPM <- rowSums(expression.scaledTPM_human[6:76])
no_reads_row_scaledTPM <- which(sum_reads_scaledTPM == 0)
expression.scaledTPM_human_clean <- expression.scaledTPM_human[-no_reads_row,]
write.table(expression.scaledTPM_human_clean,
            file=paste0(folder_salmon_lujambio_final,"all_samples_mouse_human_genes_scaledTPM_noZerorow.tsv"),
            sep="\t",quote=FALSE,row.names = FALSE)

