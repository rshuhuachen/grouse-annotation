##### In this script we will separate the gff file into regions of interest based on the great tit ########

#### Packages #####
#remotes::install_github("jasongraf1/JGmisc")
pacman::p_load(BiocManager, rtracklayer, GenomicFeatures, BiocGenerics, JGmisc, data.table)

#### Genome data ####
gff <- makeTxDbFromGFF("data/PO2979_Lyrurus_tetrix_black_grouse.annotation.gff.gz", format="gff3", organism="Lyrurus tetrix") 

## divide up ###

promoters <- promoters(gff, upstream=2000, downstream=200, columns=c("tx_name", "gene_id")) # From NIOO
TSS <- promoters(gff, upstream=300, downstream=50, columns=c("tx_name", "gene_id")) # TSS as in Laine et al., 2016. Nature Communications
downstream <- flank(genes(gff), 10000, start=FALSE, both=FALSE, use.names=TRUE)
upstream <- promoters(genes(gff), upstream=10000, downstream=0)

exons_gene <- unlist(exonsBy(gff, "gene")) # group exons by genes

introns <- unlist(intronsByTranscript(gff, use.names=TRUE))

fiveUTRs <- unlist(fiveUTRsByTranscript(gff, use.names=TRUE))
threeUTRs <- unlist(threeUTRsByTranscript(gff, use.names=TRUE))

### write out files

export(genes(gff), "output/genes.gff3")
export(promoters, "output/promoters.gff3")
export(TSS, "output/TSS.gff3")
export(downstream, "output/downstream.gff3")
export(upstream, "output/upstream.gff3")

export(exons_gene, "output/exons_gene.gff3")

introns@ranges@NAMES[is.na(introns@ranges@NAMES)]<-"Unknown"
export(introns, "output/introns_transcripts.gff3")

fiveUTRs@ranges@NAMES[is.na(fiveUTRs@ranges@NAMES)]<-"Unknown"
export(fiveUTRs, "output/fiveUTRs.gff3")

threeUTRs@ranges@NAMES[is.na(threeUTRs@ranges@NAMES)]<-"Unknown"
export(threeUTRs, "output/threeUTRs.gff3")
