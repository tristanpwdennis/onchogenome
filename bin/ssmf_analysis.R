library(data.table)
library(pegas)
library(ape)
library(RColorBrewer)

ssmf_fas <- read.dna("~/Library/Mobile Documents/com~apple~CloudDocs/Projects/onchocerca/onchogenome/data/oncho_ssmf_mito.fasta", format="fasta")
ssmf_meta <- fread("~/Library/Mobile Documents/com~apple~CloudDocs/Projects/onchocerca/onchogenome/metadata/seq_output/oncho_ssmf_mitodata.csv")
ssmf_res <- fread("~/Library/Mobile Documents/com~apple~CloudDocs/Projects/onchocerca/onchogenome/metadata/seq_output/oncho_ssmf_mitodata.csv")


#match up metadata to fasta names
netmetadata = ssmf_meta[match(rownames(ssmf_fas), ssmf_meta$mitoid), ]

#change fasta names to form id
names(ssmf_fas) <- netmetadata$loc

#make hap object from fasta
onchohaps <- pegas::haplotype(ssmf_fas)

#make table of now many individuals from each form belong to each cluster
ind.hap<-with(
  stack(setNames(attr(onchohaps, "index"), rownames(onchohaps))),
  table(hap=ind, individuals=names(ssmf_fas)[values]))

pal1 = RColorBrewer::brewer.pal("Dark2", n=6)

onchonet <- pegas::haploNet(onchohaps)
setHaploNetOptions(pie.colors.function=c("#66A61E","#E6AB02"))
plot(onchonet, 
     size = sqrt(attr(onchonet, 'freq')), 
     fast = FALSE, 
     show.mutation = 2, 
     labels= FALSE, 
     pie=ind.hap,
     cex = 0.3,
     legend = c(-40, 15))
legend(-20,15, fill=c("#66A61E","#E6AB02"), bty='n', cex=1, legend = c("Left", "Right"))


?haploNet
