library(pegas)
library(igraph)
library(tidyverse)

allmetadata <- read.csv('/Users/tristanpwdennis/Projects/DNDI/Metadata-dev/foranalysis/metadata-all.csv')
########################
#Step 1. Read in data, faff about, plot initial haplonet
#mitochondrial alignment for all samples
input <- c('~/Projects/DNDI/Data/mitochondrial/unamplified_amplified_mf_and_males/fasta/unamplified_amplified_mf_and_males.fas')
d <- read.dna(input, format = 'fasta')
mf<-read.csv('/Users/tristanpwdennis/Software/onchotoolbox/data/WS_08_19_MF_UDD.csv', stringsAsFactors = FALSE)

d <- read.dna(input, format = 'fasta')
h <- pegas::haplotype(d)
(net <- haploNet(h))
plot(net)

#get the sequence names
sn <- as.data.frame(attributes(d)$dimnames[[1]]) %>% 
    left_join(.,
            metadata,
            by = c('attributes(d)$dimnames[[1]]' = 'sample_id')
            ) %>% 
    select(., newsampleid)

#create table ind.hap that describes the individual membership of each haplotype group
ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)

#plot haplonet - individual samples are the components of the pie charts (from the ind.hap table)
plot(net, size=attr(net, "freq"), scale.ratio=0.2, pie=ind.hap) 

########################
#The plot is nice but not abundantly informative, and these plots are hard to manipulate
#The distances in the sequence are based on the hamming distance - a raw measure of SNP distances.
#I will extract the hamming distance matrix from the DNAbin object and use it to build and igraph and define clusters


mitodist <- dist.hamming(d)


























########################
#Step 2. Assign haplotypes to your site dataframe.
#This code is adapted from a blog post by Jimmy O'Donnell (jodonnellbio at gmail) with help from Kim Tenggardjaja.
#Original posting: http://jimmyodonnell.wordpress.com/2013/11/23/haplotype-networks-in-r-updated/
#Adaptations to this code by Claire Curry (curryclairem@gmail.com)

#Now, fill in the haplotypes for each site.
#You will have to do this with your real dataset.
#The site and haplotype dataframes must be in the SAME SAMPLE ORDER.

#Add empty column to sites dataframe for haplotypes
sites<-cbind(sn,
             haplotype=rep(NA,
                           nrow(sn)))
#Assign haplotypes to the samples automatically using a loop
for (i in 1:length(labels(h))){
  sites$haplotype[attr(h,
                       "index")[[i]]]<-i
}
########################
#Step 3. Convert to igraph format from haplonet.
#Here's my new adaptation to plot colors automatically in igraph, as suggested on the original blog post.
#Get igraph network
net.igraph <- pegas::as.igraph.haploNet(net,
                                        directed = FALSE, #Leave network arrows off as they are not relevant to our case.
                                        use.labels = TRUE,
                                        altlinks = FALSE) #only get main links between haplotypes or you'll end up with tons of extra lines

#Get the haplotype numbers in their weird order so you can order colors and sizes correctly.
newnums <- as.roman(as.character(V(net.igraph)$name))
V(net.igraph)$name <- newnums

plot(net.igraph)


########################
#Step 4.Count up haplotype frequencies and site frequencies.
#This code is adapted from a blog post by Jimmy O'Donnell (jodonnellbio at gmail) with help from Kim Tenggardjaja.
#Original posting: http://jimmyodonnell.wordpress.com/2013/11/23/haplotype-networks-in-r-updated/
#Adaptations to this code by Claire Curry (curryclairem@gmail.com)
#Count how many individuals for each haplotype
count.hap0<-ddply(sites,
                  .(haplotype),
                  summarise,
                  freq=length(haplotype))

#Be sure to add in our new ordering!
count.hap <- count.hap0[order(match(rownames(count.hap0),
                                    newnums)),"freq"]


#for each plot, build frequency matrix that tells how to fill the pies by site frequency.
dfcount_site<-ddply(sn,
                    .(sn, haplotype),
                    summarise,
                    freq=length(sn))
dc.site00 <- cast(haplotype ~ sn,
                  data = dfcount_site, 
                  value = "freq", 
                  fill = 0)
dc.site0 <-dc.site00
########################
#Step 5. Convert the site frequencies into the correct list format and the correct order to match the igraph object.
#Now, onto more igraph adaptations.
dc.site0$haplotype<-NULL #Remove the haplotype row before turning this into a matrix.
#Turn the haplotype frequencies into a matrix in the NEW IGRAPH ORDER.
dc.site1<-t(as.matrix(dc.site0[order(match(rownames(dc.site0),
                                           newnums)),]))
#You must do this step to match the ordering of the cast dataframe haplotypes to the ordering used in igraph,
#otherwise the colors and sizes will be off.

#Then split the matrix out into a list.
dc.site <- split(dc.site1, 
                 rep(1:ncol(dc.site1), 
                     each = nrow(dc.site1)))


#Remember, if you go back to your haplotype counts to check that colors and sizes are correct,
#that we rearranged dc.site to be in the order that the haplotypes are stored
#You can check V(net.igraph)$name to see them again.
#dc.site0's haplotype column shows the haplotype numbers, but doing order(match() then puts dc.site1 and dc.site into
#the order the haplotypes are stored in the igraph network.
#Hence, we use the following step to make it more obvious by naming each list element with the haplotype number.
names(dc.site) <- V(net.igraph)$name

########################
#Step 6. Plot the network.
plot(net.igraph, #the igraph object is always a list of 10 so far as I can tell, regardless of how many haplotypes you used.
     vertex.shape="pie",
     vertex.pie = dc.site,
     #vertex.pie specifies where you are getting the pie slice counts.
     vertex.pie.color=list(c("black",
                             "orange",
                             "blue")),
     #The colors are the number of categories you have, must be in a list, in the same order as your sites.
     vertex.size = 5*(count.hap), 
     #number of samples for each haplotype, multiplied to make bigger if desired.
     #vertex.label = NA)#, #for no vertex labels (ie no hapotype numbers added)
     #vertex.label.cex = 1)#, #if you want bigger haplotype numbers.
     vertex.label.dist = 2.5) #how far from pie you want the haplotype numbers)

#Add a legend.
legend(x= "bottomleft",
       bty="n",
       pt.bg=c("black",
               "orange",
               "blue"),
       col="black",
       pch=22,
       legend=c("Site A",
                "Site B",
                "Site C"))     