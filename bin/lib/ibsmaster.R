################
#load all our shit
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#mitochondrial genome fasta of all samples
input <- c('../data/unamplified_amplified_mf.fas')
#ibsrelate output
testpath0 = '../data/ibs.model0.results.ibspair'
#source file with ibs functions
source("ibsfunctions.R")
#source file that generates metadata
source("metadata.R")
###
source("ibslib.R")

###############################################################################################################################
#IBSMaster.R
#This script is for calculating the number of unique adult genotypes in a sample of O. volvulus microfilaria
#It takes: relatedness output from IBSRelate, and a .fasta file of mitochondrial genomes
#It uses iGraph to generate clusters of individuals by: full-sib group, and mitochodrial haplotype, and overlays them (crudely)
#over one another to derive the pattern of sibling and mtDNA connection
#Most of the work under the hood is done by the three libraries:
#ibsfunctions.R - taken from ANGSD website to calculate relatedness from IBSRelate output
#metadata.R - builds the metadata for the oncho project
#ibslib.R - assorted, rather project-specific functions that I've stowed here for readibility
#This code is meant to be rendered in the corresponding RMarkdown, but is presented here for those who want to take a look
#and see (and check!) what I have done to make my calculation
###############################################################################################################################

#################
#Build metadata (getseqandsamp is from metadata.R)
metadata <-getseqandsamp()

#################
#Calculate kinship estimates from IBSRelate output
ibsrelateoutput <-  do_derived_stats(read_ibspair_model0(testpath0))

#################
#Do some wrangling on the relatedness data:
#for each pairwise comparison, calculate whether it is "within" or "between" brood, patient, etc
rel <- relmatrixwrangle(ibsrelateoutput, metadata)
rel <- rel %>% filter(typecomp == 'mf_mf')


#################
#There's a lot of noisy samples here, so I've removed comparisons 
#that feature <0.7 of the total GL callset
#code below is for tinkering to see (fracsites is the fraction of the total callset involved in the
#estimation): 0-1

filter(rel, fracsites > 0.1) %>% 
  ggplot(., aes(x=R1, y=R0, colour = degree_bro)) +
  geom_point() +
  theme_minimal(base_size = 22)

#the following functions apply a filter of 0.5 (but mostly the same results up to 0.8)
#################
#Plot R1/R0
r1r0 <- plotr1r0(rel, 0.5)

#################
#Plot R1/R0
r1king <- plotr1king(rel, 0.5)

#################
#Based on these plots, define thresholds for different values of relatedness
#this groups individuals based on values of R0, R1 and KING
rel <- definetherelationship(rel)

######################
#plot relationships with points coloured by inference

plotinferredrel(rel)

######################
#create r0/r1 and r1/king plot

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  r1r0 + theme(legend.box.margin = margin(0, 0, 0, 12)) +
    theme(legend.position = "bottom")
)

prow <- plot_grid(
  r1r0 + theme(legend.position="none"),
  r1king + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)

plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))

#################
#create the dataframe that will be used to plot the sib network in iGraph
#outputs a list of dataframes (one for the edges, one for the vertices)
out <- makeinitialigdfs(rel, metadata)

#################
#create the iGraph object, create a color palette from the vertex dataframe, color the edges by relationship, and plot
#the curly braces on the plot function enable the legend and plot to be plotted together in .Rmd
sibnet <- graph_from_data_frame(out[[1]], out[[2]], directed = FALSE)

palstuff <- make_ig_color_pals_from_samples(out[[2]], "Set1")
#colour edges by inferred relationship
E(sibnet)$color <- as.factor(E(sibnet)$relationship)
#plot network
{plot(sibnet, vertex.label.cex = .7,  vertex.size=30, vertex.color=palstuff[[3]], vertex.label=NA) 
  legend("topleft",bty = "n",
         legend=levels(palstuff[[2]]),
         fill=palstuff[[1]], border=NA)}

#################
#Filter out half-sibs leaving only fullsibs
#because the halfsib connections are a bit confusing
out[[1]] <- out[[1]] %>% filter(., relationship == 'f-sib')

#################
#rebuild iGraph object from filtered df (palette remains the same but edges are now uncoloured)
sibnet <- graph_from_data_frame(out[[1]], out[[2]], directed = FALSE)

########plot tidier network -  palette is retained
{plot(sibnet, vertex.label.cex = .7,  vertex.size=30, vertex.color=palstuff[[3]],vertex.label=NA) 
  legend("topleft",bty = "n",
         legend=levels(palstuff[[2]]),
         fill=palstuff[[1]], border=NA)}

#################
#build haplotype network from fasta input
hnet <- makehaplotypenet(input)

#you can plot the hnet here but I haven't done much to prettify it
#I struggled with PEGAS' plotting, gave up, plotted a network in PopART
#that is rendered in the .Rmd - this one is just for extracting information from
plot(hnet)
################
#convert haplotype network into iGraph object. then to a dataframe
#whilst it is a dataframe, pairwise distance between haplotypes are mapped as edge attrributes
#a mapping is created that specifies which samples are members of the respective haplotypes too
#we end up with a list of dataframes, vertices: that contains the haplotype name, and the brood/s from which that haplotype was derived
#and edges, which contains the filtered edges: here I define a 'cluster' of a brood of samples that are  <2 different from one another
#perhaps I will revise this dependent on other intrauterine data
ig <- mito_to_haplo_igraph(metadata, input, hnet)

#create the mitochondrial igraph ibject
mitonet <- graph_from_data_frame(ig$edges, ig$vertices, directed = FALSE)
#generate vertext colours for the mitochondrial igraph
mitopal <- make_ig_color_pals_from_samples(ig$vertices, "Set1")

########plot mitochondrial network -  palette is retained
{plot(mitonet, vertex.label.cex = .7,  vertex.size=30, vertex.color=mitopal[[3]]) 
  legend("topleft",bty = "n",
         legend=levels(mitopal[[2]]),
         fill=mitopal[[1]], border=NA)}

#simple clustering of sibling network based on FS connectivity
fsibgrups<-fastgreedy.community(sibnet)
mitogrups<-fastgreedy.community(mitonet)

#plot mitochondrial clustering
plot(mitogrups, mitonet)
#plot fullsib clustering
plot(fsibgrups, sibnet)

#########################
#get individual haplotypes
h <- getindhap(dnafile = input)

##########################
#make mito map
mitomap<-makemitomap(mitonet, mitogrups, h)

##########################
#make mapping of sample to fullsib and mito group
fsibmap<-makesibmap(sibnet, mitomap, fsibgrups)




####now we totally mangle the relatedness dataframe
####this is a bit of a beast
#I join the mitochondrial group data (from mitomap) to each side of the pairs, and generate values for whether each pair belongs to the same
#mitochondrial group or not, drop a load of na values from parts of rel that are not in the final network
#then we join the full-sibling cluster information (fsibmap), selecting pairs of the DIFFERENT F-SIB GROUP THAT ARE PART OF THE SAME MITOCHONDRIAL GROUP (see the filters)
#the result is: a dataframe that comprises the adges (connecitons between members of different f-sib groups that are linked by the same mitochondrial group)
mitosibedge <- rel %>% 
  left_join(., mitomap, by = c("sample_id.x" = "sample_id")) %>% 
  left_join(., mitomap, by = c("sample_id.y" = "sample_id")) %>% 
  mutate(mitosib = ifelse(mitogroup.x == mitogroup.y, 'mitosib', 'mitononsib')) %>% 
  select(sample_id.x, sample_id.y, relationship, mitosib) %>% 
  drop_na() %>% 
  left_join(., fsibmap, by = c("sample_id.x" = "sample_id")) %>% 
  left_join(., fsibmap, by = c('sample_id.y' = 'sample_id')) %>% 
  select(., fsibgroup.x, fsibgroup.y, mitosib) %>% 
  filter(., mitosib == 'mitosib') %>% 
  filter(fsibgroup.x != fsibgroup.y)

#create vertices from membership dataframe
mitosibvertex <- fsibmap %>% 
  select(fsibgroup, newsampleid) %>% distinct()

#create graph object, fill by palette we created above, with corresponding legend
msibgraph  <- graph_from_data_frame(mitosibedge, mitosibvertex, directed = F)

mitopal <- make_ig_color_pals_from_samples(mitosibvertex, "Set1")

{plot(msibgraph, vertex.label.cex = .7,  vertex.size=30, vertex.color=mitopal[[3]]) 
  legend("topleft",bty = "n",
         legend=levels(mitopal[[2]]),
         fill=mitopal[[1]], border=NA)}

#convert this back to a dataframe for counting GTs
mitodf <- igraph::as_data_frame(msibgraph, "both")

#calculate unique adult genotypes-
#This is the number of vertices (unique adult genotypes)*2, minus the number of mitochondrial links
count(mitodf$vertices)*2 - count(mitodf$edges)

#Number of unique female genotypes - one for every f-sib cluster minus one for every shared mitochondrial haplotype
count(mitodf$vertices) - count(mitodf$edges)

metadata %>% 
  filter(short_plate_name == 'WS_08_19_MF') %>% 
  mutate(coverage = (number.of.mapped.reads*75) / 97000000) %>% 
  ggplot(aes(y=coverage)) +
  theme_half_open(12) + 
  geom_boxplot() +
  geom_jitter()


