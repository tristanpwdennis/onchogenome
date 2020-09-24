
library(pegas)
library(knitr)
library(igraph)
library(RColorBrewer)
library(tidyverse)

###############################################
#Function relmatrixwrangle
#This function takes the output of ibsrelate that has been parsed with do_derived_stats
#It joins the output to the metadata, one side for each sample (in the pair)
#It then specifies type of comparison (by pasting worm type for each sample together)
#It thencreates columns specifying whether a comparison is within brood 
#and then outputs the finished df (with only mfmf comparisons) for plotting and analysis
relmatrixwrangle <- function(ibsrelateoutput, metadata) {
#add metadata to each 'side' of the paired data on key = sample number (bamlist order)
ibsrelateoutput <- left_join(ibsrelateoutput, metadata, by = c('ind1' = 'sample_num')) %>% 
  left_join(., metadata, by = c('ind2' = 'sample_num')) 
#Here, I define what each pair represents: a caculation between mf, or between mf and a male, or an mf and a female
ibsrelateoutput$typecomp <- paste0(ibsrelateoutput$worm_type.x, "_",ibsrelateoutput$worm_type.y)
#here, I create variables specifying whether members of a pair
ibsrelateoutput <- ibsrelateoutput %>% 
  mutate(degree_bro = ifelse(newsampleid.x == newsampleid.y, 'withinbrood', 'betweenbrood')) %>% 
  mutate(degree_pat = ifelse(patient_ID.x == patient_ID.y, 'within_patient', 'between_patient')) 
#filter out non-mf comparisons for now
ibsrelateoutput <- ibsrelateoutput %>% 
  filter(., typecomp == 'mf_mf')
#get the fraction of the total GL's involved in the comparison
ibsrelateoutput$fracsites <- ibsrelateoutput$nSites / 884673
return(ibsrelateoutput)
}


###########################
#function for plotting r1 and r0

plotr1r0 <- function(rel) {
  #plot R1/R0
  filter(rel, fracsites > 0.5) %>% 
    ggplot(., aes(x=R1, y=R0, colour = degree_bro)) +
    geom_point() +
    theme_minimal(base_size = 22)
}

###########################
#function for plotting r1 and rKING robust kinship
plotr1king <- function(rel) {
  #plot R1/King ROBUST KINSHIP
  filter(rel, fracsites > 0.5) %>% 
    ggplot(., aes(x=R1, y=Kin, colour = degree_bro)) +
    geom_point() +
    theme_minimal(base_size = 22)
}

###########################
#function for calculating the relationship based on specific values of R1, King and R0
definetherelationship <- function(rel) {
rel <- rel %>% mutate(.,
          relationship = case_when(
          R0 < 0.2 & R1 < 0.5 ~ 'h-sib',
          R1 > 0.5  ~ 'f-sib',
          TRUE ~ NA_character_)
                 )
return(rel)  
}

###########################
#function for plotting r1 and and R0 with inferred rel coloured in
plotinferredrel <- function(rel) {
  rel %>% filter(., fracsites > 0.8) %>% 
    ggplot(aes(x=R1, y = R0, colour = relationship)) +
    geom_point()
}


###########################
#function for creating initial igdfs
makeinitialigdfs <- function(rel, metadata) {
  edges <- rel %>% 
    filter(fracsites > 0.7) %>% 
    select(sample_id.x, sample_id.y, relationship) %>% 
    drop_na()
  #create node (igraph calls these vertices) attrribute dataframe (selecting metadata rows that are in the relatedness frame)
  nodes <- metadata %>% select(., sample_id, newsampleid) %>%
    filter(sample_id %in% edges$sample_id.x | sample_id %in% edges$sample_id.y)
    
  edgesnodes <- list(edges, nodes)
  
  return(edgesnodes)
  
  }

###########################
#function for creating colour palettes from node sample lists
make_ig_color_pals_from_samples <- function(nodeframe, palette) {
  #define group colours
  colgroups <- as.factor(nodeframe[,2])
  n<-nlevels(colgroups)
  pal <- brewer.pal(n,palette)
  #make factor (sample): colour mapping
  vertex.col <- pal[colgroups]
  return(list(pal, colgroups, vertex.col))
}

###########################
###function for creating haplotype network object using pegas
#I have not asked PEGAS to render it as the plotting options are
#abhorrent to me
makehaplotypenet <- function(d) {
  d <- read.dna(input, format = 'fasta')
  h <- pegas::haplotype(d)
  (hnet <- pegas::haploNet(h))
  return(hnet)
}

###########################
#function for extracting individual haplotype membership from haplotype object (pegas)
getindhap <- function(h){
  d <- read.dna(input, format = 'fasta')
  h <- pegas::haplotype(d)
  ind.hap<-with(
    stack(setNames(attr(h, "index"), rownames(h))),
    table(hap=ind, pop=rownames(d)[values]))
  return(ind.hap)
}

###########################
#function for get hamming distance from fasta alignment
getdistmatrix <- function(d){
  d <- read.dna(input, format = 'fasta')
  h <- pegas::haplotype(d)
  distmatrix <- dist.hamming(h)
  distmatrix <- as.matrix(distmatrix)
  return(distmatrix)
}     

###########################
##convert haplotype network to igraph, and filter out ungrouped haplotypes
mito_to_haplo_igraph <- function(metadata, dnafile, haplonet) {
  #get individual haplotypes
  ih <- getindhap(dnafile)  
  #get individual haplotype memberships from table generated above
  membship <- as.data.frame(ih) %>% filter(Freq == 1)
  #create our map that contains haplotype information and sample information - to get per-sample haplotype info
  map <- metadata %>% 
    mutate(., knownbrood = case_when(worm_type == 'mf' ~ newsampleid,TRUE ~ 'unknown')) %>% 
    select(., sample_id, newsampleid, knownbrood, worm_type) %>% 
    left_join(., membship, by = c('sample_id' = 'pop')) %>% 
    drop_na() %>% 
    distinct()
  #derive distance matrix from haplotype object
  distmatrix <- getdistmatrix(input)
  #longform distance matrix for pairs in BOTH DIRECTIONS (important you get both so that the igraph can be found)
  #turn matrix into longform pairs
  searchdf<-
    data.frame(col=colnames(distmatrix)[col(distmatrix)], row=rownames(distmatrix)[row(distmatrix)], dist=c(distmatrix)) %>% 
    #create an id for each pair (in both directions)
    mutate(., pair = paste0(col, "_",row)) 
  
  #convert haplotype network into igraph
  net.igraph <- pegas::as.igraph.haploNet(hnet,
                                          directed = FALSE, #Leave network arrows off as they are not relevant to our case.
                                          use.labels = TRUE,
                                          altlinks = FALSE) #only get main links between haplotypes or you'll end up with tons of extra lines
  #convert graph structure to data frame
  igdf <- igraph::as_data_frame(net.igraph, "both")
  #create pair ids from igraph df
  pairdf <- igdf$edges %>% mutate(., pair = paste0(from, "_",to)) 
  #subset dist df by columns IN DIRECTION OF igraph df to get igraph df + distances (which will be edge lengths)
  pairdf <- filter(searchdf, pair %in% pairdf$pair) %>% select(col, row, dist, pair)
  colnames(pairdf) <- c("from", "to", "dist", "pairid")
  #replace original edges with our edges plus attributes(including pairid for sanity checking)
  igdf$edges <- pairdf
  #create vertex df and add attributes from metadata
  igdf$vertices <- igdf$vertices %>% left_join(map, by = c("name" = "hap")) %>% select(name, knownbrood) %>% distinct()
  #Let's set a cutoff, based on the HaploNet plot, of 2 SNPs difference - to define our clustering
  igdf$edges <- igdf$edges %>% filter(dist < 3)
  #return the graph dataframe
  return(igdf)
}


makemitomap <- function(mitonet, mitogrups, haps) {
  ##########################
  #Create map of sample:haplotype for mitochondria
  #add mitochondrial group as as a vertex attribute
  V(mitonet)$mitogroup <- membership(mitogrups)
  mitoigdf <- igraph::as_data_frame(mitonet, "both")
  map <- as.data.frame(getindhap(h)) %>% filter(Freq > 0)
  mitomap<-mitoigdf$vertices %>%
    left_join(., map, by = c("name" = "hap")) %>% 
    select(pop, name, mitogroup) %>% 
    left_join(., metadata, by = c('pop' = "sample_id")) %>% 
    select(pop, name, mitogroup)
  colnames(mitomap) <- c("sample_id", "hap", "mitogroup")
  return(mitomap)
}

makesibmap <- function(sibnet, mitomap, fsibgrups) {
  ##########################
  #convert sib netrwortk into dataframe with community membership as a vertex attribute
  #add membership as vertex attribute and make dataframe of sample id and cluster membership
  V(sibnet)$community <- membership(fsibgrups)
  mem<- igraph::as_data_frame(sibnet, "both") 
  #create map of f-sib group membership
  #extract vertices and their community characteristics
  fsibmap <- mem$vertices %>% 
    left_join(mitomap, by = c("name" = "sample_id")) %>% 
    select(name, community, newsampleid)
  colnames(fsibmap) <- c("sample_id", "fsibgroup", "newsampleid")
  return(fsibmap)
}

