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



filter(rel, fracsites > 0.1) %>% 
  ggplot(., aes(x=R1, y=R0, colour = degree_bro)) +
  geom_point() +
  theme_minimal(base_size = 22)


rel %>% select(sample_id.x, sample_id.y, parent_id.x, parent_id.y, degree_bro) 


rel %>% filter(sample_id.x == sample_id.y)

