#script for generating comperehensive metadata from sequencing plates
getseqandsamp <- function() {
#get og sample data
  sampdata <- read.csv('/Users/tristanpwdennis/Software/onchotoolbox/data/og_uhb_dec_2018.csv')
  #get table specifying bam order
  bamorder <- read.csv('/Users/tristanpwdennis/Software/onchotoolbox/data/bamlist.csv')
  bamorder <- bamorder %>% select(sample_id, sample_num)
  
  #concatenate sequencing plates
  concat <- rbind(read.csv('/Users/tristanpwdennis/Projects/DNDI/Metadata-dev/sequencing_plates/WS_08_19_MALES.csv', stringsAsFactors = FALSE), 
                  read.csv('/Users/tristanpwdennis/Projects/DNDI/Metadata-dev/sequencing_plates/WS_08_19_MF.csv', stringsAsFactors = FALSE), 
                  read.csv('/Users/tristanpwdennis/Projects/DNDI/Metadata-dev/sequencing_plates/WS_11_19_MF.csv', stringsAsFactors = FALSE))
  #bind sample data to sequencing plates
  metadata <- sampdata %>% select(newsampleid, patient_ID, nodule_id) %>% left_join(., concat, by = c('newsampleid' = 'og_sample')) %>% filter(sequenced == TRUE)
  
  #join with bam ordering scheme
  metadata <-metadata %>% left_join(bamorder)
  
#return final metadata dataframe
return(metadata)

}

getseq <- function() {
  #get table specifying bam order
  bamorder <- read.csv('/Users/tristanpwdennis/Software/onchotoolbox/data/bamlist.csv')
  #concatenate sequencing plates
  concat <- rbind(read.csv('/Users/tristanpwdennis/Projects/DNDI/Metadata-dev/sequencing_plates/WS_08_19_MALES.csv', stringsAsFactors = FALSE), 
                  read.csv('/Users/tristanpwdennis/Projects/DNDI/Metadata-dev/sequencing_plates/WS_08_19_MF.csv', stringsAsFactors = FALSE), 
                  read.csv('/Users/tristanpwdennis/Projects/DNDI/Metadata-dev/sequencing_plates/WS_11_19_MF.csv', stringsAsFactors = FALSE))
  #bind bam order to metadata
  metadata <- bamorder %>% select(sample_id, sample_num) %>% left_join(., metadata)
  #return final metadata dataframe
  return(metadata)
  
}

getseqandsamp()









