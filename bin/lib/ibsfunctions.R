# These functions read text files - the output of IBD applied to a group of at least two individuals.
# The *.ibspair has one line per pair of individuals
read_ibspair_model0 <- function(filepath){
  df = read.csv(filepath, sep ='\t')
  df['model'] = 'model0'
  df= within(df, pair <- paste(ind1, ind2, sep = '_'))
  
  
  df['fracA'] = df['pAA_AA'] + df['pCC_CC'] +df['pGG_GG'] + df['pTT_TT']
  
  df['fracB'] = (
    df['pAC_AA'] + df['pAG_AA'] + df['pAT_AA'] 
    + df['pAC_CC'] + df['pCG_CC'] + df['pCT_CC'] 
    + df['pAG_GG'] + df['pCG_GG'] + df['pGT_GG']
    + df['pAT_TT'] + df['pCT_TT'] + df['pGT_TT'])
  
  df['fracC'] = (
    df['pAA_CC'] + df['pAA_GG'] + df['pAA_TT'] 
    + df['pCC_AA'] + df['pCC_GG'] + df['pCC_TT'] 
    + df['pGG_AA'] + df['pGG_CC'] + df['pGG_TT']
    + df['pTT_AA'] + df['pTT_CC'] + df['pTT_GG'])
  
  df['fracD'] = (
    df['pAA_AC'] + df['pAA_AG'] + df['pAA_AT'] 
    + df['pCC_AC'] + df['pCC_CG'] + df['pCC_CT'] 
    + df['pGG_AG'] + df['pGG_CG'] + df['pGG_GT']
    + df['pTT_AT'] + df['pTT_CT'] + df['pTT_GT'])
  
  df['fracE'] = (
    df['pAC_AC'] + df['pAG_AG'] + df['pAT_AT'] 
    + df['pCG_CG'] + df['pCT_CT'] 
    + df['pGT_GT'])
  
  df['fracF'] = 0 # same as D
  df['fracG'] = 0 # same as C
  df['fracH'] = 0 # same as B
  df['fracI'] = 0 # same as A
  
  df['check'] = df['fracA'] + df['fracB'] + df['fracC'] +df['fracD'] +df['fracE'] # should be close to 1
  return(df)
}


read_ibspair_model1 <- function(filepath){
  df = read.csv(filepath, sep ='\t')
  df['model'] = 'model1'
  df = within(df, pair <- paste(ind1, ind2, sep = '_'))
  
  df['fracA'] = df['pAA_AA'] * 0.01
  df['fracB'] = df['pAB_AA'] * 0.01
  df['fracC'] = df['pAA_BB'] * 0.01
  df['fracD'] = df['pAA_AB'] * 0.01
  df['fracE'] = df['pAB_AB'] * 0.01
  
  df['fracF'] = 0 # same as D
  df['fracG'] = 0 # same as C
  df['fracH'] = 0 # same as B
  df['fracI'] = 0 # same as A
  
  df['check'] = df['fracA'] + df['fracB'] + df['fracC'] +df['fracD'] +df['fracE'] # should be close to 1
  
  return(df)
} 


read_ibspair_model2 <- function(filepath){
  df = read.csv(filepath, sep ='\t')
  df['model'] = 'model2'
  df = within(df, pair <- paste(ind1, ind2, sep = '_'))
  
  df['fracA'] = df['pAA_AA'] * 0.01
  df['fracB'] = df['pAB_AA'] * 0.01
  df['fracC'] = df['pAA_BB'] * 0.01
  df['fracD'] = df['pAA_AB'] * 0.01
  df['fracE'] = df['pAB_AB'] * 0.01
  # ignore estimates with more than two alleles
  
  df['fracF'] = 0 # same as D
  df['fracG'] = 0 # same as C
  df['fracH'] = 0 # same as B
  df['fracI'] = 0 # same as A
  
  df['check'] = df['fracA'] + df['fracB'] + df['fracC'] +df['fracD'] +df['fracE'] # should be close to 1
  
  return(df)
} 

do_derived_stats <- function(df){
  df['A'] = df['fracA'] * df['nSites']
  df['B'] = df['fracB'] * df['nSites']
  df['C'] = df['fracC'] * df['nSites']
  df['D'] = df['fracD'] * df['nSites']
  df['E'] = df['fracE'] * df['nSites']
  df['F'] = df['fracF'] * df['nSites']
  df['G'] = df['fracG'] * df['nSites']
  df['H'] = df['fracH'] * df['nSites']
  df['I'] = df['fracI'] * df['nSites']
  
  # Basic summary stats 
  df['HETHET'] = df['E']
  df['IBS0'] = df['C'] + df['G']
  df['IBS1'] = df['B'] + df['D'] + df['F'] + df['H']
  df['IBS2'] = df['A'] + df['E'] + df['I']
  df['fracIBS0'] = df['IBS0'] / df['nSites']
  df['fracIBS1'] = df['IBS1'] / df['nSites']
  df['fracIBS2'] = df['IBS2'] / df['nSites']
  df['fracHETHET'] = df['E'] / df['nSites']
  
  # the derived stats 
  df['R0'] = df['IBS0'] / df['HETHET']
  df['R1'] = df['HETHET'] / (df['IBS0'] +  df['IBS1'])
  # KING-robust kinship
  df['Kin'] = (df['HETHET'] - 2*(df['IBS0'])) / (df['IBS1'] + 2*df['HETHET'])
  df['Fst'] = (2*df['IBS0'] - df['HETHET']) / (2*df['IBS0'] + df['IBS1'] + df['HETHET'])
  
  # heterozygosity of each individual
  df['het_ind1'] = df['fracB'] + df['fracE']
  df['het_ind2'] = df['fracD'] + df['fracE']
  
  return(df)
}
