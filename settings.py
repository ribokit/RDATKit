import scipy.stats as stats

EXTERNAL_DIR = '/home/tsuname/Documents/lab/rhiju/rdat/external/'
RNA_STRUCTURE = EXTERNAL_DIR + 'RNAstructure/'
RNA_STRUCTURE_FOLD = EXTERNAL_DIR + 'RNAstructure/exe/Fold'
RNA_STRUCTURE_ALLSUB = EXTERNAL_DIR + 'RNAstructure/exe/AllSub'
RNA_STRUCTURE_STOCHASTIC = EXTERNAL_DIR + 'RNAstructure/exe/stochastic'
RNA_STRUCTURE_PARTITION = EXTERNAL_DIR + 'RNAstructure/exe/Partition'
RNA_STRUCTURE_CT2DOT = EXTERNAL_DIR + 'RNAstructure/exe/ct2dot'
VIENNA_RNA = EXTERNAL_DIR + 'ViennaRNA-1.4/'
MAPPING_DATABASE_PATH = '/home/tsuname/Documents/lab/rhiju/rdatkit/likelihood/databases'
dists = {}
for k in ['dangles', 'bulges', 'hairpins']:
    dists[k] = stats.gamma
for k in ['interiorloops', 'helices', '3wayjunctions', '4wayjunctions',\
          '5wayjunctions']:
    dists[k] = stats.gamma

