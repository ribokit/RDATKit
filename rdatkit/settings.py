import scipy.stats as stats

EXTERNAL_DIR = '/home/tsuname/'
RNA_STRUCTURE = EXTERNAL_DIR + 'RNAstructure/5.5/'
RNA_STRUCTURE_FOLD = RNA_STRUCTURE + 'exe/Fold'
RNA_STRUCTURE_ALLSUB = RNA_STRUCTURE + 'exe/AllSub'
RNA_STRUCTURE_STOCHASTIC = RNA_STRUCTURE + 'exe/stochastic'
RNA_STRUCTURE_ENERGY = RNA_STRUCTURE + 'exe/efn2'
RNA_STRUCTURE_PARTITION = RNA_STRUCTURE + 'exe/partition'
RNA_STRUCTURE_CT2DOT = RNA_STRUCTURE + 'exe/ct2dot'
VIENNA_RNA = EXTERNAL_DIR + 'ViennaRNA-1.4/'
VIENNA_RNA_ENERGY = '/usr/bin/RNAeval'
VIENNA_RNA_FOLD = '/usr/bin/RNAfold'
VIENNA_RNA_SUBOPT = '/usr/bin/RNAsubopt'
VARNA = '/home/tsuname/VARNA.jar'
MAPPING_DATABASE_PATH = '/home/tsuname/Documents/lab/rhiju/rdatkit/likelihood/databases'
dists = {}
for k in ['dangles', 'bulges', 'hairpins']:
    dists[k] = stats.gamma
for k in ['all', 'interiorloops', 'helices', '3wayjunctions', '4wayjunctions',\
          '5wayjunctions', 'unpaired', 'internalpairs', 'edgepairs']:
    dists[k] = stats.gamma
