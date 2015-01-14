import scipy.stats as stats

EXTERNAL_DIR = '/MATLAB_Code'
RNA_STRUCTURE = EXTERNAL_DIR + 'RNAstructure/'
RNA_STRUCTURE_FOLD = EXTERNAL_DIR + 'RNAstructure/exe/Fold'
RNA_STRUCTURE_ALLSUB = EXTERNAL_DIR + 'RNAstructure/exe/AllSub'
RNA_STRUCTURE_STOCHASTIC = EXTERNAL_DIR + 'RNAstructure/exe/stochastic'
RNA_STRUCTURE_ENERGY = EXTERNAL_DIR + 'RNAstructure/exe/efn2'
RNA_STRUCTURE_PARTITION = EXTERNAL_DIR + 'RNAstructure/exe/partition'
RNA_STRUCTURE_CT2DOT = EXTERNAL_DIR + 'RNAstructure/exe/ct2dot'
VIENNA_RNA = EXTERNAL_DIR + 'Biers/ViennaRNA-1.8.4/'
VARNA = '~/Downloads/Lab/VARNAv3-9.jar'
MAPPING_DATABASE_PATH = ''
dists = {}
for k in ['dangles', 'bulges', 'hairpins']:
    dists[k] = stats.gamma
for k in ['all', 'interiorloops', 'helices', '3wayjunctions', '4wayjunctions',\
          '5wayjunctions', 'unpaired', 'internalpairs', 'edgepairs']:
    dists[k] = stats.gamma
