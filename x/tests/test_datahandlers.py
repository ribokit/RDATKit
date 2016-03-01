from rdatkit import datahandlers
import sys
import numpy as np
print 'Creating RDATFile object'
rdat = datahandlers.RDATFile()
print 'Loading RDAT file'
rdat.load(open(sys.argv[1]))
print 'Validating file'
rdat.validate()
print 'Converting to isatab'
isatab = rdat.toISATAB()

data = np.zeros([2,2])
data[0,:] = 1
data_annotations = [{'annotation':['hello'], 'annotation2':'hello2'}, {'annotation':['bla']}]
sequence ='AAAGGGAAAUUU'

rdat = datahandlers.RDATFile()
rdat.save_a_construct('test', data, sequence, '...(((...)))', 0, {'yes':'no'}, data_annotations, 'test.rdat', comments='commenting', version=0.32)
