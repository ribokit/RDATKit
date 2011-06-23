import sys
import pdb
import scipy.stats as stats
import pickle
from rdatkit.settings import *
from matplotlib.pylab import *

db = pickle.load(open(sys.argv[1]))
params = {}
for k in db:
    if len(db[k]) > 0:
	vals = array(db[k])
	p = dists[k].fit(vals[vals >= 0])
	params[k] = p
	print 'For %s' % k
	print params[k]

pickle.dump(params, open(sys.argv[1]+'.dists','w'))
