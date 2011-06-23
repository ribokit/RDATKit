from matplotlib.pylab import *
import pdb
import pickle
import sys
import scipy.stats as stats
from rdatkit.settings import *

figure(2)
title('Combined')
db = pickle.load(open(sys.argv[1]))
params = pickle.load(open(sys.argv[1] + '.dists'))
name = sys.argv[1].replace('.db','')
cl = zip(rand(11), rand(11), rand(11))
for i, k in enumerate(db):
    if len(db[k]) > 0:
        figure(2)
        hist(db[k], 100, alpha=0.6, color=cl[i])
	figure(1)
	clf()
	title(k)
	if dists[k] == stats.gamma:
	    #dist = dists[k](params[k][0], loc=params[k][1], scale=params[k][2])
	    dist = dists[k](params[k][0], scale=params[k][2])
	if dists[k] == stats.expon or dists[k] == stats.norm:
	    dist = dists[k](params[k][0], params[k][1])
	print 'Plotting %s as %s' % (k, dists[k])
        #plot([dist.pdf(x) for x in frange(min(db[k]), max(db[k]), 0.01)])
        n, bins, patches = hist(db[k], 100, normed=1, alpha=0.6, color=cl[i])
        plot(bins, dist.pdf(bins), c='r')
	savefig('%s_%s.png' % (name, k))
    else:
	print 'Skipping %s, no data found' % k
figure(2)
savefig('%s_combined.png' % name)
