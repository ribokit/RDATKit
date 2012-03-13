from matplotlib.pylab import *
import pdb
import pickle
import sys
import scipy.stats as stats
from rdatkit.settings import *

figure(2)
title('Combined')
db = pickle.load(open(sys.argv[1]))
params = pickle.load(open(sys.argv[1].replace('.db', '.dists')))
scale = 0.6
if 'dms' in sys.argv[1]:
    scale = 0.5
name = sys.argv[1].replace('.db','')
cl = zip(rand(13), rand(13), rand(13))
all_distobjects = {}
for i, k in enumerate(db):
    if len(db[k]) > 0:
        figure(2)
        if k not in ('all', 'unpaired'):
	    hist(db[k], 100, alpha=0.6, color=cl[i], label=k)
	figure(1)
	clf()
	title(k)
	if dists[k] == stats.gamma:
	    dist = dists[k](params[k][0], loc=params[k][1], scale=params[k][2])
	    #dist = dists[k](params[k][0], scale=params[k][2])
	if dists[k] == stats.expon or dists[k] == stats.norm:
	    dist = dists[k](params[k][0], params[k][1])
        if dists[k] == stats.beta:
	    dist = dists[k](params[k][0], params[k][1], loc=params[k][2], scale=params[k][3])
        all_distobjects[k] = dist
	print 'Plotting %s as %s' % (k, dists[k])
        #plot([dist.pdf(x) for x in frange(min(db[k]), max(db[k]), 0.01)])
        n, bins, patches = hist(db[k], 100, alpha=0.6, color=cl[i], label=k)
        #n, bins, patches = hist(db[k], 100, normed=1, alpha=0.6, color=cl[i])
        #plot(bins, dist.pdf(bins)/dist.pdf(bins).max(), c='r')
	savefig('%s_%s.png' % (name, k))
    else:
	print 'Skipping %s, no data found' % k
clf()
title('Paired vs unpaired')
npaired, bin, patches = hist(db['helices'], 100, normed=1, alpha=0.6, color='b', label='Paired')
nunpaired, bin, patches = hist(db['unpaired'], 100, normed=1, alpha=0.6, color='r', label='Unpaired')
ycap = min((100, max(npaired), max(nunpaired)))
ylim(0, ycap)
xlim(0,3.5)
legend()
savefig('%s_paired_vs_unpaired.png' % name)
clf()
title('Edge pairs vs internal pairs')
npaired, bin, patches = hist(db['internalpairs'], 100, alpha=0.6, color='b', label='Internal Pairs')
nunpaired, bin, patches = hist(db['edgepairs'], 100, alpha=0.6, color='r', label='Edge Pairs')
ycap = min((100, max(npaired), max(nunpaired)))
ylim(0, ycap)
xlim(0,2)
legend()
savefig('%s_edge_vs_internal.png' % name)

figure(2)
ylim(0,ycap)
xlim(0,2.3)
legend()
savefig('%s_combined.png' % name)
paireddist = lambda x: dists['helices'](params['helices'][0], loc=params['helices'][1], scale=params['helices'][2]).pdf(x)
unpaireddist = lambda x: dists['unpaired'](params['unpaired'][0], loc=params['unpaired'][1], scale=params['unpaired'][2]).pdf(x)
figure(1)
clf()
title('Unpaired/Paired likelihood ratio')
plot(arange(0,3,0.01), [scale*log(unpaireddist(x)/paireddist(x)) for x in arange(0,3,0.01)], color='b', label='log-likelihood ratio')
plot(arange(0,3,0.01), [log(x + 1)*2.6 - 0.8 for x in arange(0,3,0.01)], color='r', label='pseudoenergy')
legend()
savefig('%s_unpaired_vs_paired_ratio.png' % name)
figure(10)
clf()
title('Unpaired vs paired distributions')
plot(arange(0,3,0.01), [unpaireddist(x) for x in arange(0,3,0.01)], color='r', label='Unpaired distribution')
plot(arange(0,3,0.01), [paireddist(x) for x in arange(0,3,0.01)], color='b', label='Paired distribution')
legend()
savefig('%s_unpaired_vs_paired_dists.png' % name)
