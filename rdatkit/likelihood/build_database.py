#!/usr/bin/python

import sys
import os
import pickle
from numpy import *
import argparse
from rdatkit.datahandlers import RDATFile
from rdatkit.secondary_structure import SecondaryStructure
from helpers import normalize
from scipy.stats.stats import scoreatpercentile
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('rdatdir', type=str)
parser.add_argument('outfile', type=str)
parser.add_argument('--normalize', dest='normalize', const=True, default=False, action='store_const')
parser.add_argument('--nooutliers', dest='nooutliers', const=True, default=False, action='store_const')

args = parser.parse_args()

fragtypes = ['all', 'helices', 'interiorloops', 'hairpins', 'dangles', 'bulges',\
	     '2wayjunctions', '3wayjunctions', '4wayjunctions', '5wayjunctions', 'unpaired', 'edgepairs', 'internalpairs']
db = {}
dbidx = {}
for t in fragtypes:
    db[t] = []
    dbidx[t] = {}
for filename in os.listdir(args.rdatdir):
    if not os.path.isdir(args.rdatdir+'/'+filename):
	print filename
	rdat = RDATFile()
	rdat.load(open(args.rdatdir+'/'+filename))
	for cname in rdat.constructs:
	    construct = rdat.constructs[cname]
	    struct = SecondaryStructure(construct.structure)
	    frags = struct.explode()
	    for data in construct.data:
		if (('mutation' not in data.annotations) or \
		   ('mutation' in data.annotations and \
		 'WT' in data.annotations['mutation'])):
                    if 'modifier' in data.annotations:
			if args.normalize:
			    normvals = normalize(data.values)
			else:
			    normvals = data.values
			iqr = scoreatpercentile(normvals, 75) - scoreatpercentile(normvals, 25)
			for fragtype in frags:
                            db['all'].extend(normvals)
			    dbidx['all'] = dict([((construct.name, construct.seqpos[i]), v) for i, v in enumerate(normvals)])
			    fraglist = frags[fragtype]
			    for frag in fraglist:
				vals = []
				pos = []
				for idx in frag:
				    try:
				        iddx = construct.seqpos.index(idx + construct.offset + 1)
					if ('DMS' in data.annotations['modifier'] and construct.sequence[idx].upper() not in ['A', 'C']) or\
					   ('CMCT' in data.annotations['modifier'] and construct.sequence[idx].upper() not in ['G', 'U']) or\
					     (args.nooutliers and (normvals[iddx] < 0)):
					     #(args.nooutliers and (normvals[iddx] > 1.5*iqr or normvals[iddx] < 0)):
					    continue
					if construct.structure[idx] == '.':
					    db['unpaired'].append(normvals[iddx])
					    dbidx['unpaired'][(construct.name, idx + construct.offset + 1)] = normvals[iddx]
					if construct.structure[idx] in (')', '('):
					    db['helices'].append(normvals[iddx])
					    dbidx['helices'][(construct.name, idx + construct.offset + 1)] = normvals[iddx]
					    if '.' in (construct.structure[idx-1], construct.structure[idx+1]):
						db['edgepairs'].append(normvals[iddx])
						dbidx['edgepairs'][(construct.name, idx + construct.offset + 1)] = normvals[iddx]
					    else:
						db['internalpairs'].append(normvals[iddx])
						dbidx['internalpairs'][(construct.name, idx + construct.offset + 1)] = normvals[iddx]
					val = normvals[iddx]
					if not isnan(val):
					    vals.append(val)
					    pos.append(idx + construct.offset + 1)
				    except ValueError:
					pass
				if len(vals) > 0 and fragtype != 'helices':
				    db[fragtype].extend(vals)
				    for i, v in enumerate(vals):
					dbidx[fragtype][(construct.name, pos[i])] = v

pickle.dump(db, open(args.outfile,'w'))
pickle.dump(dbidx, open(args.outfile + '.idx', 'w' ))
