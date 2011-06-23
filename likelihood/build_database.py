#!/usr/bin/python

import sys
import os
import pickle
from numpy import *
from datahandlers import RDATFile
from rdatkit.secondary_structure import SecondaryStructure
from helpers import normalize

fragtypes = ['helices', 'interiorloops', 'hairpins', 'dangles', 'bulges',\
	     '2wayjunctions', '3wayjunctions', '4wayjunctions', '5wayjunctions']
rdatdir = sys.argv[1]
if len(sys.argv) > 2 and sys.argv[2] == '--normalize':
    do_normalize =True
else:
    do_normalize = False
db = {}
for t in fragtypes:
    db[t] = []
for filename in os.listdir(rdatdir):
    if not os.path.isdir(rdatdir+'/'+filename):
	print filename
	rdat = RDATFile()
	rdat.load(open(rdatdir+'/'+filename))
	for cname in rdat.constructs:
	    construct = rdat.constructs[cname]
	    struct = SecondaryStructure(construct.structure)
	    frags = struct.explode()
	    for data in construct.data:
		if (('mutation' not in data.annotations) or \
		   ('mutation' in data.annotations and \
		 'WT' in data.annotations['mutation'])):
                    if 'modifier' in data.annotations and 'SHAPE' in data.annotations['modifier']:
			if do_normalize:
			    normvals = normalize(data.values)
			else:
			    normvals = data.values
			for fragtype in frags:
			    fraglist = frags[fragtype]
			    for frag in fraglist:
				vals = []
				for idx in frag:
				    try:
					val = normvals[construct.seqpos.index(idx + construct.offset) - 1]
					if not isnan(val):
					    vals.append(val)
				    except ValueError:
					pass
				if len(vals) > 0:
				    db[fragtype].extend(vals)


pickle.dump(db, open(sys.argv[2],'w'))
