import settings
import secondary_structure

class RNA:
    def __init__(self, sequence=''):
        self.sequence = sequence

    def bootstrap(self, mapping_data, nboot, nsamples=-1, algorithm='rnastructure', replacement=False):
        print 'Starting bootstrap...'
	print 'Folding RNA with complete data'
        if nsamples < 0:
	    nsamples = len(self.sequence)
        full_bps = secondary_structure.fold(self.sequence, algorithm=algorithm, mapping_data=mapping_data)[0].base_pairs()
	bpdict = dict([(bp, 0) for bp in full_bps])
	for i in range(nboot):
	    print 'Doing bootstrap iteration %s' % i
	    md = mapping_data.sample(nsamples, replacement=replacement)
	    bps = secondary_structure.fold(self.sequence, algorithm=algorithm, mapping_data=md)[0].base_pairs()
	    for bp in bps:
		if bp in bpdict:
		    bpdict[bp] += 1
		else:
		    bpdict[bp] = 1
        for bp in bpdict:
	    bpdict[bp] *= 100./nboot
	return bpdict

    
