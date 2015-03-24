import settings
import secondary_structure
from numpy import array, indices, zeros
from numpy.random import randint
from random import choice

class RNA:
    def __init__(self, sequence=''):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def bootstrap(self, mapping_data, nboot, nsamples=-1, algorithm='rnastructure', bonus2d=False, replacement=False, fold_opts=''):
        print 'Starting bootstrap...'
        print 'Folding RNA with complete data'
        if bonus2d:
            mapping_data = array(mapping_data)
        if nsamples < 0:
            nsamples = len(self.sequence)

        # print len(self.sequence), mapping_data.shape()
        full_bps = secondary_structure.fold(self.sequence, algorithm=algorithm, mapping_data=mapping_data, fold_opts='', bonus2d=bonus2d)
        # print len(full_bps)
        # print full_bps
        full_bps = full_bps[0].base_pairs()
        bpdict = dict([(bp, 0) for bp in full_bps])
        for i in range(nboot):
            print 'Doing bootstrap iteration %s' % i
            if bonus2d:
                md = zeros(mapping_data.shape)
                nres = mapping_data.shape[0]
                for x in xrange(nres):
                    ncontrib = sum(randint(0, nres, nres) ==  x)
                    md[:,x] = mapping_data[:,x]*ncontrib
                """
                grid = indices(mapping_data.shape)
                all_indices = zip(grid[0].ravel(), grid[1].ravel())
                sampled_indices = [choice(all_indices) for x in xrange(len(all_indices))]
                for j,k in sampled_indices:
                md[j,k] += mapping_data[j,k]
                """
            else:
                md = mapping_data.sample(nsamples, replacement=replacement)
            struct = secondary_structure.fold(self.sequence, algorithm=algorithm, mapping_data=md, fold_opts='', bonus2d=bonus2d)[0]
            bps = struct.base_pairs()
            for bp in bps:
                if bp in bpdict:
                    bpdict[bp] += 1
                else:
                    bpdict[bp] = 1

        for bp in bpdict:
            bpdict[bp] *= 100./nboot
        return bpdict


