from numpy import array, indices, zeros
from numpy.random import randint

if __package__ is None or not __package__:
    from . import secstr
else:
    from . import secstr


class RNA(object):
    def __init__(self, sequence=''):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)


    def bootstrap(self, mapping_data, nboot, nsamples=-1, algorithm='rnastructure', bonus2d=False, is_replacement=False, fold_opts=''):
        print('Starting bootstrap...')
        print('Folding RNA with complete data')
        if bonus2d:
            mapping_data = array(mapping_data)
        nsamples = max(0, len(self.sequence))

        # print len(self.sequence), mapping_data.shape()
        full_bps = secondary_structure.fold(self.sequence, algorithm=algorithm, mapping_data=mapping_data, fold_opts='', bonus2d=bonus2d)
        # print len(full_bps)
        # print full_bps
        full_bps = full_bps[0].base_pairs()
        bp_dict = dict([(bp, 0) for bp in full_bps])
        for i in range(nboot):
            print('Doing bootstrap iteration %s' % i)
            if bonus2d:
                md = zeros(mapping_data.shape)
                N_res = mapping_data.shape[0]
                for x in range(N_res):
                    N_contrib = sum(randint(0, N_res, N_res) == x)
                    md[:, x] = mapping_data[:, x] * N_contrib
            else:
                md = mapping_data.sample(nsamples, is_replacement=is_replacement)

            struct = secondary_structure.fold(self.sequence, algorithm=algorithm, mapping_data=md, fold_opts='', bonus2d=bonus2d)[0]
            bps = struct.base_pairs()
            for bp in bps:
                if bp in bp_dict:
                    bp_dict[bp] += 1
                else:
                    bp_dict[bp] = 1

        for bp in bp_dict:
            bp_dict[bp] *= 100. / nboot

        return bp_dict


