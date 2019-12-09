from numpy import *
import os
import pickle
from random import *
from scipy.stats import gamma
import subprocess
import tempfile

if __package__ is None or not __package__:
    from .util import *
    from . import mapping
else:
    from .util import *
    from . import mapping

debug = False


class SecondaryStructure(object):

    def __init__(self, dbn=''):
        self.dbn = dbn

    def __len__(self):
        return len(self.dbn)

    def __str__(self):
        return self.dbn


    def _get_base_pairs(self, l_char, r_char):
        stack = []
        bps = []
        for i, s in enumerate(self.dbn):
            if s == l_char:
                stack.append(i)
            elif s == r_char:
                bps.append((i, stack.pop()))
        return bps

    def base_pairs(self):
        bps = self._get_base_pairs('(', ')')
        bps += self._get_base_pairs('{', '}')
        bps += self._get_base_pairs('[', ']')
        return bps


    def _get_base_pair_dict(self, l_char, r_char):
        stack = []
        bps = {}
        for i, s in enumerate(self.dbn):
            if s == l_char:
                stack.append(i)
            elif s == r_char:
                j = stack.pop()
                bps[j] = i
                bps[i] = j
        return bps

    def base_pair_dict(self):
        bps = self._get_base_pair_dict('(', ')')
        bps.update(self._get_base_pair_dict('{', '}'))
        bps.update(self._get_base_pair_dict('[', ']'))
        return bps


    def _get_helices(self, l_char, r_char):
        stack = []
        helices = []
        curr_helix = []
        for i, s in enumerate(self.dbn):
            if s == l_char:
                stack.append(i)
            if s == r_char:
                prev_base = stack.pop()
                if len(curr_helix):
                    if curr_helix[-1][0] == prev_base + 1 and curr_helix[-1][1] == i - 1:
                        curr_helix.append((prev_base, i))
                    else:
                        helices.append(curr_helix)
                        curr_helix = [(prev_base, i)]
                else:
                    curr_helix.append((prev_base, i))
            if s == '.':
                if len(curr_helix):
                    helices.append(curr_helix)
                    curr_helix = []
        if len(curr_helix):
            helices.append(curr_helix)
        return helices

    def helices(self):
        helices = self._get_helices('(', ')')
        helices += self._get_helices('{', '}')
        helices += self._get_helices('[', ']')
        return helices


    def explode(self):
        hstack = []  # for helices
        statestack = []  # for single stranded regions
        jlist = []  # for junctions
        wstack = []  # for junction ways
        jstartstack = []
        fragments = {'helices': [], 'interiorloops': [], 'hairpins': [], 'dangles': [], 'bulges': []}
        helices = self.helices()
        nhelices = []
        """
        For junction closing
        """
        def _close_junction(junctions, pos):
            junction = []
            ways = 0
            startj = len(junctions)
            for i, jun in enumerate(junctions):
                if len(jun) and jun[0] - 1 == pos:
                    startj = i
                    break
            for i in range(startj, len(junctions)):
                ways += 1
                junction += [x for x in junctions[i] if self.dbn[x] == '.' or self.dbn[x].lower() == 'a']
            newjunctions = junctions[:startj]
            return (junction, ways, newjunctions)

        def _get_ssregion(start, i):
            if start < 0:
                return []
            else:
                return list(range(start, i))


        for h in helices:
            nhelices.append([])
            for b in h:
                nhelices[-1].append(b[0])
                nhelices[-1].append(b[1])
        fragments['helices'] = nhelices

        prev_state = '-'
        for i, s in enumerate(self.dbn):
            if s == '(':
                hstack.append(i)
                s_state = ''
                if len(statestack):
                    s_state = statestack.pop()
                if prev_state == '.' and s_state == '-':
                    fragments['dangles'].append(list(range(i)))
                if s_state != '-':
                    if (prev_state == '.' or prev_state == ')') and s_state != '-':
                        jlist.append(_get_ssregion(ss_region_start, i))
                        jstartstack.append(ss_region_start-1)
                        wstack.append(1)
                        ss_region_start = -1
                prev_state = '('
            elif s == '.' or s == 'a':
                if prev_state != '.':
                    ss_region_start = i
                    statestack.append(prev_state)
                prev_state = '.'
            elif s == ')':
                prev_base = hstack.pop()
                if prev_state == '.':
                    (junction, ways, jlist) = _close_junction(jlist, prev_base)
                    if ways > 0:
                        for w in range(ways):
                            jstartstack.pop()
                        if ways == 1:
                            fragments['interiorloops'].append(junction + _get_ssregion(ss_region_start, i))
                        else:
                            key = '%dwayjunctions' % (ways + 1)
                            if key not in fragments:
                                fragments[key] = []
                            fragments[key].append(junction + _get_ssregion(ss_region_start, i))
                    else:
                        s_state = statestack.pop()
                        if s_state == '(':
                            fragments['hairpins'].append(_get_ssregion(ss_region_start, i))
                        if s_state == ')':
                            fragments['bulges'].append(_get_ssregion(ss_region_start, i))
                    ss_region_start = -1

                elif prev_state == ')':
                    if len(jstartstack) and jstartstack[-1] == prev_base:
                        jstartstack.pop()
                        bulge = jlist.pop()
                        fragments['bulges'].append(bulge)
                        wstack.pop()
                        junction_closed = True
                    else:
                        (junction, ways, jlist) = _close_junction(jlist, prev_base)
                        if ways > 0:
                            for w in range(ways):
                                jstartstack.pop()
                                if ways == 1:
                                    fragments['interiorloops'].append(junction + _get_ssregion(ss_region_start, i))
                                else:
                                    key = '%dwayjunctions' % (ways + 1)
                                    if key not in fragments:
                                        fragments[key] = []
                                    fragments[key].append(junction + _get_ssregion(ss_region_start, i))
                                ss_region_start = -1
                prev_state = ')'

        if self.dbn[-1] == '.':
            fragments['dangles'].append(_get_ssregion(ss_region_start, len(self.dbn)))

        # Eliminate spurious adds for single stranded regions and re-assign aptamer regions
        aptamers = []
        for k, v in fragments.items():
            if k != 'helices':
                for i, ntlist in enumerate(v):
                    v[i] = [x for x in ntlist if self.dbn[x] == '.']
                    apt = [x for x in ntlist if self.dbn[x].lower() == 'a']
                    if len(apt):
                        aptamers.append(apt)
            # Get rid of empty lists
            fragments[k] = [x for x in v if len(x)]
        fragments['aptamers'] = aptamers

        prev_i = -1
        sstrand_region = []
        fragments['sstrand'] = []
        for i in range(len(self.dbn)):
            found = False
            for k, ntlists in fragments.items():
                for ntlist in ntlists:
                    if i in ntlist:
                        found = True
                        break

            if not found:
                if i != prev_i + 1 and prev_i >= 0:
                    fragments['sstrand'].append(sstrand_region)
                    sstrand_region = []
                sstrand_region.append(i)
                prev_i = i

        if len(sstrand_region):
            fragments['sstrand'].append(sstrand_region)

        return fragments


    def likelihood(self, mapping_data, database='default.db.dists', db_obj=None):
        frags = self.explode()
        if not db_obj:
            db = pickle.load(open('%s/models/%s' % (MAPPING_DATABASE_PATH, database)))
        else:
            db = db_obj

        data = mapping.normalize(mapping_data)
        probs = array([1.] * len(self.dbn))
        for k in frags:
            if len(frags[k]) and k in db:
                g = gamma(db[k][0], db[k][1], db[k][2])
                for frag in frags[k]:
                    for i in frag:
                        if g.pdf(data[i]) > 0:
                            probs[i] = g.pdf(data[i])

        return probs, prod(probs)


def remove_file(f):
    if isinstance(f, str):
        os.remove(os.path.abspath(f))
    else:
        os.remove(os.path.abspath(f.name))

def to_seqfile(sequence, name='placeholder'):
    seqfile = tempfile.NamedTemporaryFile(delete=False)
    seqfile.write(';\n%(name)s\n%(sequence)s1' % {'name': name, 'sequence': sequence})
    return seqfile


def _prepare_ct_and_seq_files(sequence):
    ctfile = tempfile.NamedTemporaryFile(delete=False)
    ct_name = ctfile.name
    ctfile.close()
    seqfile = to_seqfile(sequence)
    seq_name = seqfile.name
    seqfile.close()
    return (seq_name, ct_name)

def _prepare_fasta_file(sequence):
    fasta_file = tempfile.NamedTemporaryFile(delete=False)
    fasta_file.write('>bla\n%s\n' % sequence)
    return fasta_file.name

def _get_fasta_structures(fname):
    structures = []
    for l in open(fname).readlines():
        l = l.strip()
        if l[0] not in ['.', '(', ')']:
            continue
        structures.append(SecondaryStructure(dbn=l.split()[0]))
    return structures

def _get_dot_structs(ct_name, N_structs, unique=False):
    structs = []
    dbns = []
    for i in range(N_structs):
        dbnfile = tempfile.NamedTemporaryFile(delete=False)
        dbnname = dbnfile.name
        dbnfile.close()
        subprocess.check_call(PATH_RNA_STRUCTURE_CT2DOT + ' %s %d %s ' % (ct_name, i + 1, dbnname), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        dbn = open(dbnname).readlines()[-1].strip()
        # Append only non trivial structures
        if '(' in dbn:
            if unique:
                if dbn not in dbns:
                    dbns.append(dbn)
                    structs.append(SecondaryStructure(dbn=dbn))
            else:
                structs.append(SecondaryStructure(dbn=dbn))
        remove_file(dbnfile)

    return structs


def _to_ct_file(sequence, struct, filename):
    f = open(filename, 'w')
    if isinstance(struct, list):
        structs = struct
    else:
        structs = [struct]

    for s in structs:
        f.write('%s bla\n' % len(sequence))
        bps = s.base_pair_dict()
        for i in range(len(sequence)):
            pair = bps[i] + 1 if i in bps else 0
            n = 0 if i == len(sequence) - 1 else i + 1
            f.write('%s %s %s %s %s %s\n' % (i+1, sequence[i], i, i+2, pair, n))
    f.close()

def _to_fasta_file(sequence, structures, fastaname):
    fasta_file = open(fastaname, 'w')

    if isinstance(structures, list):
        for struct in structures:
            fasta_file.write('>bla\n%s\n%s\n' % (sequence, str(struct)))
    else:
        fasta_file.write('>bla\n%s\n%s\n' % (sequence, structures))
    fasta_file.close()


def get_boltzmann_weight(sequence, structure, algorithm='rnastructure'):
    if algorithm == 'rnastructure':
        (seq_name, ct_name) = _prepare_ct_and_seq_files(sequence)
        _to_ct_file(sequence, structure, ct_name)

        energy = get_energies(ct_name)[0]
        cmd = PATH_RNA_STRUCTURE_PARTITION + ' %s %s ' % (seq_name, '/dev/null')
        ensemble_energy = float(subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].strip().split('\n')[-1])

        remove_file(seq_name)
        remove_file(ct_name)
        kT = 0.5905  # TODO Need to generalize/correct this
        weight = exp(-energy / kT) / exp(-ensemble_energy / kT)
    return weight


def mea_structure(sequence, algorithm='rnastructure', N_structs=1, gamma=1.0, opts='', returnct=False):
    if algorithm == 'rnastructure':
        (seq_name, ct_name) = _prepare_ct_and_seq_files(sequence)
        cmd = PATH_RNA_STRUCTURE + 'exe/MaxExpect -g %s --sequence %s %s %s' % (gamma, seq_name, ct_name, opts)
        subprocess.check_call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        structs = _get_dot_structs(ct_name, N_structs)
        remove_file(seq_name)
        remove_file(ct_name)

    if returnct:
        return structs, ct_name
    else:
        return structs


def fold(sequence, algorithm='rnastructure', modifier='shape', mapping_data=[], N_structs=1, fold_opts='', bonus2d=False):
    if algorithm == 'rnastructure':
        (seq_name, ct_name) = _prepare_ct_and_seq_files(sequence)
        CMD = PATH_RNA_STRUCTURE_FOLD + ' %s %s ' % (seq_name, ct_name)
        tmp = None

        if len(mapping_data) > 0:
            if bonus2d:
                tmp = tempfile.NamedTemporaryFile(delete=False)
                savetxt(tmp, mapping_data)
                tmp.close()
                CMD += '-x %s ' % tmp.name
            else:
                tmp = tempfile.NamedTemporaryFile(delete=False)
                if modifier.lower() == 'dms':
                    for pos in mapping_data.seqpos:
                        s = sequence[pos].lower()
                        if s == 'a' or s == 'c':
                            tmp.write('%d %.3f\n' % (pos + 1, float(mapping_data[pos])))
                    CMD += '-dms %s ' % tmp.name
                elif modifier.lower() == 'cmct':
                    for pos in mapping_data.seqpos:
                        s = sequence[pos].lower()
                        if s == 'g' or s == 'u':
                            tmp.write('%d %.3f\n' % (pos + 1, float(mapping_data[pos])))
                    CMD += '-cmct %s ' % tmp.name
                else:
                    CMD += '-sh %s ' % tmp.name
                    tmp.write(str(mapping_data))
                tmp.close()
        CMD += fold_opts
        subprocess.check_call(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        structs = _get_dot_structs(ct_name, N_structs)
        remove_file(seq_name)
        remove_file(ct_name)
        if tmp is not None:
            remove_file(tmp)

    elif algorithm == 'viennarna':
        fastaname = _prepare_fasta_file(sequence)
        CMD = PATH_VIENNA_RNA_FOLD + ' < %s > %s.out' % (fastaname, fastaname)
        subprocess.check_call(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        structs = [SecondaryStructure(dbn=l.split()[0].strip()) for l in open(fastaname + '.out').readlines() if l.strip()[0] in ['(', ')', '.']]
        remove_file(fastaname)

    return structs


def partition(sequence, algorithm='rnastructure', mapping_data=None, fold_opts='', bonus2d=False):
    if algorithm == 'rnastructure':
        (seq_name, ct_name) = _prepare_ct_and_seq_files(sequence)
        CMD = PATH_RNA_STRUCTURE_PARTITION + ' %s %s ' % (seq_name, ct_name)
        if mapping_data:
            CMD += _get_mapping_data_file(mapping_data, bonus2d=bonus2d)
            CMD += fold_opts
        subprocess.check_call(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    bppm = loadtxt('bpp.txt')
    for i in range(bppm.shape[0]):
        for j in range(i, bppm.shape[1]):
            if bppm[i, j] != 0:
                bppm[j, i] = bppm[i, j]
    remove_file(seq_name)
    remove_file(ct_name)
    return bppm


def _get_mapping_data_file(mapping_data, bonus2d=False):
    tmp = tempfile.NamedTemporaryFile(delete=False)
    if bonus2d:
        savetxt(tmp, mapping_data)
        tmp.close()
        opt = ' -x %s ' % tmp.name
    else:
        tmp.write(str(mapping_data))
        tmp.close()
        opt = ' -sh %s ' % tmp.name
    return opt


def get_structure_energies(sequence, structures, mapping_data=None, algorithm='rnastructure'):
    if algorithm == 'rnastructure':
        energies = []
        ctfile = tempfile.NamedTemporaryFile(delete=False)
        ct_name = ctfile.name
        _to_ct_file(sequence, structures, ct_name)
        energies = get_energies(ct_name, mapping_data=mapping_data)
        remove_file(ctfile)
    if algorithm == 'viennarna':
        fasta_file = tempfile.NamedTemporaryFile(delete=False)
        fastaname = fasta_file.name
        _to_fasta_file(sequence, structures, fastaname)
        energies = get_energies(fastaname, format='fasta')
        remove_file(fasta_file)
    return energies


def get_energies(fname, format='ct', mapping_data=None):
    if format == 'ct':
        energyfile = tempfile.NamedTemporaryFile(delete=False)
        EFN2CMD = PATH_RNA_STRUCTURE_ENERGY + ' %s %s' % (fname, energyfile.name)
        if mapping_data:
            EFN2CMD += _get_mapping_data_file(mapping_data)
        subprocess.check_call(EFN2CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        energyfile.seek(0)
        energies = [float(line.split(' ')[-1]) for line in energyfile.readlines()]
        remove_file(energyfile)
        return energies

    elif format == 'fasta':
        RNAEVALCMD = PATH_VIENNA_RNA_ENERGY + '< %s > %s.out' % (fname, fname)
        subprocess.check_call(RNAEVALCMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        energies = []
        for l in open(fname + '.out').readlines():
            l = l.strip()
            if l[0] not in ['.', '(', ')']:
                continue
            energies.append(float(l.split()[-1].strip(')(')))
        remove_file(fname + '.out')
        return energies


def sample(sequence, algorithm='rnastructure', mapping_data=None, N_structs=1000, unique=False, energies=False):
    if algorithm == 'rnastructure':
        (seq_name, ct_name) = _prepare_ct_and_seq_files(sequence)
        CMD = PATH_RNA_STRUCTURE_STOCHASTIC + ' -e %s --sequence %s %s ' % (N_structs, seq_name, ct_name)
        subprocess.check_call(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        structs = _get_dot_structs(ct_name, N_structs, unique=unique)
        energies = get_energies(ct_name)
        remove_file(seq_name)
        remove_file(ct_name)

    if energies:
        return structs, energies
    else:
        return structs



def subopt(sequence, algorithm='rnastructure', mapping_data=None, fraction=0.05, N_structs=1000, energies=False):
    if algorithm == 'rnastructure':
        (seq_name, ct_name) = _prepare_ct_and_seq_files(sequence)
        CMD = PATH_RNA_STRUCTURE_ALLSUB + ' -p %d %s %s ' % (fraction*100, seq_name, ct_name)
        subprocess.check_call(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        structs = _get_dot_structs(ct_name, N_structs)
        energies = get_energies(ct_name)
        remove_file(seq_name)
        remove_file(ct_name)
    elif algorithm == 'viennarna':
        fastaname = _prepare_fasta_file(sequence)
        CMD = PATH_VIENNA_RNA_SUBOPT + ' -p %s < %s > %s.out' % (N_structs, fastaname, fastaname)
        subprocess.check_call(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        structs = _get_fasta_structures(fastaname + '.out')
        energies = get_structure_energies(sequence, structs, algorithm='viennarna')
        remove_file(fastaname)
        remove_file(fastaname + '.out')
    if energies:
        return structs, energies
    return structs


def random(N_structs, length, nbp):
    structs = []
    for i in range(N_structs):
        dbnlist = ['.'] * length
        for j in range(nbp):
            b2 = randint(0, length - 1)
            b1 = randint(0, length - 1)
            if b1 <= b2:
                dbnlist[b1] = '('
                dbnlist[b2] = ')'
            else:
                dbnlist[b2] = '('
                dbnlist[b1] = ')'
        structs.append(SecondaryStructure(dbn=''.join(dbnlist)))
    return structs


def base_pair_fractions_in_structures(reference, structures, factors=None):
    ref_bp = reference.base_pairs()
    if factors is None:
        factors = [1] * len(structures)

    bp_dict = dict([(bp, 0) for bp in ref_bp])
    for i, s in enumerate(structures):
        bps = s.base_pairs()
        for bp in bps:
            if bp in bp_dict:
                bp_dict[bp] = bp_dict[bp] + factors[i]
            if bp[::-1] in bp_dict:
                bp_dict[bp[::-1]] = bp_dict[bp[::-1]] + factors[i]

    for bp in list(bp_dict.keys()):
        if factors is None:
            bp_dict[bp] *= 100. / len(structures)
        else:
            bp_dict[bp] *= 100.

    for bp in ref_bp:
        bp_dict[bp[::-1]] = bp_dict[bp]

    return bp_dict

