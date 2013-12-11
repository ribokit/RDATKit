import pdb
import pickle
import os
import tempfile
import scipy.stats as stats
from numpy import *
from random import *
from rdatkit import settings
from rdatkit import mapping

debug = False
class SecondaryStructure:
    def __init__(self, dbn=''):
        self.dbn = dbn

    def __len__(self):
        return len(self.dbn)

    def __str__(self):
        return self.dbn

    def base_pairs(self):
        stack = []
        bps = []
        for i, s in enumerate(self.dbn):
            if s == '(':
                stack.append(i)
            if s == ')':
                bps.append((i, stack.pop()))
        return bps

    def base_pair_dict(self):
        stack = []
        bps = {}
        for i, s in enumerate(self.dbn):
            if s == '(':
                stack.append(i)
            if s == ')':
                j = stack.pop()
                bps[j] = i
                bps[i] = j
        return bps


    def helices(self):
        stack = []
        helices = []
        currhelix = []
        for i, s in enumerate(self.dbn):
            if s == '(':
                stack.append(i)
            if s == ')':
                prevbase = stack.pop()
                if len(currhelix) > 0:
                    if currhelix[-1][0] - 1 == prevbase:
                        currhelix.append((prevbase, i))
                    else:
                        helices.append(currhelix)
                        currhelix = [(prevbase, i)]
                else:
                    currhelix.append((prevbase,i))
            if s == '.':
                if len(currhelix) > 0:
                    helices.append(currhelix)
                    currhelix = []
            if len(currhelix) > 0:
                helices.append(currhelix)
        return helices


    def explode(self):
        hstack = [] # for helices
        statestack = [] # for single stranded regions
        jlist = [] # for junctions
        wstack = [] # for junction ways
        jstartstack = []
        fragments = {'helices':[], 'interiorloops':[], 'hairpins':[], 'dangles':[], 'bulges':[]}
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
                if jun[0] - 1 == pos:
                    startj = i
                    break
            for i in range(startj, len(junctions)):
                ways += 1
                junction += junctions[i]
            newjunctions = junctions[:startj]
            return junction, ways, newjunctions

        """
        """
        for h in helices:
            nhelices.append([])
            for b in h:
                nhelices[-1].append(b[0])
                nhelices[-1].append(b[1])
        fragments['helices'] = nhelices
        prevstate = '-'
        for i, s in enumerate(self.dbn):
            if s == '(':
                hstack.append(i)
                sstate = ''
                if len(statestack) > 0:
                    sstate = statestack.pop()
                if prevstate == '.' and sstate == '-':
                    fragments['dangles'].append(range(i))
                if (prevstate == '.' or prevstate == ')') and sstate != '-':
                    jlist.append(range(ssregion_start, i))
                    jstartstack.append(ssregion_start-1)
                    wstack.append(1)
                prevstate = '('
            if s == '.':
                if prevstate != '.':
                    ssregion_start = i
                    statestack.append(prevstate)
                prevstate = '.'
            if s == ')':
                prevbase = hstack.pop()
                if prevstate == '.':
                    junction, ways, jlist = _close_junction(jlist, prevbase)
                    if ways > 0:
                        for w in range(ways):
                            jstartstack.pop()
                        if ways == 1:
                            fragments['interiorloops'].append(junction + range(ssregion_start, i))
                        else:
                            key = '%dwayjunctions' % (ways + 1)
                            if key not in fragments:
                                fragments[key] = []
                            fragments[key].append(junction + range(ssregion_start, i))
                    else:
                        sstate = statestack.pop()
                        if sstate == '(':
                            fragments['hairpins'].append(range(ssregion_start, i))
                        if sstate == ')':
                            fragments['bulges'].append(range(ssregion_start, i))
                if prevstate == ')':
                    if len(jstartstack) > 0 and jstartstack[-1] == prevbase:
                        jstartstack.pop()
                        bulge = jlist.pop()
                        fragments['bulges'].append(bulge)
                        wstack.pop()
                        junction_closed = True
                    else:
                        junction, ways, jlist = _close_junction(jlist, prevbase)
                        if ways > 0:
                            for w in range(ways):
                                jstartstack.pop()
                                if ways == 1:
                                    fragments['interiorloops'].append(junction + range(ssregion_start, i))
                                else:
                                    key = '%dwayjunctions' % (ways + 1)
                                    if key not in fragments:
                                        fragments[key] = []
                                    fragments[key].append(junction + range(ssregion_start, i))
                prevstate = ')'
        if self.dbn[-1] == '.':
            fragments['dangles'].append(range(ssregion_start, len(self.dbn)))
        return fragments


    def likelihood(self, mapping_data, database='default.db.dists', db_obj=None):
        frags = self.explode()
        if not db_obj:
            db = pickle.load(open('%s/models/%s' % (settings.MAPPING_DATABASE_PATH, database)))
        else:
            db = db_obj
        data = mapping.normalize(mapping_data)
        probs = array([1.]*len(self.dbn))
        for k in frags:
            if len(frags[k]) > 0 and k in db:
                g = stats.gamma(db[k][0], db[k][1], db[k][2])
                for frag in frags[k]:
                    for i in frag:
                        if g.pdf(data[i]) > 0:
                            probs[i] = g.pdf(data[i])
        return probs, prod(probs)


def removefile(f):
    if type(f) == str:
        os.remove(os.path.abspath(f))
    else:
        os.remove(os.path.abspath(f.name))

def to_seqfile(sequence, name='placeholder'):
    seqfile = tempfile.NamedTemporaryFile(delete=False)
    seqfile.write(';\n%(name)s\n%(sequence)s1' % {'name':name, 'sequence':sequence})
    return seqfile

def _prepare_ct_and_seq_files(sequence):
    ctfile = tempfile.NamedTemporaryFile(delete=False)
    ctname = ctfile.name
    ctfile.close()
    seqfile = to_seqfile(sequence)
    seqname = seqfile.name
    seqfile.close()
    return seqname, ctname

def _prepare_fasta_file(sequence):
    fastafile = tempfile.NamedTemporaryFile(delete=False)
    fastafile.write('>bla\n%s\n' % sequence)
    return fastafile.name

def _get_fasta_structures(fname):
    structures = []
    for l in open(fname).readlines():
        l = l.strip()
        if l[0] not in  ['.', '(', ')']:
            continue
        structures.append(SecondaryStructure(dbn=l.split()[0]))
    return structures

def _get_dot_structs(ctname, nstructs, unique=False):
    structs = []
    dbns = []
    for i in range(nstructs):
        dbnfile = tempfile.NamedTemporaryFile(delete=False)
        dbnname = dbnfile.name
        dbnfile.close()
        os.popen(settings.RNA_STRUCTURE_CT2DOT + ' %s %d %s ' % \
             (ctname, i+1, dbnname))
        if debug: print(settings.RNA_STRUCTURE_CT2DOT + ' %s %d %s ' % \
             (ctname, i+1, dbnname));
        dbn = open(dbnname).readlines()[-1].strip()
        # Append only non trivial structures
        if '(' in dbn:
            if unique:
                if dbn not in dbns:
                    dbns.append(dbn)
                    structs.append(SecondaryStructure(dbn=dbn))
            else:
                structs.append(SecondaryStructure(dbn=dbn))
        removefile(dbnfile)
    return structs

def _to_ct_file(sequence, struct, filename):
    f = open(filename, 'w')
    if type(struct) != list:
	structs = [struct]
    else:
	structs = struct
    for s in structs:
        f.write('%s bla\n' % len(sequence))
	bps = s.base_pair_dict()
	for i in range(len(sequence)):
	    if i in bps:
		pair = bps[i]+1
	    else:
		pair = 0
	    if i == len(sequence) - 1:
		n = 0
	    else:
		n = i+1
	    f.write('%s %s %s %s %s %s\n' % (i+1, sequence[i], i, i+2, pair, n))
    f.close()

def _to_fasta_file(sequence, structures, fastaname):
    fastafile = open(fastaname, 'w')
    if type(structures) == list:
        for struct in structures:
            fastafile.write('>bla\n%s\n%s\n' % (sequence, str(struct)))
    else:
        fastafile.write('>bla\n%s\n%s\n' % (sequence, structures))
    fastafile.close()

def get_boltzmann_weight(sequence, structure, algorithm='rnastructure'):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        _to_ct_file(sequence, structure, ctname)
        energy = get_energies(ctname)[0]
        PARTCMD = settings.RNA_STRUCTURE_PARTITION + ' %s %s ' % (seqname, '/dev/null')
        if debug: print PARTCMD;
        ensemble_energy = float(os.popen(PARTCMD).read().split('\n')[-1])
        removefile(seqname)
        removefile(ctname)
        kT = 0.5905 #TODO Need to generalize/correct this
        weight = exp(-energy/kT)/exp(-ensemble_energy/kT)
    return weight

def mea_structure(sequence, algorithm='rnastructure', nstructs=1, gamma=1.0, opts='', returnct=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE + '/exe/MaxExpect -g %s --sequence %s %s %s' % (gamma, seqname, ctname, opts)
        if debug: print CMD;
        os.popen(CMD)
        structs = _get_dot_structs(ctname, nstructs)
        removefile(seqname)
        removefile(ctname)
    if returnct:
        return structs, ctname
    else:
        return structs

def fold(sequence, algorithm='rnastructure', mapping_data=[], nstructs=1, fold_opts='', bonus2d=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE_FOLD + ' %s %s ' % (seqname, ctname)
        tmp = None
        if len(mapping_data) > 0:
            if bonus2d:
                tmp = tempfile.NamedTemporaryFile(delete=False)
                savetxt(tmp, mapping_data)
                tmp.close()
                CMD += '-x %s ' % tmp.name
            else:
                tmp = tempfile.NamedTemporaryFile(delete=False)
                tmp.write(str(mapping_data))
                tmp.close()
                CMD += '-sh %s ' % tmp.name
                CMD += fold_opts
        if debug: print(CMD);
        os.popen(CMD)
        structs = _get_dot_structs(ctname, nstructs)
        removefile(seqname)
        removefile(ctname)
        if tmp != None:
            removefile(tmp)
    if algorithm == 'viennarna':
        fastaname = _prepare_fasta_file(sequence)
        CMD = settings.VIENNA_RNA_FOLD + ' < %s > %s.out' % (fastaname, fastaname)
        if debug: print(CMD);
        os.popen(CMD)
        structs = [SecondaryStructure(dbn=l.split()[0].strip()) for l in open(fastaname + '.out').readlines() if l.strip()[0] in ['(', ')', '.']]
        removefile(fastaname)
    
    return structs

def partition(sequence, algorithm='rnastructure', mapping_data=None, fold_opts='', bonus2d=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE_PARTITION + ' %s %s ' % (seqname, ctname)
        if mapping_data:
            CMD += _get_mapping_data_file(mapping_data, bonus2d=bonus2d)
            CMD += fold_opts
        if debug: print(CMD);
        os.popen(CMD)
    bppm = loadtxt('bpp.txt')
    for i in range(bppm.shape[0]):
        for j in range(i,bppm.shape[1]):
            if bppm[i,j] != 0:
                bppm[j,i] = bppm[i,j]
    removefile(seqname)
    removefile(ctname)
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
        ctname = ctfile.name
        _to_ct_file(sequence, structures, ctname)
        energies = get_energies(ctname, mapping_data=mapping_data)
        removefile(ctfile)
    if algorithm == 'viennarna':
        fastafile = tempfile.NamedTemporaryFile(delete=False)
        fastaname = fastafile.name
        _to_fasta_file(sequence, structures, fastaname)
        energies = get_energies(fastaname, format='fasta')
        removefile(fastafile)
    return energies

def get_energies(fname, format='ct', mapping_data=None):
    if format == 'ct':
        energyfile = tempfile.NamedTemporaryFile(delete=False)
        EFN2CMD = settings.RNA_STRUCTURE_ENERGY + ' %s %s' % (fname, energyfile.name)
        if mapping_data:
            EFN2CMD += _get_mapping_data_file(mapping_data)
        os.popen(EFN2CMD)
        energyfile.seek(0)
        energies = [float(line.split(' ')[-1]) for line in energyfile.readlines()]
        removefile(energyfile)
        return energies
    if format == 'fasta':
        RNAEVALCMD = settings.VIENNA_RNA_ENERGY + '< %s > %s.out' % (fname, fname)
        os.popen(RNAEVALCMD)
        energies = []
        for l in open(fname + '.out').readlines():
            l = l.strip()
            if l[0] not in  ['.', '(', ')']:
                continue
            energies.append(float(l.split()[-1].strip(')(')))
        removefile(fname + '.out')
        return energies

def sample(sequence, algorithm='rnastructure', mapping_data=None, nstructs=1000, unique=False, energies=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE_STOCHASTIC + ' -e %s --sequence %s %s ' % (nstructs, seqname, ctname)
        os.popen(CMD)
        structs = _get_dot_structs(ctname, nstructs, unique=unique)
        energies = get_energies(ctname)
        removefile(seqname)
        removefile(ctname)
    if energies:
        return structs, energies
    else:
        return structs



def subopt(sequence, algorithm='rnastructure', mapping_data=None, fraction=0.05, nstructs=1000, energies=False):
    if algorithm == 'rnastructure':
        seqname, ctname = _prepare_ct_and_seq_files(sequence)
        CMD = settings.RNA_STRUCTURE_ALLSUB + ' -p %d %s %s ' % (fraction*100, seqname, ctname)
        os.popen(CMD)
        structs = _get_dot_structs(ctname, nstructs)
        energies = get_energies(ctname)
        removefile(seqname)
        removefile(ctname)
    if algorithm == 'viennarna':
        fastaname = _prepare_fasta_file(sequence)
        CMD = settings.VIENNA_RNA_SUBOPT + ' -p %s < %s > %s.out' % (nstructs, fastaname, fastaname)
        if debug: print CMD;
        os.popen(CMD)
        structs = _get_fasta_structures(fastaname + '.out')
        energies = get_structure_energies(sequence, structs, algorithm='viennarna')
        removefile(fastaname)
        removefile(fastaname + '.out')
    if energies:
        return structs, energies
    return structs

def random(nstructs, length, nbp):
    structs = []
    for i in xrange(nstructs):
        dbnlist = ['.']*length
        for j in xrange(nbp):
            b2 = randint(0,length-1)
            b1 = randint(0,length-1)
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
    if factors == None:
        factors = [1]*len(structures)
    bpdict = dict([(bp, 0) for bp in ref_bps])
    for s in structures:
        for bp in bps:
            if bp in bpdict:
                bpdict[bp] += 1 * factors[s]
            else:
                bpdict[bp] = 1 * factors[s]
    for bp in bpdict:
        bpdict[bp] *= 1/len(structures)
    return bpdict

