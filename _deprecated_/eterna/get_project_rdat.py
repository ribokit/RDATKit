import argparse
import pdb

parser = argparse.ArgumentParser()

parser.add_argument('input', type=argparse.FileType('r'))
parser.add_argument('projectname', type=str)
parser.add_argument('--m2', default=False, action='store_true')
parser.add_argument('--offset', type=int, default=0)
parser.add_argument('--refseq', type=str, default='')
parser.add_argument('--maxmut', type=int, default=1)
parser.add_argument('--skip', type=int, default=0)
parser.add_argument('--exclude', type=str, default='')

args = parser.parse_args()
if args.m2:
    obligatory_tags = ['VERSION', 'COMMENT', 'ANNOTATION\t']
else:
    obligatory_tags = ['VERSION', 'SEQPOS', 'SEQUENCE', 'COMMENT', 'ANNOTATION\t']
def get_seq(line):
    for anno in line.strip().split('\t'):
        if 'sequence' in anno:
            return anno.replace('sequence:','')
                 
line = args.input.readline()
skipcount = 0
annoidx = 0
oblglines = []
reactlines = []
annolines = []
reactidces = []
offsetindices = {}
finished = False
read_prev = False
args.projectname = args.projectname.replace('\\t', '\t')
if len(args.refseq) > 0:
    wtseq = args.refseq
else:
    wtseq = ''
while line:
    addstr = ''
    for tag in obligatory_tags:
        if tag in line:
            oblglines.append(line)
            if tag == 'VERSION':
                oblglines.append('NAME\t%s\n' % args.projectname)
    if args.projectname in line and not finished:
        if len(args.exclude) == 0 or args.exclude not in line:
            if skipcount < args.skip:
                skipcount += 1
            else:
                read_prev = True
                idx = int(line.strip().replace('ANNOTATION_DATA:', '').split('\t')[0])
                reactidces.append(str(idx))
                if annoidx == 0:
                    annooffset = idx - 1
                offsetindices[idx] = idx - annooffset
                idx -= annooffset
                if args.m2:
                    writefirstmutant = False
                    if annoidx == 0:
                        if len(wtseq) == 0:
                            wtseq = get_seq(line)
                        else:
                            writefirstmutant = True
                        seqpos = ['%s%s' % (s, i+1+args.offset) for i, s in enumerate(wtseq)]
                        oblglines.append('SEQUENCE\t%s\n' % wtseq)
                    if annoidx != 0  or writefirstmutant:
                        seq = get_seq(line)
                        count = 0
                        for i in xrange(len(seq)):
                            if seq[i] != wtseq[i]:
                                addstr += '\tmutation:%s%s%s' % (wtseq[i], i+1+args.offset, seq[i])
                                count += 1
                                if count > args.maxmut:
                                    break
                annoidx += 1
                annoline = 'ANNOTATION_DATA:%d\t' % idx
                annoline += line[line.find('\t'):].strip() + addstr + '\n'
                annolines.append(annoline)
    else:
        if read_prev:
            finished = True
    if 'REACTIVITY' in line:
        if 'ERROR' in line:
            tag = 'REACTIVITY_ERROR'
        else:
            tag = 'REACTIVITY'
        for idx in reactidces:
            if tag+':'+idx + ' ' in line or tag + ':' + '\t' in line:
                reacts = line.strip().replace('%s:%s' % (tag, idx), '').strip().split('\t')
                
                if args.m2:
                    reactlines.append('%s:%s\t%s\n' % (tag, offsetindices[int(idx)], '\t'.join([reacts[i] for i in range(len(seqpos))])))
                else:
                    reactlines.append('%s:%s\t%s\n' % (tag, offsetindices[int(idx)], '\t'.join(reacts)))
    if args.m2:
        if 'SEQPOS' in line:
            inseqpos = [x[1:] for x in line.replace('SEQPOS\t', '').strip().split('\t')]
            seqpos = [s for i, s in enumerate(seqpos) if i < len(inseqpos)]
            oblglines.append('SEQPOS\t%s\n' % '\t'.join(seqpos))

    if 'OFFSET' in line:
        oblglines.append('OFFSET\t%s\n' % args.offset)
    line = args.input.readline()

def writelines(lines):
    for line in lines:
        print line.strip() + '\n'

writelines(oblglines)
writelines(annolines)
writelines(reactlines)
