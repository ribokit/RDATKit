"""
datahandlers is a module that contains classes for representing structure mapping file formats in python.
@author Pablo Sanchez Cordero
"""
from collections import defaultdict
import copy
# from itertools import chain
import os
import xlwt
import xlrd

if __package__ is None or not __package__:
    from util import *
else:
    from .util import *


def split(s, delims=None):
    if ',' in delims:
        delims = delims.split(',')
    process = False
    for d in delims:
        if d in s:
            process = True
            break

    if process:
        for d in delims[1:]:
            s = s.replace(d, delims[0])
        return [x for x in s.split(delims[0]) if x]
    return [s]


"""
For the RDAT file format for sharing chemical footprinting experiments
"""

class RDATSection(object):
    def __init__(self, attr_list=[], attr_str=[]):
        for attr in attr_list:
            setattr(self, attr, [])
        for attr in attr_str:
            setattr(self, attr, '')


class RDATFile(object):
    def __init__(self):
        self.constructs = defaultdict(list)

        self.values = defaultdict(list)
        self.errors = defaultdict(list)
        self.traces = defaultdict(list)
        self.reads = defaultdict(list)

        self.comments = ''
        self.annotations = defaultdict(list)
        self.data_types = defaultdict(list)
        self.mutpos = defaultdict(list)
        self.xsels = defaultdict(list)

        self.filename = None
        self.version = None
        self.loaded = False


    def _append_new_data_section(self, this_construct):
        d = RDATSection(['seqpos', 'values', 'errors', 'trace', 'reads', 'xsel'])
        d.annotations = {}
        self.constructs[this_construct].data.append(d)

    def _parse_data_block(self, line, key, start_idx=0):
        if not isinstance(key, list):
            key = [key]
        for k in key:
            attheader = k + ':' if ':' in line else k
            line = line.replace(attheader, '')

        fields = split(line.strip('\n ,'), delims='\t, ')
        data_idx = int(fields[0]) - 1 if start_idx else None
        data = [float(x) if ':' not in x else float(x[:x.find(':')]) for x in fields[start_idx:]]
        return (data, data_idx)

    def _parse_annotations(self, s):
        d = {}
        if self.version == 0.1:
            token = ';'
            s = s.split(',')
        else:
            token = ':'

        for item in s:
            if item:
                pair = item.split(token)
                if pair[0].strip() in d:
                    d[pair[0].strip()].append(':'.join(pair[1:]))
                else:
                    d[pair[0].strip()] = [':'.join(pair[1:])]
        return d

    def _annotation_str(self, a, delim):
        s = ''
        for k in a:
            if isinstance(a[k], list):
                for i in range(len(a[k])):
                    s += k + ':' + a[k][i] + delim
            else:
                s += k + ':' + str(a[k]) + delim
        return s


    def load(self, file):
        self.filename = file.name

        # only used for self.version == 0.1:
        current_section = 'general'
        fill_data_types = False
        # data_dict = {}

        lines = file.readlines()
        for line in lines:
            line = line.strip()

            if line.startswith('VERSION:'):
                self.version = float(line.replace('VERSION:', ''))
                continue
            elif line.startswith('RDAT_VERSION') or line.startswith('VERSION'):
                self.version = float(line.replace('RDAT_VERSION', '').replace('VERSION', ''))
                continue

            if self.version == 0.1:
                if line.startswith('COMMENTS:'):
                    self.comments = line.replace('COMMENTS:', '').strip()

                elif line.startswith('ANNOTATION:'):
                    if current_section == 'general':
                        self.annotations = self._parse_annotations(line.replace('ANNOTATION:', ''))
                    elif current_section == 'construct':
                        annotations = self._parse_annotations(line.replace('ANNOTATION:', ''))
                        self.constructs[this_construct].annotations = annotations

                        if 'modifier' in annotations:
                            self.data_types[this_construct].append(annotations['modifier'][0])
                            fill_data_types = True
                    elif current_section == 'data':
                        annotations = self._parse_annotations(line.replace('ANNOTATION:', ''))
                        self.constructs[this_construct].data[data_idx].annotations = annotations

                        if 'modifier' in annotations:
                            self.data_types[this_construct].append(annotations['modifier'][0])
                        if 'mutation' in annotations:
                            try:
                                self.mutpos[this_construct][-1] = int(annotations['mutation'][0][1:-1])
                            except ValueError:
                                pass
                    else:
                        print('Attribute :%s does not belong to a valid section' % line)

                elif 'CONSTRUCT' in line:
                    current_section = 'construct'
                    if fill_data_types:
                        self.data_types[this_construct] = [self.data_types[this_construct][0]] * len(self.values[this_construct])
                        fill_data_types = False

                    this_construct = file.readline().strip().replace('NAME:', '').strip()
                    data_idx = -1
                    self.constructs[this_construct] = RDATSection(['seqpos', 'data', 'xsel'], ['structure'])
                    self.constructs[this_construct].name = this_construct
                    self.constructs[this_construct].structure = ''

                elif line.startswith('SEQUENCE:'):
                    self.constructs[this_construct].sequence = line.replace('SEQUENCE:', '').strip()

                elif line.startswith('STRUCTURE:'):
                    self.constructs[this_construct].structure = line.replace('STRUCTURE:', '').strip()

                elif line.startswith('WELLS:'):
                    self.constructs[this_construct].wells = line.replace('WELLS:', '').strip().split(',')

                elif line.startswith('OFFSET:'):
                    self.constructs[this_construct].offset = int(line.replace('OFFSET:', ''))

                elif line.startswith('DATA'):
                    current_section = 'data'
                    data_idx += 1
                    d = RDATSection(['seqpos', 'xsel'])
                    self.mutpos[this_construct].append('WT')
                    self.constructs[this_construct].data.append(d)

                elif line.startswith('SEQPOS:'):
                    if current_section == 'construct':
                        self.constructs[this_construct].seqpos = [int(x) for x in line.replace('SEQPOS:', '').strip(' ,').split(',')]
                    else:
                        self.constructs[this_construct].data[data_idx].seqpos = [int(x) for x in line.replace('SEQPOS:', '').strip(' ,').split(',')]

                elif line.startswith('VALUES'):
                    self.constructs[this_construct].data[data_idx].values = [float(x) for x in line.replace('VALUES:', '').strip(' ,').split(',')]
                    self.values[this_construct].append(self.constructs[this_construct].data[data_idx].values)

                elif line.startswith('TRACE'):
                    self.constructs[this_construct].data[data_idx].trace = [float(x) for x in line.replace('TRACE:', '').strip(' ,').split(',')]
                    self.traces[this_construct].append(self.constructs[this_construct].data[data_idx].trace)

                elif line.startswith('XSEL:'):
                    if current_section == 'construct':
                        self.constructs[this_construct].xsel = [float(x) for x in line.replace('XSEL:', '').strip(' ,').split(',')]
                    else:
                        self.constructs[this_construct].data[data_idx].xsel = [float(x) for x in line.replace('XSEL:', '').strip(' ,').split(',')]
                        self.xsels[this_construct].append(self.constructs[this_construct].data[data_idx].xsel)
                else:
                    if line.strip():
                        raise AttributeError('Invalid section: ' + line)

            elif self.version >= 0.2 and self.version < 0.4:
                if line.startswith('COMMENT'):
                    parsed_line = line
                    for sep in ' \t':
                        parsed_line = parsed_line.replace('COMMENTS' + sep, '').replace('COMMENT' + sep, '')
                    self.comments += parsed_line + '\n'

                elif line.startswith('ANNOTATION') and not line.startswith('ANNOTATION_DATA'):
                    self.annotations = self._parse_annotations(split(line.replace('ANNOTATION', ''), delims='\t'))

                elif 'CONSTRUCT' in line or line.startswith('NAME'):
                    if 'CONSTRUCT' in line:
                        line = file.readline().strip()  # Advance to 'NAME' line.

                    this_construct = line.replace('NAME', '').strip()
                    data_idx = -1
                    self.constructs[this_construct] = RDATSection(['seqpos', 'data', 'xsel'], ['sequence', 'structure'])
                    self.constructs[this_construct].name = this_construct
                    self.constructs[this_construct].annotations = {}
                    self.constructs[this_construct].structures = defaultdict(str)
                    self.constructs[this_construct].sequences = defaultdict(str)

                elif line.startswith('SEQUENCE'):
                    attheader = 'SEQUENCE:' if ':' in line else 'SEQUENCE'
                    line = line.replace(attheader, '')
                    if len(line.split()) > 1:
                        seqidx, seq = line.strip().split()
                        self.constructs[this_construct].sequences[int(seqidx)] = seq.strip()
                        self.constructs[this_construct].sequence = seq.strip()
                    else:
                        seq = line
                        self.constructs[this_construct].sequence = seq.strip()
                        self.constructs[this_construct].sequences[0] = seq.strip()

                elif line.startswith('STRUCTURE'):
                    attheader = 'STRUCTURE:' if ':' in line else 'STRUCTURE'
                    line = line.replace(attheader, '')
                    if len(line.split()) > 1:
                        structidx, struct = line.strip().split()
                        self.constructs[this_construct].structures[int(structidx)] = struct.strip()
                        self.constructs[this_construct].structure = struct.strip()
                    else:
                        struct = line
                        self.constructs[this_construct].structure = struct.strip()
                        self.constructs[this_construct].structures[0] = struct.strip()

                elif line.startswith('OFFSET'):
                    self.constructs[this_construct].offset = int(line.replace('OFFSET', ''))

                elif line.startswith('DATA_TYPE'):
                    self.data_types[this_construct] = split(line.replace('DATA_TYPE', '').strip(), delims='\t')

                elif line.startswith('SEQPOS'):
                    seqpos_tmp = split(line.replace('SEQPOS', '').strip(), delims='\t, ')
                    if self.version >= 0.32:
                        self.constructs[this_construct].seqpos = [int(x[1:]) for x in seqpos_tmp]
                    else:
                        self.constructs[this_construct].seqpos = [int(x) for x in seqpos_tmp]

                elif line.startswith('MUTPOS'):
                    self.mutpos[this_construct] = [x.strip() for x in split(line.replace('MUTPOS', '').strip(), delims='\t')]

                elif line.startswith('ANNOTATION_DATA'):
                    fields = split(line.replace('ANNOTATION_DATA:', '').replace('ANNOTATION_DATA ', '').strip(), delims='\t')
                    if len(fields) < 2:
                        fields = split(fields[0], delims=' ')
                    data_idx = int(fields[0]) - 1
                    annotations = self._parse_annotations(fields[1:])
                    for l in range(data_idx - len(self.constructs[this_construct].data) + 1):
                        self._append_new_data_section(this_construct)
                    self.constructs[this_construct].data[data_idx].annotations = annotations
                    if 'mutation' in annotations:
                        try:
                            if len(self.mutpos[this_construct]) > 0:
                                self.mutpos[this_construct][-1] = int(annotations['mutation'][0][1:-1])
                            else:
                                self.mutpos[this_construct].append(int(annotations['mutation'][0][1:-1]))
                        except ValueError:
                            pass

                elif line.startswith('AREA_PEAK') or line.startswith('REACTIVITY:'):
                    (peaks, data_idx) = self._parse_data_block(line, ['AREA_PEAK', 'REACTIVITY'], 1)
                    if (data_idx >= len(self.constructs[this_construct].data)):
                        self._append_new_data_section(this_construct)
                    self.constructs[this_construct].data[data_idx].values = peaks
                    self.values[this_construct].append(self.constructs[this_construct].data[data_idx].values)

                elif line.startswith('AREA_PEAK_ERROR') or line.startswith('REACTIVITY_ERROR:'):
                    (errors, data_idx) = self._parse_data_block(line, ['AREA_PEAK_ERROR', 'REACTIVITY_ERROR'], 1)
                    self.constructs[this_construct].data[data_idx].errors = errors
                    self.errors[this_construct].append(self.constructs[this_construct].data[data_idx].errors)

                elif line.startswith('TRACE'):
                    (trace, data_idx) = self._parse_data_block(line, 'TRACE', 1)
                    if data_idx < len(self.constructs[this_construct].data):
                        self.constructs[this_construct].data[data_idx].trace = trace
                        self.traces[this_construct].append(self.constructs[this_construct].data[data_idx].trace)

                elif line.startswith('READS'):
                    (reads, data_idx) = self._parse_data_block(line, 'READS', 1)
                    if data_idx < len(self.constructs[this_construct].data):
                        self.constructs[this_construct].data[data_idx].reads = reads
                        self.reads[this_construct].append(self.constructs[this_construct].data[data_idx].reads)

                elif line.startswith('XSEL_REFINE'):
                    (xsel, data_idx) = self._parse_data_block(line, 'XSEL_REFINE', 1)
                    self.constructs[this_construct].data[data_idx].xsel = xsel
                    self.xsels[this_construct].append(self.constructs[this_construct].data[data_idx].xsel)

                elif line.startswith('XSEL'):
                    (xsel, data_idx) = self._parse_data_block(line, 'XSEL', 0)
                    self.constructs[this_construct].xsel = xsel

                else:
                    if line.strip():
                        raise AttributeError('Invalid section: ' + line)

            elif self.version >= 0.4:
                if line.startswith('COMMENT'):
                    parsed_line = line
                    for sep in ' \t':
                        parsed_line = parsed_line.replace('COMMENTS' + sep, '').replace('COMMENT' + sep, '')
                    self.comments += parsed_line + '\n'

                elif line.startswith('ANNOTATION'):
                    self.annotations = self._parse_annotations(split(line.replace('ANNOTATION', ''), delims='\t'))

                elif line.startswith('NAME'):
                    this_construct = line.replace('NAME', '').strip()
                    data_idx = -1
                    self.constructs[this_construct] = RDATSection(['seqpos', 'data', 'xsel'], ['sequence', 'structure'])
                    self.constructs[this_construct].name = this_construct
                    self.constructs[this_construct].annotations = {}
                    self.constructs[this_construct].structures = defaultdict(str)
                    self.constructs[this_construct].sequences = defaultdict(str)

                elif line.startswith('SEQUENCE'):
                    attheader = 'SEQUENCE:' if ':' in line else 'SEQUENCE'
                    line = line.replace(attheader, '')
                    self.constructs[this_construct].sequence = line.strip()
                    self.constructs[this_construct].sequences[0] = line.strip()

                elif line.startswith('STRUCTURE'):
                    attheader = 'STRUCTURE:' if ':' in line else 'STRUCTURE'
                    line = line.replace(attheader, '')
                    self.constructs[this_construct].structure = line.strip()
                    self.constructs[this_construct].structures[0] = line.strip()

                elif line.startswith('OFFSET'):
                    self.constructs[this_construct].offset = int(line.replace('OFFSET', ''))

                elif line.startswith('SEQPOS'):
                    seqpos_tmp = split(line.replace('SEQPOS', '').strip(), delims='\t, ')
                    self.constructs[this_construct].seqpos = [int(x[1:]) for x in seqpos_tmp]

                elif line.startswith('DATA_ANNOTATION:'):
                    fields = split(line.replace('DATA_ANNOTATION:', '').strip(), delims='\t')

                    if len(fields) < 2:
                        fields = split(fields[0], delims=' ')
                    data_idx = int(fields[0]) - 1
                    annotations = self._parse_annotations(fields[1:])
                    for l in range(data_idx - len(self.constructs[this_construct].data) + 1):
                        self._append_new_data_section(this_construct)
                    self.constructs[this_construct].data[data_idx].annotations = annotations

                elif line.startswith('DATA:'):
                    (data, data_idx) = self._parse_data_block(line, 'DATA', 1)
                    if (data_idx >= len(self.constructs[this_construct].data)):
                        self._append_new_data_section(this_construct)
                    self.constructs[this_construct].data[data_idx].values = data
                    self.values[this_construct].append(self.constructs[this_construct].data[data_idx].values)

                elif line.startswith('XSEL_REFINE'):
                    (xsel, data_idx) = self._parse_data_block(line, 'XSEL_REFINE', 1)
                    self.constructs[this_construct].data[data_idx].xsel = xsel
                    self.xsels[this_construct].append(self.constructs[this_construct].data[data_idx].xsel)

                elif line.startswith('XSEL'):
                    (xsel, data_idx) = self._parse_data_block(line, 'XSEL', 0)
                    self.constructs[this_construct].xsel = xsel

                else:
                    if line.strip():
                        raise AttributeError('Invalid section: ' + line)

            else:
                raise ValueError('Wrong version number %s' % version)

            if self.version == 0.1 and fill_data_types:
                self.data_types[this_construct] = [self.data_types[this_construct][0]] * len(self.values[this_construct])

        if self.version >= 0.2:
            self.comments = self.comments[:-1]
        self.loaded = True


    def save_construct(self, construct, data, sequence, structure, offset, annotations, data_annotations, filename, comments='', version=0.32, seqpos=None, errors=[]):
        self.version = version
        self.constructs[construct] = RDATSection()
        self.constructs[construct].sequences = defaultdict(int)
        self.constructs[construct].structures = defaultdict(int)
        self.constructs[construct].name = construct
        self.comments = comments

        self.mutpos[construct] = []
        if seqpos is not None:
            self.constructs[construct].seqpos = seqpos
        else:
            self.constructs[construct].seqpos = [i + offset for i in range(len(sequence))]
        self.constructs[construct].sequence = sequence
        self.constructs[construct].structure = structure
        self.constructs[construct].offset = offset
        self.constructs[construct].annotations = {}
        self.constructs[construct].xsel = []

        self.constructs[construct].data = []
        self.constructs[construct].errors = []
        if isinstance(data_annotations, dict):
            self.values[construct] = [data]
            self.errors[construct] = [errors]
            self._append_new_data_section(construct)
            self.constructs[construct].data[0].values = data
            self.constructs[construct].data[0].annotations = data_annotations
            self.constructs[construct].data[0].errors = errors
        else:
            for i, data_annotation in enumerate(data_annotations):
                self.values[construct] = data
                if errors is not None: self.errors[construct] = errors
                self._append_new_data_section(construct)
                self.constructs[construct].data[-1].values = data[i, :]
                self.constructs[construct].data[-1].annotations = data_annotation
                if len(errors) > 0 : self.constructs[construct].data[-1].errors = errors[i, :]

        self.loaded = True
        self.save(filename)


    def save(self, filename, version=None, delim='\t'):
        if not version:
            version = 0.3 if not self.version else self.version

        if not self.loaded:
            raise UnboundLocalError('Data not loaded yet ...')
        else:
            f = open(filename, 'w')

            if version == 0.1:
                f.write('VERSION: %s\n' % str(self.version))
                f.write('COMMENTS: %s\n' % str(self.comments))
                f.write('ANNOTATION: %s\n' % self._annotation_str(self.annotations, delim))

                for name in self.constructs:
                    construct = self.constructs[name]
                    f.write('CONSTRUCT\n')
                    f.write('NAME: %s\n' % name)
                    f.write('SEQUENCE: %s\n' % construct.sequence)
                    f.write('WELLS: %s\n' % ','.join([x for x in construct.wells]))
                    f.write('OFFSET: %s\n' % str(construct.offset))
                    f.write('ANNOTATION: %s\n' % self._annotation_str(construct.annotations, delim))
                    for d in construct.data:
                        f.write('DATA\n')
                        f.write('SEQPOS: %s\n' % ','.join([str(x) for x in d.seqpos]))
                        f.write('ANNOTATION: %s\n' % self._annotation_str(d.annotations, delim))
                        f.write('VALUES: %s\n' % ','.join([str(x) for x in d.values]))

            elif version == 0.2 or version == 0.21:
                f.write('VERSION %s\n' % str(version))
                f.write('COMMENTS %s\n' % str(self.comments))

                for name in self.constructs:
                    construct = self.constructs[name]
                    f.write('CONSTRUCT\n')
                    f.write('NAME %s\n' % name)
                    f.write('SEQUENCE %s\n' % construct.sequence)
                    f.write('STRUCTURE %s\n' % construct.structure)
                    f.write('OFFSET %s\n' % str(construct.offset))
                    if construct.annotations:
                        f.write('ANNOTATION %s\n' % self._annotation_str(construct.annotations, delim))
                    f.write('MUTPOS %s\n' % delim.join([str(x) for x in self.mutpos[name]]))
                    f.write('SEQPOS %s\n' % delim.join([str(x) for x in construct.seqpos]))
                    f.write('DATA_TYPE %s\n' % ' '.join(self.data_types[name]))

                    for i, d in enumerate(construct.data):
                        f.write('ANNOTATION_DATA %s %s\n' % (i + 1, self._annotation_str(d.annotations, delim)))
                    for i, row in enumerate(self.values[name]):
                        f.write('AREA_PEAK %s %s\n' % (i + 1, delim.join([str(x) for x in row])))
                    for i, row in enumerate(self.traces[name]):
                        f.write('TRACE %s %s\n' % (i + 1, delim.join([str(x) for x in row])))

                    if self.errors:
                        for i, row in enumerate(self.errors[name]):
                            f.write('AREA_PEAK_ERRORS %s %s\n' % (i + 1, delim.join([str(x) for x in row])))
                    if construct.xsel:
                        f.write('XSEL %s\n' % ' '.join([str(x) for x in construct.xsel]))

                    for i, row in enumerate(self.xsels[name]):
                        f.write('XSEL_REFINE %s %s\n' % (i + 1, delim.join([str(x) for x in row])))

            elif version >= 0.24 and version < 0.4:
                f.write('RDAT_VERSION%s%s\n' % (delim, str(version)))

                if len(self.constructs) == 1:
                    name = list(self.constructs.keys())[0]
                    construct = self.constructs[name]
                    f.write('NAME%s%s\n' % (delim, name))
                    f.write('SEQUENCE%s%s\n' % (delim, construct.sequence))
                    # if 'sequences' in construct.__dict__:
                    #     for k, v in construct.sequences.iteritems():
                    #         f.write('SEQUENCE:%s%s%s\n' % (k, delim, v))
                    f.write('STRUCTURE%s%s\n' % (delim, construct.structure))
                    # if 'structures' in construct.__dict__:
                    #     for k, v in construct.structures.iteritems():
                    #         f.write('STRUCTURE:%s%s%s\n' % (k, delim, v))
                    f.write('OFFSET%s%s\n\n' % (delim, str(construct.offset)))

                    if self.comments:
                        for com in self.comments.split('\n'):
                            f.write('COMMENT%s%s\n' % (delim, com))
                        f.write('\n')

                    if version < 0.32:
                        if name in self.mutpos:
                            f.write('MUTPOS%s%s\n' % (delim, delim.join([str(x) for x in self.mutpos[name]])))
                        else:
                            f.write('MUTPOS%s%s\n' % (delim, 'WT ' * len(construct.data)))
                    if version >= 0.32:
                        f.write('SEQPOS%s%s\n\n' % (delim, delim.join([construct.sequence[x - construct.offset - 1] + str(x) for x in construct.seqpos])))
                    else:
                        f.write('SEQPOS%s%s\n\n' % (delim, delim.join([str(x + 1) for x in construct.seqpos])))

                    f.write('ANNOTATION%s%s\n\n' % (delim, self._annotation_str(self.annotations, delim)))
                    if construct.annotations:
                        f.write('ANNOTATION%s%s\n' % (delim, self._annotation_str(construct.annotations, delim)))

                    for i, d in enumerate(construct.data):
                        f.write('ANNOTATION_DATA:%s%s%s\n' % (i + 1, delim, self._annotation_str(d.annotations, delim)))
                    f.write('\n')

                    for i, row in enumerate(self.values[name]):
                        f.write('REACTIVITY:%s%s%s\n' % (i + 1, delim, delim.join(['%.4f' % x for x in row])))
                    for i, row in enumerate(self.traces[name]):
                        if len(row):
                            f.write('TRACE:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))
                    for i, row in enumerate(self.reads[name]):
                        if len(row):
                            f.write('READS:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))
                    for i, row in enumerate(self.errors[name]):
                        if len(row):
                            f.write('REACTIVITY_ERROR:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))
                    if construct.xsel:
                        f.write('XSEL%s%s\n' % (delim, delim.join([str(x) for x in construct.xsel])))
                    for i, row in enumerate(self.xsels[name]):
                        if len(row):
                            f.write('XSEL_REFINE:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))

                else:
                    if self.comments:
                        for com in self.comments.split('\n'):
                            f.write('COMMENT%s%s\n' % (delim, com))
                        f.write('\n')

                    f.write('ANNOTATION%s%s\n\n' % (delim, self._annotation_str(self.annotations, delim)))

                    for name in self.constructs:
                        construct = self.constructs[name]
                        f.write('NAME%s%s\n' % (delim, name))

                        if construct.sequence:
                            f.write('SEQUENCE%s%s\n' % (delim, construct.sequence))
                        # if 'sequences' in construct.__dict__:
                        #     for k, v in construct.sequences.iteritems():
                        #         f.write('SEQUENCE:%s%s%s\n' % (k, delim, v))
                        if construct.structure:
                            f.write('STRUCTURE%s%s\n' % (delim, construct.structure))
                        # if 'structures' in construct.__dict__:
                        #     for k, v in construct.structures.iteritems():
                        #         f.write('STRUCTURE:%s%s%s\n' % (k, delim, v))
                        f.write('OFFSET%s%s\n\n' % (delim, str(construct.offset)))

                        if construct.annotations:
                            f.write('ANNOTATION%s%s\n' % (delim, self._annotation_str(construct.annotations, delim)))
                        if version < 0.32:
                            if name in self.mutpos:
                                f.write('MUTPOS%s%s\n' % (delim, delim.join([str(x) for x in self.mutpos[name]])))
                            else:
                                f.write('MUTPOS%s%s\n' % (delim, 'WT ' * len(construct.data)))
                        if version >= 0.32:
                            f.write('SEQPOS%s%s\n\n' % (delim, delim.join([construct.sequence[x - construct.offset - 1] + str(x) for x in construct.seqpos])))
                        else:
                            f.write('SEQPOS%s%s\n\n' % (delim, delim.join([str(x + 1) for x in construct.seqpos])))

                        for i, d in enumerate(construct.data):
                            f.write('ANNOTATION_DATA:%s%s%s\n' % (i + 1, delim, self._annotation_str(d.annotations, delim)))
                        f.write('\n')

                        if name in self.values:
                            for i, row in enumerate(self.values[name]):
                                f.write('REACTIVITY:%s%s%s\n' % (i + 1, delim, delim.join(['%.4f' % x for x in row])))
                        if name in self.traces:
                            for i, row in enumerate(self.traces[name]):
                                if len(row):
                                    f.write('TRACE:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))
                        if name in self.reads:
                            for i, row in enumerate(self.reads[name]):
                                if len(row):
                                    f.write('READS:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))
                        if name in self.errors:
                            for i, row in enumerate(self.errors[name]):
                                if len(row):
                                    f.write('REACTIVITY_ERRORS:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))
                        if construct.xsel:
                            f.write('XSEL%s%s\n' % (delim, delim.join([str(x) for x in construct.xsel])))
                        if name in self.xsels:
                            for i, row in enumerate(self.xsels[name]):
                                if len(row):
                                    f.write('XSEL_REFINE:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))

            elif version >= 0.4:
                f.write('RDATVERSION%s%s\n' % (delim, str(version)))

                if len(self.constructs) == 1:
                    name = list(self.constructs.keys())[0]
                    construct = self.constructs[name]
                    f.write('NAME%s%s\n' % (delim, name))

                    if construct.sequence:
                        f.write('SEQUENCE%s%s\n' % (delim, construct.sequence))
                    if construct.structure:
                        f.write('STRUCTURE%s%s\n' % (delim, construct.structure))
                    f.write('OFFSET%s%s\n' % (delim, str(construct.offset)))
                    f.write('SEQPOS%s%s\n\n' % (delim, delim.join([construct.sequence[x - construct.offset - 1] + str(x) for x in construct.seqpos])))

                    if self.comments:
                        for com in self.comments.split('\n'):
                            f.write('COMMENT%s%s\n' % (delim, com))
                        f.write('\n')

                    f.write('ANNOTATION%s%s\n\n' % (delim, self._annotation_str(self.annotations, delim)))
                    for i, d in enumerate(construct.data):
                        f.write('DATA_ANNOTATION:%s%s%s\n' % (i + 1, delim, self._annotation_str(d.annotations, delim)))
                    f.write('\n')

                    if name in self.values:
                        for i, row in enumerate(self.values[name]):
                            f.write('DATA:%s%s%s\n' % (i + 1, delim, delim.join(['%.4f' % x for x in row])))

                    if construct.xsel:
                        f.write('\nXSEL%s%s\n' % (delim, delim.join([str(x) for x in construct.xsel])))
                    for i, row in enumerate(self.xsels[name]):
                        if len(row) > 0:
                            f.write('XSEL_REFINE:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))

                else:
                    if self.comments:
                        for com in self.comments.split('\n'):
                            f.write('COMMENT%s%s\n' % (delim, com))
                        f.write('\n')

                    f.write('ANNOTATION%s%s\n\n' % (delim, self._annotation_str(self.annotations, delim)))

                    for name in self.constructs:
                        construct = self.constructs[name]
                        f.write('NAME%s%s\n' % (delim, name))

                        if construct.sequence:
                            f.write('SEQUENCE%s%s\n' % (delim, construct.sequence))
                        if construct.structure:
                            f.write('STRUCTURE%s%s\n' % (delim, construct.structure))
                        f.write('OFFSET%s%s\n\n' % (delim, str(construct.offset)))

                        if construct.annotations:
                            f.write('ANNOTATION%s%s\n' % (delim, self._annotation_str(construct.annotations, delim)))
                        f.write('SEQPOS%s%s\n\n' % (delim, delim.join([construct.sequence[x - construct.offset - 1] + str(x) for x in construct.seqpos])))

                        for i, d in enumerate(construct.data):
                            f.write('DATA_ANNOTATION:%s%s%s\n' % (i + 1, delim, self._annotation_str(d.annotations, delim)))
                        f.write('\n')

                        if name in self.values:
                            for i, row in enumerate(self.values[name]):
                                f.write('DATA:%s%s%s\n' % (i + 1, delim, delim.join(['%.4f' % x for x in row])))

                        if construct.xsel:
                            f.write('\nXSEL%s%s\n' % (delim, delim.join([str(x) for x in construct.xsel])))
                        if name in self.xsels:
                            for i, row in enumerate(self.xsels[name]):
                                if len(row) > 0:
                                    f.write('XSEL_REFINE:%s%s%s\n' % (i + 1, delim, delim.join([str(x) for x in row])))

            else:
                raise ValueError('Wrong version number %s' % version)
            f.close()


    def validate(self):
        messages = []
        for name in self.constructs:
            c = self.constructs[name]
            if len(name) == 0:
                messages.append('WARNING! Must give a name!')
            if len(c.sequence) == 0:
                messages.append('WARNING! Must supply sequence!')
            if 'T' in c.sequence:
                messages.append('WARNING! Warning: you have a T instead of a U in the sequence!!')
            if min(c.seqpos) - c.offset < 1:
                messages.append('WARNING! Offset/seqpos does not look right -- at least one index is too low for sequence')
            if max(c.seqpos) - c.offset > len(c.sequence):
                messages.append('WARNING! Offset/seqpos does not look right -- at least one index is too high for sequence')
            if len(c.data[0].values) != len(c.seqpos):
                messages.append('WARNING! Number of bands in area_peak [%s] does not match len of seqpos [%s]' % (len(c.data[0].values), len(c.seqpos)))

            for i, d in enumerate(c.data):
                if not hasattr(d, 'annotations'):
                    messages.append('WARNING! Data for index %s has no annotations' % i)
                if not hasattr(d, 'values'):
                    messages.append('WARNING! Data for index %s has no values for area peaks' % i)
                if not hasattr(d, 'trace'):
                    messages.append('WARNING! Data for index %s has no trace' % i)
                if len(self.xsels) > 0 and (not hasattr(d, 'xsel')):
                    messages.append('WARNING! Data for index %s has no xsel refine' % i)
                if hasattr(d, 'xsel'):
                    if len(c.xsel) != 0 and len(c.xsel) != len(d.values):
                        messages.append('WARNING! Number of bands in construct xsel [%s] does not match number of bands in values area peak [%s] of data indexed %s' % (len(c.xsel), len(d.values), i))
                    if len(d.xsel) != 0 and len(d.xsel) != len(d.values):
                        messages.append('WARNING! Number of bands in xsel indexed %s [%s] does not match number of bands in values area peak [%s]' % (i, len(d.xsel), len(d.values)))
            return messages


    def toISATAB(self):

        def parse_concentration(s):
            for i, ch in enumerate(s):
                if not ch in [str(x) for x in range(10)] + ['.']:
                    idx = i
                    break
            return s[:idx], s[idx:]

        isatabfile = ISATABFile()
        j = 0
        (general_factors, general_protocol) = ([], '')
        (protocols, chemicals) = (set(), set())
        isatabfile.investigation_dict['Study File Name'].append(self.filename + ' (check entry ' + self.filename[:self.filename.find('.')] + ' at https://rmdb.stanford.edu/deposit/validate/ for details)')
        """
        if 'pmid' in self.annotations:
            pmid = self.annotations['pmid'][0]
            h = Entrez.efetch(db='pubmed', id=[pmid], rettype='medline', retmode='text')
            records = Medline.parse(h)
            for r in records:
                record = r
            doi = ''
            for item in record.get('AID', '?'):
                if '[doi]' in item:
                    doi = item.replace('[doi]','')
            isatabfile.investigation_dict['Study Title'].append(record.get('TI', '?'))
            isatabfile.investigation_dict['Study Publication DOI'].append(doi)
            isatabfile.investigation_dict['Study Publication Title'].append(record.get('TI', '?'))
            isatabfile.investigation_dict['Study Publication Author list'].append(','.join(record.get('AU', '?')))
            isatabfile.investigation_dict['Study Description'].append(record.get('AB', '?'))
            isatabfile.investigation_dict['Study Public Release Date'].append(record.get('DP', '?'))
            isatabfile.investigation_dict['Study PubMed ID'].append(pmid)
        isatabfile.investigation_dict['Study Publication Status'].append('indexed in pubmed')
        """

        for k in ['chemical', 'salt', 'buffer', 'temperature']:
            if k in self.annotations:
                isatabfile.investigation_dict['Study Factor Name'].append('%s %s (constant for all assays)' % (k, self.annotations[k]))
                if k != 'temperature':
                    isatabfile.investigation_dict['Study Factor Type'].append('Compound')
                else:
                    isatabfile.investigation_dict['Study Factor Type'].append('Temperature')

        if 'technology' in self.annotations:
            tech = self.annotations['technology'][0]
        else:
            tech = 'capillary electrophoresis'
        if 'modifier' in self.annotations:
            if self.annotations['modifier'][0] in ONTOLOGY.MODIFIER_PROTOCOL:
                general_protocol = ONTOLOGY.MODIFIER_PROTOCOL[self.annotations['modifier'][0]]
            else:
                general_protocol = self.annotations['modifier'][0]
            protocols.add(general_protocol)

        """
        for k in ['chemical', 'salt', 'buffer']:
            if k in self.annotations:
                name, concentration = self.annotations[k][0].split(':')
                isatabfile.investigation_dict['Study Factor Name'].append(name)
                isatabfile.investigation_dict['Study Factor Name'].append(name + ' concentration')
                term = ONTOLOGY.CHEMICAL[name]
                concentration, units = parse_concentration(concentration)
                general_factors.append(['Factor Value[%s]' % name, name])
                general_factors.append(['Term Source REF[%s]' % name, term.split(':')[0]])
                general_factors.append(['Term Accession Number[%s]' % name, term])
                general_factors.append(['Factor Value[%s concentration]' % name, name])
                general_factors.append(['Unit[%s concentration]' % name, units])
                isatabfile.assays_keys.append('Factor Value[%s]' % name)
                isatabfile.assays_keys.append('Term Source REF[%s]' % name)
                isatabfile.assays_keys.append('Term Accession Number[%s]' % name)
                isatabfile.assays_keys.append('Factor Value[%s concentration]' % name)
                isatabfile.assays_keys.append('Unit[%s concentration]' % name)
        """
        for cname in self.constructs:
            construct = self.constructs[cname]
            for i, d in enumerate(construct.data):
                name = cname.replace(' ', '-')
                seq = ''
                for j in construct.seqpos:
                    seq += construct.sequence[j - construct.offset - 1]
                if 'mutation' in d.annotations:
                    for j in range(len(d.annotations['mutation'])):
                        mutlabel = d.annotations['mutation'][j].strip('\t')
                        name = name + '_' + mutlabel
                        if mutlabel != 'WT':
                            try:
                                index = int(mutlabel[1:-1]) - construct.offset
                                seq = seq[:index-1] + mutlabel[-1] + seq[index:]
                            except ValueError:
                                # The mutation label is not in standard format, default to normal sequence and make a note
                                seq += ' Note, mutation=%s' % mutlabel
                else:
                    name = name + '_WT'

                idname = name + '_' + str(i + 1)
                isatabfile.assays_dict['Assay Name'].append(idname)
                isatabfile.sample_id_name_map[idname] = name
                isatabfile.data[idname] = d.values
                isatabfile.data_id_order.append(idname)
                isatabfile.assays_dict['Source Name'].append(name)
                isatabfile.assays_dict['Characteristics[Nucleotide Sequence]'].append(seq)
                isatabfile.assays_dict['Characteristics[Nucleotide Type]'].append('RNA')

                if 'production' in d.annotations:
                    isatabfile.assays_dict['Characteristics[RNA Production]'].append(d.annotations['production'][0].replace('-', ' '))
                else:
                    isatabfile.assays_dict['Characteristics[RNA Production]'].append('in vitro synthesis')

                if 'modifier' in d.annotations:
                    modifier = d.annotations['modifier'][0]
                    if modifier in ONTOLOGY.MODIFIER_PROTOCOL:
                        isatabfile.assays_dict['Protocol REF'].append(ONTOLOGY.MODIFIER_PROTOCOL[d.annotations['modifier'][0]].replace('-', ' '))
                        protocol = ONTOLOGY.MODIFIER_PROTOCOL[d.annotations['modifier'][0]]
                    else:
                        isatabfile.assays_dict['Protocol REF'].append(modifier)
                        protocol = modifier
                    protocols.add(protocol)
                elif 'modifier' in self.annotations:
                    modifier = self.annotations['modifier'][0]
                    if modifier in ONTOLOGY.MODIFIER_PROTOCOL:
                        isatabfile.assays_dict['Protocol REF'].append(ONTOLOGY.MODIFIER_PROTOCOL[self.annotations['modifier'][0]].replace('-', ' '))
                        protocol = ONTOLOGY.MODIFIER_PROTOCOL[self.annotations['modifier'][0]]
                    else:
                        isatabfile.assays_dict['Protocol REF'].append(modifier)
                        protocol = modifier
                    protocols.add(protocol)
                elif general_protocol:
                    isatabfile.assays_dict['Protocol REF'].append(general_protocol)
                    protocol = general_protocol
                else:
                    isatabfile.assays_dict['Protocol REF'].append('')
                    protocol = ''

                isatabfile.assays_dict['Parameter Value[Data Starts at Sequence Position]'].append(str(min(construct.seqpos) - construct.offset))
                isatabfile.assays_dict['Raw Data File'].append('datamatrix.txt')
                if 'performer' in d.annotations:
                    isatabfile.assays_dict['Performer'].append(d.annotations['performer'][0].replace('-', ' '))
                elif 'performer' in self.annotations:
                    isatabfile.assays_dict['Performer'].append(self.annotations['performer'][0].replace('-', ' '))
                else:
                    isatabfile.assays_dict['Performer'].append('')

                if 'date' in d.annotations:
                    isatabfile.assays_dict['Date'].append(d.annotations['date'][0])
                elif 'date' in self.annotations:
                    isatabfile.assays_dict['Date'].append(self.annotations['date'][0])
                else:
                    isatabfile.assays_dict['Date'].append('')
                isatabfile.assays_dict['Term Source REF'].append('OBI')

                if protocol in ONTOLOGY.PROTOCOL:
                    isatabfile.assays_dict['Term Accession Number'].append(ONTOLOGY.PROTOCOL[protocol])
                else:
                    isatabfile.assays_dict['Term Accession Number'].append(protocol)

                for k in ['chemical', 'folding-salt', 'buffer', 'salt']:
                    if k in d.annotations:
                        chemical, concentration = d.annotations[k][0].split(':')
                        concentration, units = parse_concentration(concentration)
                        if chemical in ONTOLOGY.CHEMICAL:
                            term = ONTOLOGY.CHEMICAL[chemical]
                        else:
                            term = chemical
                        if not k in isatabfile.assays_factors:
                            isatabfile.assays_factors[k] = {}
                            isatabfile.assays_factors[k]['value'] = []
                            isatabfile.assays_factors[k]['ref'] = []
                            isatabfile.assays_factors[k]['concentration'] = []
                            isatabfile.assays_factors[k]['unit'] = []
                            isatabfile.assays_factors[k]['accession'] = []
                        isatabfile.assays_factors[k]['value'].append(chemical)
                        isatabfile.assays_factors[k]['ref'].append(term.split(':')[0])
                        isatabfile.assays_factors[k]['concentration'].append(concentration)
                        isatabfile.assays_factors[k]['unit'].append(units)
                        isatabfile.assays_factors[k]['accession'].append(term)
                        chemicals.add(k.replace('-', ' '))
                        chemicals.add(k.replace('-', ' ') + ' concentration')
                    else:
                        if k in chemicals:
                            isatabfile.assays_factors[k]['value'].append('')
                            isatabfile.assays_factors[k]['ref'].append('')
                            isatabfile.assays_factors[k]['concentration'].append('')
                            isatabfile.assays_factors[k]['unit'].append('')
                            isatabfile.assays_factors[k]['accession'].append('')

                for p in general_factors:
                    if p[0] in isatabfile.assays_dict:
                        isatabfile.assays_dict[p[0]].append(p[1])
                    else:
                        isatabfile.assays_dict[p[0]] = [p[1]]

            for p in protocols:
                if p in ONTOLOGY.PROTOCOL:
                    term = ONTOLOGY.PROTOCOL[p]
                else:
                    term = p
                isatabfile.investigation_dict['Study Protocol Name'].append(p.replace('-', ' '))
                isatabfile.investigation_dict['Study Protocol Type Term Source REF'].append('OBI')
                isatabfile.investigation_dict['Study Assay Measurement Type'].append(p.replace('-', ' '))
                isatabfile.investigation_dict['Study Assay Measurement Type Term Accession Number'].append(term)
                isatabfile.investigation_dict['Study Assay Measurement Type Term Source REF'].append(term.split(':')[0])
                isatabfile.investigation_dict['Study Assay Technology Type'].append(tech)
                isatabfile.investigation_dict['Study Assay File Name'].append('study-assay.txt')
            for ch in chemicals:
                isatabfile.investigation_dict['Study Factor Name'].append(ch)
                isatabfile.investigation_dict['Study Factor Type'].append('Compound')
                isatabfile.investigation_dict['Study Factor Type Term Source REF'].append('CHEBI')

        return isatabfile



"""
For the ISATAB format for chemical footprinting experiments
"""

class ISATABFile(object):
    def __init__(self):
        self.assays_keys = copy.deepcopy(ISATAB_ASSAY_KEYS)
        self.investigation_dict = copy.deepcopy(ISATAB_INVEST_DICT)
        self.assays_dict = copy.deepcopy(ISATAB_ASSAY_DICT)

        self.sample_id_name_map = {}
        self.data_id_order = []
        self.assays_factors = {}
        self.data = {}


    def load(self, name, type='xls'):
        self.name = name

        if type == 'dir':
            investigation_file = open(name + '/investigation.txt')
            for l in investigation_file.readlines():
                fields = l.strip().split('\t')
                self.investigation_dict[fields[0]] = fields[1:]

            assays_file = open(name+'/'+self.investigation_dict['Study Assay File Name'][0])
            assays_keys = assays_file.readline().strip().split('\t')
            l = assays_file.readline()

            while l:
                for i, f in enumerate(l.strip().split('\t')):
                    if assays_keys[i] in self.assays_dict:
                        self.assays_dict[assays_keys[i]].append(f)
                    else:
                        self.assays_dict[assays_keys[i]] = [f]
                l = assays_file.readline()

            datamatrix_file = open(name + '/' + self.assays_dict['Raw Data File'][0])
            data_keys = datamatrix_file.readline().strip().split('\t')
            for i, k in enumerate(data_keys):
                self.data[k + '_' + str(i+1)] = []
                self.sample_id_name_map[k + '_' + str(i + 1)] = k
                self.data_id_order.append(k + '_' + str(i + 1))
            l = datamatrix_file.readline()

            while l:
                for i, f in enumerate(l.strip().split('\t')):
                    if f:
                        self.data[data_keys[i] + '_' + str(i + 1)].append(float(f))
                l = datamatrix_file.readline()

            assays_file.close()
            datamatrix_file.close()
            investigation_file.close()

        elif type == 'xls':
            wb = xlrd.open_workbook(name)
            investigation_sh = wb.sheet_by_name('investigation')

            for j in range(investigation_sh.nrows):
                fields = investigation_sh.row_values(j)
                self.investigation_dict[fields[0]] = fields[1:]

            try:
                assays_sh = wb.sheet_by_name(self.investigation_dict['Study Assay File Name'][0].replace('.txt', ''))
            except xlrd.biffh.XLRDError:
                assays_sh = wb.sheet_by_name(self.investigation_dict['Study File Name'][0].replace('.txt', ''))

            assays_keys = assays_sh.row_values(0)
            for j in range(1, assays_sh.nrows):
                l = assays_sh.row_values(j)
                for i, f in enumerate(l):
                    if assays_keys[i] in self.assays_dict:
                        self.assays_dict[assays_keys[i]].append(f)
                    else:
                        self.assays_dict[assays_keys[i]] = [f]

            datamatrix_sh = wb.sheet_by_name(self.assays_dict['Raw Data File'][0].replace('.txt', ''))
            data_keys = datamatrix_sh.row_values(0)
            for i, k in enumerate(data_keys):
                self.data[k + '_' + str(i + 1)] = []
                self.sample_id_name_map[k + '_' + str(i + 1)] = k
                self.data_id_order.append(k + '_' + str(i + 1))

            for j in range(1, datamatrix_sh.nrows):
                l = datamatrix_sh.row_values(j)
                for i, f in enumerate(l):
                    if f:
                        self.data[data_keys[i] + '_' + str(i + 1)].append(float(f))
        else:
            raise TypeError('Unrecognized type %s for loading isatab file' % type)


    def save(self, name, type='xls'):
        self.name = name

        if type == 'dir':
            if not os.path.exists(name):
                os.mkdir(name)
            investigation_file = open(name + '/investigation.txt', 'w')
            assays_file = open(name + '/study-assay.txt', 'w')
            datamatrix_file = open(name + '/datamatrix.txt', 'w')

            for k in self.assays_keys:
                assays_file.write(k + '\t')
            for k in self.assays_factors:
                assays_file.write('Factor Value[%s]\t' % k.replace('-', ' '))
                assays_file.write('Term Source REF\t')
                assays_file.write('Term Accession Number\t')
                assays_file.write('Factor Value[%s concentration]\t' % k.replace('-', ' '))
                assays_file.write('Unit\t')
            assays_file.write('\n')

            for i in range(len(list(self.assays_dict.values())[0])):
                line = ''
                for k in self.assays_keys:
                    if k in self.assays_dict:
                        if len(self.assays_dict[k]) <= i:
                            line += '\t'
                        else:
                            line += self.assays_dict[k][i] + '\t'
                for k in self.assays_factors:
                    line += self.assays_factors[k]['value'][i] + '\t'
                    line += self.assays_factors[k]['ref'][i] + '\t'
                    line += self.assays_factors[k]['accession'][i] + '\t'
                    line += self.assays_factors[k]['concentration'][i] + '\t'
                    line += self.assays_factors[k]['unit'][i] + '\t'
                assays_file.write(line.strip('\t') + '\n')

            for k in ISATAB_INVEST_KEYS:
                line = k + '\t'
                for i in range(len(self.investigation_dict[k])):
                    line += self.investigation_dict[k][i] + '\t'
                investigation_file.write(line.strip('\t') + '\n')

            maxlen = max([len(x) for x in list(self.data.values())])
            datamatrix_file.write('\t'.join(self.assays_dict['Source Name']) + '\n')
            for i in range(maxlen):
                datamatrix_file.write('\t'.join(['' if i >= len(self.data[k]) else str(self.data[k][i]) for k in self.data_id_order]) + '\n')

            assays_file.close()
            datamatrix_file.close()
            investigation_file.close()

        elif type == 'xls':
            (assay_row, inv_row, data_row) = (0, 0, 0)

            wb = xlwt.Workbook()
            investigation_sh = wb.add_sheet('investigation')
            assays_sh = wb.add_sheet('study-assay')
            datamatrix_sh = wb.add_sheet('datamatrix')

            for i, k in enumerate(self.assays_keys):
                assays_sh.write(assay_row, i, k)
            for k in self.assays_factors:
                i += 1
                assays_sh.write(assay_row, i, 'Factor Value[%s]\t' % k.replace('-', ' '))
                i += 1
                assays_sh.write(assay_row, i, 'Term Source REF\t')
                i += 1
                assays_sh.write(assay_row, i, 'Term Accession Number\t')
                i += 1
                assays_sh.write(assay_row, i, 'Factor Value[%s concentration]\t' % k.replace('-', ' '))
                i += 1
                assays_sh.write(assay_row, i, 'Unit\t')

            assay_row += 1
            for i in range(len(list(self.assays_dict.values())[0])):
                line = []
                for k in self.assays_keys:
                    if k in self.assays_dict:
                        if len(self.assays_dict[k]) <= i:
                            line.append('')
                        else:
                            line.append(self.assays_dict[k][i])
                if i < len(self.assays_factors):
                    for k in self.assays_factors:
                        line .append(self.assays_factors[k]['value'][i])
                        line .append(self.assays_factors[k]['ref'][i])
                        line .append(self.assays_factors[k]['accession'][i])
                        line .append(self.assays_factors[k]['concentration'][i])
                        line .append(self.assays_factors[k]['unit'][i])
                for j in range(len(line)):
                    assays_sh.write(assay_row, j, line[j])
                assay_row += 1

            for i, k in enumerate(ISATAB_INVEST_KEYS):
                line = [k]
                for j in range(len(self.investigation_dict[k])):
                    line.append(self.investigation_dict[k][j])
                for j in range(len(line)):
                    investigation_sh.write(inv_row, j, line[j])
                inv_row += 1

            maxlen = max([len(x) for x in list(self.data.values())])
            for i, k in enumerate(self.assays_dict['Source Name']):
                datamatrix_sh.write(data_row, i, k)
            data_row += 1

            for i in range(maxlen):
                for j, k in enumerate(['' if i >= len(self.data[k]) else str(self.data[k][i]) for k in self.data_id_order]):
                    datamatrix_sh.write(data_row, j, k)
                data_row += 1
            wb.save(name)

        else:
            raise TypeError('Unrecognized type %s for saving isatab file' % type)


    def validate(self):
        messages = []

        def check_terms(d, prefix, ontdict):
            m = []
            for i, t in enumerate(d[prefix]):
                if not t.replace(' ', '-') in ontdict:
                    if t.strip() == '':
                        pass
                    else:
                        m.append('WARNING! For %s, term %s is unknown for its respective ontology' % (prefix, t))
                else:
                    r = d[prefix + ' Term Accession Number'][i].strip().replace('_', ':')
                    if ontdict[t.replace(' ', '-')] != r:
                        m.append('WARNING! For %s, term %s and accession number %s do not match' % (prefix, t, r))

            for i, t in enumerate(d[prefix + ' Term Accession Number']):
                if not d[prefix + ' Term Source REF'][i] in t:
                    m.append('WARNING! For %s, accession number and REF ontology do not match.' % prefix)
            return m

        if len(self.data) != len(self.assays_dict['Source Name']):
            messages.append('WARNING! Number of samples in assays and data file do not match')
        for k in self.data:
            if not self.sample_id_name_map[k] in self.assays_dict['Source Name']:
                messages.append('WARNING! Sample %s in data file not referenced from assays file' % k)
        for i, k in enumerate(self.assays_dict['Source Name']):
            if not k + '_' + str(i) in self.data:
                messages.append('WARNING! Sample %s in assays file is missing from data file' % k)
        for p in self.investigation_dict['Study Protocol Name']:
            if p.strip() != '' and (p.strip().lower() not in [x.strip().lower() for x in self.assays_dict['Protocol REF']]):
                messages.append('WARNING! Protocol %s in investigation file is not in assays file' % p)
        # Checks for ontology term consistency
        #messages = messages + check_terms(self.investigation_dict, 'Study Assay Measurement Type', ONTOLOGY.PROTOCOL)
        #messages = messages + check_terms(self.investigation_dict, 'Study Factor Type', ONTOLOGY.CHEMICAL)
        for i, seq in enumerate(self.assays_dict['Characteristics[Nucleotide Sequence]']):
            if i >= len(self.assays_dict['Source Name']):
                messages.append('ERROR! Cannot continue validation as list of source names and sequences do not match in number!')
                return messages
            for k in ['Source Name', 'Parameter Value[Data Starts at Sequence Position]', 'Characteristics[Nucleotide Type]']:
                stop = False
                if len(self.assays_dict[k]) == 0:
                    messages.append('ERROR! Cannot continue validation as "%s" column is not present in assays tab! (check any spelling inconsistencies in column names)' % k)
                    stop = True
                if stop:
                    return messages

            sn = self.assays_dict['Source Name'][i]
            idn = sn + '_' + str(i)
            m = int(self.assays_dict['Parameter Value[Data Starts at Sequence Position]'][i])
            if idn in self.data and len(self.data[idn]) != len(seq) - m + 1:
                messages.append('WARNING! Length of data lane %s [%s] and length of respective sequence [%s] do not match' % (sn, len(self.data[idn]), len(seq) - m + 1))

            chartype = self.assays_dict['Characteristics[Nucleotide Type]'][i]
            if chartype == 'RNA' and 'T' in seq:
                messages.append('WARNING! Sequence for %s specified as RNA but looks like DNA.' % n)
            if chartype == 'DNA' and 'U' in seq:
                messages.append('WARNING! Sequence for %s specified as DNA but looks like RNA.' % n)
        return messages


    def toRDAT(self):
        rdatfile = RDATFile()
        general_annotations = defaultdict(list)

        for i, name in enumerate(self.assays_dict['Assay Name']):
            d = RDATSection(['seqpos', 'errors', 'trace', 'xsel'])
            d.annotations = {}
            d.values = self.data[name]
            rdatfile.values[name] = [d.values]

            c = RDATSection(['values', 'traces', 'mutpos', 'data_types', 'xsel'])
            c.name = name
            c.annotations = {}
            for k, v in general_annotations.items():
                c._annotation_strtations[k] = v

            c.sequence = self.assays_dict['Characteristics[Nucleotide Sequence]'][i]
            c.seqpos = list(range(len(c.sequence)))
            c.data = [d]
            c.offset = 0
            c.structure = '.' * len(c.sequence)
            rdatfile.constructs[name] = c

        rdatfile.loaded = True
        return rdatfile

