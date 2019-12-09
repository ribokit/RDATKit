import subprocess

if __package__ is None or not __package__:
    from .util import *
    from .secstr import SecondaryStructure
else:
    from .util import *
    from .secstr import SecondaryStructure


class VARNA(object):

    def __init__(self, sequences=[], structures=[], mapping_data=[]):
        self.sequences = sequences
        self.structures = structures
        self.mapping_data = mapping_data

        (self.rows, self.columns, self.width, self.height) = (1, 1, 522, 200)
        (self.annotation_font_size, self.annotation_color) = (9, '#FF0000')


    @classmethod
    def get_colormap(self, values, sequence=''):
        is_incomplete = (len(sequence) > len(values))
        is_negative = any([x < 0 for x in values])
        if is_incomplete or is_negative:
            return '-0.001:#C0C0C0,0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000'
        else:
            return '0:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000'

    @classmethod
    def run(self, sequence, structure, outfile, options={}):
        option_str = ''
        for key in options:
            val = options[key]
            if isinstance(val, list):
                argval = str(val).strip('[]').replace('L', '').replace('u', '')
            else:
                argval = str(val)
            option_str += '-%s "%s" ' % (key, argval)
        cmd = 'java -cp %s fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN "%s" %s -o %s' % (PATH_VARNA_JAR, sequence, structure, option_str, outfile)
        subprocess.check_call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


    def get_values(self, att):
        if att == 'structures':
            return [struct.dbn for struct in self.structures]
        return getattr(self, att)

    def get_frames(self, overlap_structures):
        res = -1
        for key, val in list(self.__dict__.items()):
            if overlap_structures and key == 'structures':
                continue
            if isinstance(val, list):
                res = max(res, len(val))
        return res


    def _get_option_string(self, options):
        option_str = ''
        for key in options:
            val = options[key]
            if isinstance(val, list):
                arg_val = str(val).strip('[]').replace('L', '').replace('u', '')
            else:
                arg_val = str(val)
            option_str += '-%s "%s" ' % (key, arg_val)
        return option_str


    def _get_base_annotation_string(self, base_annotations, annotation_by_helix=False, helix_function=(lambda x, y: x), is_applet=False, base_offset=0, stype='L', helix_side=0):
        base_annotation_string = ''
        for i, ba in enumerate(base_annotations):
            idx = str(i + 1) if len(base_annotations) > 1 else ''
            if is_applet:
                base_annotation_string += '<param name="annotations%s" value=\n\t"' % idx
            if annotation_by_helix:
                for helix in self.structures[i].helices():
                    anchor = helix[0][helix_side] + 1 + (helix[-1][helix_side] - helix[0][helix_side])/2 + base_offset
                    if helix[0] in ba:
                        annotation_value = ba[helix[0]]
                    elif helix[0][::-1] in ba:
                        annotation_value = ba[helix[0][::-1]]
                    else:
                        annotation_value = 0

                    for bp in helix:
                        if bp in ba:
                            nextval = ba[bp]
                            annotation_value = helix_function(annotation_value, nextval)
                    base_annotation_string += '%s:type=%s,anchor=%d,size=%d,color=%s;' % (annotation_value, stype, anchor, self.annotation_font_size, self.annotation_color)
            else:
                for bp in base_annotations:
                    base_annotation_string += '%s:type=B,anchor=%d,size=%d,color=%s;' % (base_annotations[bp], bp[0] + base_offset, self.annotation_font_size, self.annotation_color)
            base_annotation_string = base_annotation_string.strip()
        if is_applet:
            base_annotation_string += '"/>\n'
        return base_annotation_string


    def render(self, base_annotations=[], annotation_by_helix=False, helix_function=(lambda x, y: x), map_data_function=(lambda x: x), overlap_structures=False, reference_structure=SecondaryStructure(), annotation_def_val='', output='applet', cmd_options={}):
        struct_string = ''
        if output == 'applet':
            applet_string = '<applet code="VARNA.class" codebase="https://rmdb.stanford.edu/site_src/"\n'
            applet_string += 'archive="VARNA.jar" width="%d" height="%d">\n' % (self.width, self.height)
            frames = self.get_frames(overlap_structures=overlap_structures)

            if overlap_structures:
                struct_string += '<param name="structureDBN" value="%s"/>\n' % self.structures[0].dbn
                bps = self.structures[0].base_pairs()
                struct_string += '<param name="auxBPs" value="'

                for i in range(1, len(self.structures)):
                    for bp in self.structures[i].base_pairs():
                        if bp not in bps:
                            bps.append(bp)
                            struct_string += '(%s,%s):edge5=s,edge3=h,stericity=cis;' % (bp[0] + 1, bp[1] + 1)
                struct_string += '"/>\n'

            applet_string += self._get_base_annotation_string(base_annotations, annotation_by_helix=annotation_by_helix, helix_function=helix_function, is_applet=True)

            if self.rows * self.columns != frames:
                rows = frames
                columns = 1
            else:
                rows = self.rows
                columns = self.columns

            param_string = '<param name="rows" value="%d"/>\n' % rows
            param_string += '<param name="columns" value="%d"/>\n' % columns
            is_default_colormap = not hasattr(self, 'colorMapStyle')

            for att in self.__dict__:
                if overlap_structures and att == 'structures':
                    param_string += struct_string
                if isinstance(self.get_values(att), list):
                    if att == 'structures':
                        if overlap_structures:
                            continue
                        name = 'structureDBN'

                        if reference_structure:
                            for i in range(len(self.structures)):
                                frame_idx = str(i + 1) if len(self.structures) > 1 else ''
                                param_string += '<param name="%s%s" value="%s" />\n' % (name, frame_idx, self.structures[i].dbn)
                                auxbps_string = '<param name="auxBPs%s" value="' % frame_idx
                                refbps = reference_structure.base_pairs()
                                bps = self.structures[i].base_pairs()

                                for bp in bps:
                                    if bp not in refbps:
                                        auxbps_string += '(%s,%s):color=#00CCFF;' % (bp[0] + 1, bp[1] + 1)
                                for bp in refbps:
                                    if bp not in bps:
                                        auxbps_string += '(%s,%s):color=#FF6600;' % (bp[0] + 1, bp[1] + 1)
                                auxbps_string += '"/>\n'
                                param_string += auxbps_string
                            continue

                    elif att == 'sequences':
                        name = 'sequenceDBN'
                    elif att == 'mapping_data':
                        name = 'colorMap'
                    else:
                        name = att

                    for i, val in enumerate(self.get_values(att)):
                        frame_idx = str(i + 1) if frames > 1 else ''
                        if name == 'colorMap':
                            idx = i if len(self.sequences) == len(self.mapping_data) else 0
                            if is_default_colormap:
                                param_string += '<param name="colorMapStyle%s" value="%s" />\n' % (frame_idx, VARNA.get_colormap(val, sequence=self.sequences[idx]))
                            val = [map_data_function(val[x]) if (val[x] is not None) and val[x] > 0 else -0.001 for x in range(len(self.sequences[idx]))]
                            val = str(val).strip('[]').replace(' ', '')
                        param_string += '<param name="%s%s" value="%s" />\n' % (name, frame_idx, val)

            applet_string += param_string
            applet_string += '</applet>'
            return applet_string

        else:
            all_options = {}
            all_options.update(cmd_options)
            base_annotation_string = self._get_base_annotation_string(base_annotations, annotation_by_helix=annotation_by_helix, helix_function=helix_function)
            if 'annotations' in all_options:
                all_options['annotations'] = all_options['annotations'].strip('"') + base_annotation_string.strip('"')
            else:
                all_options['annotations'] = base_annotation_string
            option_str = self._get_option_string(all_options)

            i = 0
            cmd = ''
            for structure, sequence in zip(self.structures, self.sequences):
                if len(self.structures) > 1:
                    if i > 0:
                        cmd += '&& '
                    cmd += 'java -cp %s  fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN "%s" %s -o %s_%s' % (PATH_VARNA_JAR, sequence, structure.dbn, option_str, i, output)
                else:
                    cmd = 'java -cp %s  fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN "%s" %s -o %s' % (PATH_VARNA_JAR, sequence, structure.dbn, option_str, output)
                i += 1
            return cmd


