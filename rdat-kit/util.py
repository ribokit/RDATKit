# from scipy.stats import gamma
# if __package__ is None or not __package__:
#     from path import *
# else:
#     from .path import *


# PATH_RNA_STRUCTURE_FOLD = PATH_RNA_STRUCTURE + 'exe/Fold'
# PATH_RNA_STRUCTURE_ALLSUB = PATH_RNA_STRUCTURE + 'exe/AllSub'
# PATH_RNA_STRUCTURE_STOCHASTIC = PATH_RNA_STRUCTURE + 'exe/stochastic'
# PATH_RNA_STRUCTURE_ENERGY = PATH_RNA_STRUCTURE + 'exe/efn2'
# PATH_RNA_STRUCTURE_PARTITION = PATH_RNA_STRUCTURE + 'exe/partition'
# PATH_RNA_STRUCTURE_CT2DOT = PATH_RNA_STRUCTURE + 'exe/ct2dot'

# # these do not work. orphans by creator. no longer supported
# PATH_VIENNA_RNA_FOLD = PATH_VIENNA_RNA + 'RNAfold'
# PATH_VIENNA_RNA_ENERGY = PATH_VIENNA_RNA + 'RNA2Dfold'
# PATH_VIENNA_RNA_SUBOPT = PATH_VIENNA_RNA + 'RNAsubopt'


# DIST = {}
# for k in ['all', 'dangles', 'bulges', 'hairpins', 'interiorloops', 'helices', '3wayjunctions', '4wayjunctions', '5wayjunctions', 'unpaired', 'internalpairs', 'edgepairs']:
#     DIST[k] = gamma


class Ontology(object):
    def __init__(self):
        self.MODIFIER = {
            'DMS': 'CHEBI:59050',
            'DEPC': 'CHEBI:59051',
            'kethoxal': 'CHEBI:59052',
            'CMCT': 'CHEBI:59053',
            'NMIA': 'CHEBI:59054',
            'Fe(II)-BABE': 'CHEBI:59055',
            'methidiumpropyl-EDTA.FE(II)': 'CHEBI:59056',
            'ethylnitrosourea': 'CHEBI:23995',
            'lead': 'CHEBI:27889',
            'rhodium': 'CHEBI:33359',
            'ruthenium': 'CHEBI:30682',
            'terbium': 'CHEBI:33376',
            'DNAse-1': 'PRO:000006592',
            'RNAse-CL3': 'PRO:000025478',
            'nuclease-S1': 'PRO:000025471',
            'RNAse-T1': 'PRO:000025467',
            'RNAse-T2': 'PRO:000014060',
            'RNAse-U2': 'PRO:000025475',
            'RNAse-V1': 'PRO:000025477',
            'OH-radical': 'CHEBI:29191',
            'inline-probing': 'N/A',
            '1M7': 'CHEBI:60343',
            'RNAse-1': 'PRO:000014042'
        }

        self.CHEMICAL = {
            'MgCl2': 'CHEBI:6636',
            'Na-HEPES': 'NCIM:C0887237'
        }

        self.PROTOCOL = {
            'DMS-structure-mapping-assay': 'OBI:0001015',
            'DEPC-structure-mapping-assay': 'OBI:0000897',
            'kethoxal-structure-mapping-assay': 'OBI:0001013',
            'CMCT-structure-mapping-assay': 'OBI:0001006',
            'NMIA-RNA-structure-mapping-assay': 'OBI:0001026',
            'Fe-BABE-RNA-structure-mapping-assay': 'OBI:0001023',
            'MPE-Fe(II)-structure-mapping-assay': 'OBI:0001008',
            'ENU-structure-mapping-assay': 'OBI:0001011',
            'lead-structure-mapping-assay': 'OBI:0001020',
            'rhodium-DNA-structure-mapping-assay': 'OBI:0001018',
            'ruthenium-DNA-structure-mapping-assay': 'OBI:0001038',
            'terbium-RNA-structure-mapping-assay': 'OBI:0001027',
            'DNASE-1-structure-mapping-assay': 'OBI:0001016',
            'RNASE-CL3-structure-mapping-assay': 'OBI:0001005',
            'nuclease-S1-structure-mapping-assay': 'OBI:0001035',
            'RNASE-T1-structure-mapping-assay': 'OBI:0001030',
            'RNASE-T2-structure-mapping-assay': 'OBI:0001021',
            'RNASE-U2-structure-mapping-assay': 'OBI:0001024',
            'RNASE-V1-structure-mapping-assay': 'OBI:0001012',
            'OH-radical-structure-mapping-assay': 'OBI:0001029',
            'inline-Probing': 'OBI:0001039'
        }

        self.MODIFIER_PROTOCOL = {
            'DMS': 'DMS-structure-mapping-assay',
            'DEPC': 'DEPC-structure-mapping-assay',
            'kethoxal': 'kethoxal-structure-mapping-assay',
            'CMCT': 'CMCT-structure-mapping-assay',
            'NMIA': 'NMIA-RNA-structure-mapping-assay',
            'SHAPE': 'NMIA-RNA-structure-mapping-assay',
            'Fe(II)-BABE': 'Fe-BABE-RNA-structure-mapping-assay',
            'methidiumpropyl-EDTA.Fe(II)': 'MPE-Fe(II)-structure-mapping-assay',
            'ethylnitrosourea': 'ENU-structure-mapping-assay',
            'lead': 'Lead-structure-mapping-assay',
            'rhodium': 'Rhodium-DNA-structure-mapping-assay',
            'ruthenium': 'Ruthenium-DNA-structure-mapping-assay',
            'terbium': 'Terbium-RNA-structure-mapping-assay',
            'DNAse-1': 'DNASE-1-structure-mapping-assay',
            'RNAse-CL3': 'RNASE-CL3-structure-mapping-assay',
            'nuclease-S1': 'Nuclease-S1-structure-mapping-assay',
            'RNAse-T1': 'RNASE-T1-structure-mapping-assay',
            'RNAse-T2': 'RNASE-T2-structure-mapping-assay',
            'RNAse-U2': 'RNASE-U2-structure-mapping-assay',
            'RNAse-V1': 'RNASE-V1-structure-mapping-assay',
            'OH-radical': 'OH-radical-structure-mapping-assay',
            'inline-Probing': 'Inline-Probing'
        }

ONTOLOGY = Ontology()


ISATAB_INVEST_KEYS = ['ONTOLOGY SOURCE REFERENCE', 'Term Source Name', 'Term Source File', 'Term Source Version', 'Term Source Description', 'INVESTIGATION', 'Investigation Identifier', 'Investigation Title', 'Investigation Description', 'Investigation Submission Date', 'Investigation Public Release Date', 'INVESTIGATION PUBLICATIONS', 'Investigation PubMed ID', 'Investigation Publication DOI', 'Investigation Publication Author list', 'Investigation Publication Title', 'Investigation Publication Status', 'Investigation Publication Status Term Accession Number', 'Investigation Publication Status Term Source REF', 'INVESTIGATION CONTACTS', 'Investigation Person Last Name', 'Investigation Person First Name', 'Investigation Person Mid Initials', 'Investigation Person Email', 'Investigation Person Phone', 'Investigation Person Fax', 'Investigation Person Address', 'Investigation Person Affiliation', 'Investigation Person Roles', 'Investigation Person Roles Term Accession Number', 'Investigation Person Roles Term Source REF', 'STUDY', 'Study Identifier', 'Study Title', 'Study Submission Date', 'Study Public Release Date', 'Study Description', 'Study File Name', 'STUDY DESIGN DESCRIPTORS', 'Study Design Type', 'Study Design Type Term Accession Number', 'Study Design Type Term Source REF', 'STUDY PUBLICATIONS', 'Study PubMed ID', 'Study Publication DOI', 'Study Publication Author list', 'Study Publication Title', 'Study Publication Status', 'Study Publication Status Term Accession Number', 'Study Publication Status Term Source REF', 'STUDY FACTORS', 'Study Factor Name', 'Study Factor Type', 'Study Factor Type Term Accession Number', 'Study Factor Type Term Source REF', 'Study Assay Measurement Type', 'Study Assay Measurement Type Term Accession Number', 'Study Assay Measurement Type Term Source REF', 'Study Assay Technology Type', 'Study Assay Technology Type Term Accession Number', 'Study Assay Technology Type Term Source REF', 'Study Assay Technology Platform', 'Study Assay File Name', 'STUDY PROTOCOLS', 'Study Protocol Name', 'Study Protocol Type', 'Study Protocol Type Term Accession Number', 'Study Protocol Type Term Source REF', 'Study Protocol Description', 'Study Protocol URI', 'Study Protocol Version', 'Study Protocol Parameters Name', 'Study Protocol Parameters Name Term Accession Number', 'Study Protocol Parameters Name Term Source REF', 'Study Protocol Components Name', 'Study Protocol Components Type', 'Study Protocol Components Type Term Accession Number', 'Study Protocol Components Type Term Source REF', 'STUDY CONTACTS', 'Study Person Last Name', 'Study Person First Name', 'Study Person Mid Initials', 'Study Person Email', 'Study Person Phone', 'Study Person Fax', 'Study Person Address', 'Study Person Affiliation', 'Study Person Roles', 'Study Person Roles Term Accession Number', 'Study Person Roles Term Source REF']

ISATAB_ASSAY_OPTIONS = ['Factor Value[chemical]', 'Term Source REF[chemical]', 'Term Accession Number[chemical]', 'Factor Value[chemical concentration]', 'Unit[chemical]']

ISATAB_ASSAY_KEYS = ['Source Name', 'Characteristics[Nucleotide Sequence]', 'Characteristics[Nucleotide Type]', 'Characteristics[RNA Production]', 'Protocol REF', 'Parameter Value[Data Starts at Sequence Position]', 'Term Source REF', 'Term Accession Number', 'Assay Name', 'Raw Data File', 'Performer', 'Date']

ISATAB_INVEST_DICT = {
    'ONTOLOGY SOURCE REFERENCE': [],
    'Term Source Name': ['OBI', 'CHEBI', 'UO'],
    'Term Source File': [],
    'Term Source Version': ['v 1.26', 'v 1.26', 'v 1.26'],
    'Term Source Description': ['Ontology for Biomedical Investigations', 'Chemical Entity of Biological Interest', 'Unit Ontology'],
    'INVESTIGATION': [],
    'Investigation Identifier': [],
    'Investigation Title': [],
    'Investigation Description': [],
    'Investigation Submission Date': [],
    'Investigation Public Release Date': [],
    'INVESTIGATION PUBLICATIONS': [],
    'Investigation PubMed ID': [],
    'Investigation Publication DOI': [],
    'Investigation Publication Author list': [],
    'Investigation Publication Title': [],
    'Investigation Publication Status': [],
    'Investigation Publication Status Term Accession Number': [],
    'Investigation Publication Status Term Source REF': [],
    'INVESTIGATION CONTACTS': [],
    'Investigation Person Last Name': [],
    'Investigation Person First Name': [],
    'Investigation Person Mid Initials': [],
    'Investigation Person Email': [],
    'Investigation Person Phone': [],
    'Investigation Person Fax': [],
    'Investigation Person Address': [],
    'Investigation Person Affiliation': [],
    'Investigation Person Roles': [],
    'Investigation Person Roles Term Accession Number': [],
    'Investigation Person Roles Term Source REF': [],
    'STUDY': [],
    'Study Identifier': [],
    'Study Title': [],
    'Study Submission Date': [],
    'Study Public Release Date': [],
    'Study Description': [],
    'Study File Name': [],
    'STUDY DESIGN DESCRIPTORS': [],
    'Study Design Type': ['SNRNASM'],
    'Study Design Type Term Accession Number': ['OBI'],
    'Study Design Type Term Source REF': [],
    'STUDY PUBLICATIONS': [],
    'Study PubMed ID': [],
    'Study Publication DOI': [],
    'Study Publication Author list': [],
    'Study Publication Title': [],
    'Study Publication Status': [],
    'Study Publication Status Term Accession Number': [],
    'Study Publication Status Term Source REF': [],
    'STUDY FACTORS': [],
    'Study Factor Name': [],
    'Study Factor Type': [],
    'Study Factor Type Term Accession Number': [],
    'Study Factor Type Term Source REF': [],
    'Study Assay Measurement Type': [],
    'Study Assay Measurement Type Term Accession Number': [],
    'Study Assay Measurement Type Term Source REF': [],
    'Study Assay Technology Type': [],
    'Study Assay Technology Type Term Accession Number': [],
    'Study Assay Technology Type Term Source REF': [],
    'Study Assay Technology Platform': [],
    'Study Assay File Name': [],
    'STUDY PROTOCOLS': [],
    'Study Protocol Name': [],
    'Study Protocol Type': [],
    'Study Protocol Type Term Accession Number': [],
    'Study Protocol Type Term Source REF': [],
    'Study Protocol Description': [],
    'Study Protocol URI': [],
    'Study Protocol Version': [],
    'Study Protocol Parameters Name': ['data starts at sequence position'],
    'Study Protocol Parameters Name Term Accession Number': [],
    'Study Protocol Parameters Name Term Source REF': [],
    'Study Protocol Components Name': [],
    'Study Protocol Components Type': [],
    'Study Protocol Components Type Term Accession Number': [],
    'Study Protocol Components Type Term Source REF': [],
    'STUDY CONTACTS': [],
    'Study Person Last Name': [],
    'Study Person First Name': [],
    'Study Person Mid Initials': [],
    'Study Person Email': [],
    'Study Person Phone': [],
    'Study Person Fax': [],
    'Study Person Address': [],
    'Study Person Affiliation': [],
    'Study Person Roles': [],
    'Study Person Roles Term Accession Number': [],
    'Study Person Roles Term Source REF': []
}

ISATAB_ASSAY_DICT = {
    'Source Name': [],
    'Characteristics[Nucleotide Sequence]': [],
    'Characteristics[Nucleotide Type]': [],
    'Characteristics[RNA Production]': [],
    'Protocol REF': [],
    'Parameter Value[Data Starts at Sequence Position]': [],
    'Assay Name': [],
    'Raw Data File': [],
    'Performer': [],
    'Date': [],
    'Term Source REF': [],
    'Term Accession Number': []
}

