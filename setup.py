
from distutils.core import setup

setup(name='RDATkit',
        version='0.5',
        description='RNA dataset toolkit',
        author='Pablo Cordero',
        author_email='tsuname@stanford.edu',
        url='http://rmdb.stanford.edu/rdatkit',
        packages=['rdatkit', 'rdatkit.likelihood', 'rdatkit.mutate_and_map'],
        py_modules=['secondary_structure', 'mapping', 'datahandlers', 'rna', 'view', 'ontology'],
        requires=['xlwt', 'xlrd', 'biopython']
        )
