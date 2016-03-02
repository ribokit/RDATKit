from setuptools import setup, find_packages

from rdatkit.__init__ import __version__

setup(
    name='rdatkit',
    description='RNA Dataset Toolkit',
    keywords='RNA Data RDAT Kit',
    version=__version__,

    author='Siqi Tian, Pablo Cordero, Rhiju Das',
    author_email='rhiju@stanford.edu',

    url='https://github.com/hitrace/rdatkit/',
    license='https://rmdb.stanford.edu/rdatkit',

    packages=find_packages(),
    install_requires=open('requirements.txt', 'r').readlines(),
    classifiers=(
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7'
    )
)
