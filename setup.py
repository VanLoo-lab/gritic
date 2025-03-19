from setuptools import setup

setup(
    name='gritic',
    version='0.0.5',
    packages=['gritic'],
    author='Toby Baker',
    author_email='tobybaker@mednet.ucla.edu',
    maintainer='Toby Baker',
    maintainer_email='tobybaker@mednet.ucla.edu',
    license='GNU AFFERO GENERAL PUBLIC LICENSE Version 3',
    url='https://github.com/VanLoo-lab/gritic',
    download_url='https://github.com/VanLoo-lab/gritic/archive/refs/tags/v0.0.5.tar.gz',
    description='A module to time complex copy number gains in cancer genomes',
    long_description='''\
A tool for timing complex copy number gains in cancer. Provides gain timing estimates for segments with a total copy number of up to 9. Only copy number segments with 10 or more SNVs will be timed.

Each gain timing is measured in mutation time, a scale that ranges from 0 to 1. A timing of 0 indicates that the gain occured close to conception and 1 that the gain occurred very close to the emergence of the tumour's most recent common ancestor.
''',
    long_description_content_type='text/plain',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: OS Independent',
    ],
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'scipy',
        'scikit-learn',
        'networkx',
        'numba'
    ]
)

