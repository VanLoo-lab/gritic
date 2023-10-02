from setuptools import setup, find_packages

setup(
    name='gritic',
    version='0.0.1',
    packages=['gritic'],
    author='Toby Baker',
    author_email='toby.baker@crick.ac.uk',
    maintainer='Toby Baker',
    maintainer_email='toby.baker@crick.ac.uk',
    license='GNU AFFERO GENERAL PUBLIC LICENSE Version 3',
    url='https://github.com/VanLoo-lab/gritic',
    download_url='hh
    description='A module to time complex copy number gains in cancer genomes',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU AFFERO GENERAL PUBLIC LICENSE Version 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'scipy',
        'sklearn',
        'networkx',
        'numba'
    ]
)

