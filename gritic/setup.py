from setuptools import setup, find_packages

setup(
    name='gritic',
    version='0.1',
    packages=find_packages(),
    author='Toby Baker',
    author_email='toby.baker@crick.ac.uk',
    maintainer='Toby Baker',
    maintainer_email='toby.baker@crick.ac.uk',
    description='A module for GRITIC analysis',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
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

