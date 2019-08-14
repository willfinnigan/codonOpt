from setuptools import setup, find_packages

setup(
    name='codonOpt',
    version='0.23',

    packages=find_packages(),
    package_data={'codonOpt': ['Examples']},

    url='https://github.com/willfinnigan/codonOpt',
    license='MIT',
    author='William Finnigan',
    author_email='wjafinnigan@gmail.com',
    description='A package for codon optimisation',
    keywords = ['DNA', 'Codon Optimisation', 'Gene Design'])
