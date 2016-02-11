from setuptools import setup

setup(name='shand',
      version='0.1',
      description='A pipeline for investigating cospeciation in microbiomes',
      scripts=['scripts/shand'],
      url='http://github.com/ryneches/Shand',
      author='Russell Neches',
      author_email='ryneches@ucdavis.edu',
      license='BSD',
      packages=['shand'],
      install_requires=[
        'pandas',
        'screed',
        'hat_trie',
        'skbio'
      ],      
      zip_safe=False)
