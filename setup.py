__author__ = 'dtgillis'

from distutils.core import setup

setup(name='cnvpy',
      version='1.0-dev',
      description='Copy Number Variant Python',
      author='Daniel Gillis',
      author_email='dtgillis@renci.org',
      url='www.somewhere.com',
      py_modules=['cnvpy.cnv_caller'],
      packages=['cnvpy', 'cnvpy.depth_coverage',
                'cnvpy.samtools_utils', 'cnvpy.util_objects'],
      requires=['numpy', 'scipy', 'pysam'])


