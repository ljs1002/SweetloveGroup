
from distutils.core import setup
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

setup(name='SweetloveGroup',
      version='1.0',
      description='',
      author='Sanu Shameer',
      author_email='sanushameer@gmail.com',
      url='https://github.com/sshameer/SweetloveGroup/',
      packages=[dir_path+'../SweetloveGroup'],
     )
