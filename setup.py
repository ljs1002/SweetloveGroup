
from distutils.core import setup
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = dir_path+"/SweetloveGroup/"

setup(name='SweetloveGroup',
      version='1.0',
      description='',
      author='Sanu Shameer',
      author_email='sanushameer@gmail.com',
      url='https://github.com/sshameer/SweetloveGroup/',
      py_modules=[dir_path+'/FVA',dir_path+'/FBA'],
      license='GNU v3',
      packages=['SweetloveGroup'],
      install_requires=[
	      'cobra',
	      ]
     )
