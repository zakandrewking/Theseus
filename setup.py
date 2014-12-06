try:
    from setuptools import setup, Command
except:
    from distutils.core import setup, Command

setup(name='Theseus',
      version=0.1,
      author='Zachary King',
      packages=['theseus'],
      package_data={'escher': ['data/models', 'data/model_pickles']})
