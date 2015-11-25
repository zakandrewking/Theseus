try:
    from setuptools import setup, Command
except:
    from distutils.core import setup, Command

setup(
    name='theseus',
    version=0.2,
    author='Zachary King',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    packages=['theseus'],
    install_requires=['cobra>=0.4.0b4'],
)
