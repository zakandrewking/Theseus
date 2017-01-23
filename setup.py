from setuptools import setup, Command, find_packages

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
    packages=find_packages(),
    install_requires=['cobra>=0.5.0', 'cloudpickle>=0.2.2'],
)
