
from setuptools import setup, find_packages
setup(
    name="Seqlib",
    version="0.1",
    packages=find_packages(),
    entry_points={
        'console_scripts': ['Seqlib = Seqlib.__main__:main']
        }
)