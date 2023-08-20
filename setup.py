from setuptools import setup, find_packages

setup(
    name='ROCKET',
    version='1.0',
    packages=find_packages(),
    description="Rational Oligonucleotide design Calculated with Kinetic parameter for Enhanced in vitro Transcription",
    author="Teppei_Matsuda",
    author_email="",
    license="",
    url="",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=['Biopython','argparse']

)
