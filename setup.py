from setuptools import setup, find_packages

setup(
    name='ROCKET',
    version='1.0',
    packages=find_packages(),
    description="Rational Oligonucleotide design Calculated with Kinetic parameter for Enhanced in vitro Transcription",
    author="Teppei_Matsuda",
    author_email="teppei.matsuda.lab@gmail.com",
    license="MIT",
    url="https://github.com/TEPPEI-MAT/ROCKET",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={"console_scripts":["ROCKET=ROCKET.main:run"]},
    install_requires=['Biopython','argparse']

)
