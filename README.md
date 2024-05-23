# ROCKET
All kinds of RNA molecules can be produced by in vitro transcription using T7 RNA polymerase using DNA templates obtained by solid-phase chemical synthesis, primer extension, PCR, or DNA cloning. The oligonucleotide design, however, is a challenge to nonexperts as this relies on a set of rules that have been established empirically over time. Here, we describe a Python program to facilitate the rational design of oligonucleotides, calculated with kinetic parameters for enhanced in vitro transcription (ROCKET). The Python tool uses thermodynamic parameters, performs folding-energy calculations, and selects oligonucleotides suitable for the polymerase extension reaction. These oligonucleotides improve yields of template DNA. With the oligonucleotides selected by the program, the tRNA transcripts can be prepared by a one-pot reaction of the DNA polymerase extension reaction and the transcription reaction. Also, the ROCKET-selected oligonucleotides provide greater transcription yields than that from oligonucleotides selected by Primerize, a leading software for designing oligonucleotides for in vitro transcription, due to the enhancement of template DNA synthesis. Apart from over 50 tRNA genes tested, an in vitro transcribed self-cleaving ribozyme was found to have catalytic activity. In addition, the program can be applied to the synthesis of mRNA, demonstrating the wide applicability of the ROCKET software.
>(Original article: Matsuda T, Hori H, Yamagami R. Rational design of oligonucleotides for enhanced in vitro transcription of small RNA. RNA. 2024 May 16;30(6):710-727. doi: 10.1261/rna.079923.123. PMID: 38423625.)

__Please cite it if you use !!!__

# Instalation
```python
$ git clone https://github.com/TEPPEI-MAT/ROCKET.git

$ python setup.py develop
```
`Pypy3 environment recomended`
## Usage

### Python
```python
import ROCKET

### just return Forward and Reverse primer.
ROCKET("GCCGGGGUGGUGUAGCCUGGUUAGCACAGGGGACUGUGGAUCCCCUAGCCCGGGUUCAAAUCCCGGCCCCGGCCCCA")

#>Forward
#GCCTAATACGACTCACTATAGCCGGGGTGGTGTAGCCTGGTTAGCACAGGGGACTGTGGATCCCCTAGCCCGGGTTCAAATCCC
#>Reverse
#TGGGGCCGGGGCCGGGATTTGAACCCGGGCTAGGGGATCCAC

### if you need to add -1G, CCA end and precursor.
ROCKET("GCCGGGGUGGUGUAGCCUGGUUAGCACAGGGGACUGUGGAUCCCCUAGCCCGGGUUCAAAUCCCGGCCCCGGCCCCA",g=True, c=True, p=True )

#>Forward
#GCCTAATACGACTCACTATAGGGAGACCACAACGGTTTCCCTCTAGAGGCCGGGGTGGTGTAGCCTGGTTAGCACAGGGGACTGTGGATCCCCTAGC
#>Reverse
#TGGTGGGGCCGGGGCCGGGATTTGAACCCGGGCTAGGGGATCCACAGTCCCCTGTGCTAAC
```
* `g`, `c` and `p` correspond to addition of `-1G`, `CCA end` and `precursor` respectively.
* You can select length of oligonucleotide by using `min_nt` and `max_nt`.
* Also using `min_tm` and `max_tm`, you can decide Tm of oligonucleotides yourself.
 
### CLI
```
usage:usage: ROCKET [-h] [-s --sequence] [-d --dir_path] [-T maximum_tm] [-t minimum_tm] [-M --maximum_primer] [-m --minimum_primer] [-p --precursor] [-g --g_addition] [-c --cca_addition]

Rational Oligonucleotide design Calculated with Kinetic parameter for Enhanced in vitro Transcription.

optional arguments:
  -h, --help            show this help message and exit
  -s, --sequence        tRNA sequence to transcribe
  -d, --dir_path        Path to fasta file (Output as 'output.txt' in the same directory)
  -T, --maximum_tm      maximum Tm (default is 72 degrees C)
  -t, --minimum_tm      minimum Tm (default is 68 degrees C)
  -M, --maximum_primer  maximum length of primer (default is 100 nt)
  -m, --minimum_primer  minimum length of primer (default is 7 nt)
  -p, --precursor       Add precursor to the 5' end
  -g, --g_addition      Add a guanine base to the 5' end
  -c, --cca_addition    Add CCA sequence to the 3' end
``` 
### Example (When typing a sequence on the CLI.)
```python
python ROCKET -s GCCGGGGUGGUGUAGCCUGGUUAGCACAGGGGACUGUGGAUCCCCUAGCCCGGGUUCAAAUCCCGGCCCCGGCCCCA

#>Forward
#GCCTAATACGACTCACTATAGCCGGGGTGGTGTAGCCTGGTTAGCACAGGGGACTGTGGATCCCCTAGCCCGGGTTCAAATCCC
#>Reverse
#TGGGGCCGGGGCCGGGATTTGAACCCGGGCTAGGGGATCCAC
#Completed!!
```
### Example (When specifying a fasta file.)
```python
python ROCKET -d C:\Users\username\document\folder\data.fasta
#Completed!!
```
* If you need to add `-1G`, `CCA end` and `precursor`, please use `-g`, `-c` and `-p` respectively.
* If you need DNA oligonucleotides **longer than 100 oligos** and **less than 7 oligos**, please use `-M` and `-m`.
