# ROCKET
In vitro transcription from a template DNA by T7 RNA polymerase can synthesize any class of RNA molecule such as mRNA, rRNA, tRNA, small RNA, or long non-coding RNA. The technique is broadly used in RNA-related research including RNA therapeutics, biochemistry, structural biology, and molecular biology. The template DNA for the transcription reaction can be synthesized from DNA oligonucleotides by a polymerase extension reaction. This is particularly useful for preparing templates encoding small RNA genes. The oligonucleotide design, however, is a challenge to non-experts as this relies on a set of rules that have been established empirically over time. Here, we describe a Python program to facilitate the Rational design of Oligonucleotides, Calculated with Kinetic parameters for Enhanced in vitro Transcription (ROCKET). The Python tool employs thermodynamic parameters, performs folding-energy calculations, and selects oligonucleotides suitable for in vitro transcription.  These oligonucleotides improve yields of template DNA. With optimized reaction conditions and oligonucleotides singled out by the tool, the yield of tRNA transcription increased up to ~5-fold. Apart from over 50 tRNA genes tested, an in vitro transcribed self-cleaving ribozyme was found to have catalytic activity, demonstrating the wide applicability of the ROCKET software, which has been made available at https://github.com/TEPPEI-MAT/ROCKET.
(Teppei Matsuda, Hiroyuki Hori*, and Ryota Yamagami*:under submission)

## Installation

### Conda
```python
```

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
* if you need to add `-1G`, `CCA end` and `precursor`, please use `-g`, `-c` and `-p` respectively.
* if you need DNA oligonucleotides **longer than 100 oligos** and **less than 7 oligos**, please use `-M` and `-m`.
