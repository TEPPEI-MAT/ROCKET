# ROCKET
Rational Oligonucleotide design Calculated with Kinetic parameter for Enhanced in vitro Transcription.
ROCKET designed DNA oligonucleotide for in vitro transcription....

## Installation

### pypy3
```python
pypy3 -m pip ROCKET
```
### Conda
```python
conda install -c biocinda ROCKET
```


## Usage

### Python
```python
import ROCKET

### just return Forward and Reverse primer.
ROCKET("GCCGGGGUGGUGUAGCCUGGUUAGCACAGGGGACUGUGGAUCCCCUAGCCCGGGUUCAAAUCCCGGCCCCGGCCCCA")

### if you need to add -1G, CCA end and precursor.
ROCKET("GCCGGGGUGGUGUAGCCUGGUUAGCACAGGGGACUGUGGAUCCCCUAGCCCGGGUUCAAAUCCCGGCCCCGGCCCCA",g=True, c=True, p=True )
```
* g, c and p correspond to addition of -1G, CCA end and precursor respectively.


### CLI
```
usage:usage: ROCKET [-h] [-s --sequence] [-d --dir_path] [-M --maximum_primer] [-m --minimum_primer] [-p --precursor] [-g --g_addition] [-c --cca_addition]

Rational Oligonucleotide design Calculated with Kinetic parameter for Enhanced in vitro Transcription.

optional arguments:
  -h, --help            show this help message and exit
  -s , --sequence       tRNA sequence to transcribe
  -d , --dir_path       Path to fasta file (Output as 'output.txt' in the same directory)
  -M , --maximum_primer  maximum length of primer (default is 100 nt)
  -m , --minimum_primer minimum length of primer (default is 7 nt)
  -p, --precursor       Add precursor to the 5' end
  -g, --g_addition      Add a guanine base to the 5' end
  -c, --cca_addition    Add CCA sequence to the 3' end
```
### Example (When typing a sequence on the CLI.)
```python
python ROCKET -s GCCGGGGUGGUGUAGCCUGGUUAGCACAGGGGACUGUGGAUCCCCUAGCCCGGGUUCAAAUCCCGGCCCCGGCCCCA

#Forward
#GCCTAATACGACTCACTATAGCCGGGGTGGTGTAGCCTGGTTAGCACAGGGGACTGTGGATCCCCTAGCCCGGGTTCAAATCCC
#Reverse
#TGGGGCCGGGGCCGGGATTTGAACCCGGGCTAGGGGATCCAC
#Completed!!
```
### Example (When specifying a fasta file.)
```python
python ROCKET -d C:\Users\username\document\folder\data.fasta
#Completed!!
```
* if you need to add -1G, CCA end and precursor, please use -g, -c and -p respectively.
* if you need sequence longer than 100 oligos and less than 7 oligos, please use -M and -m.
