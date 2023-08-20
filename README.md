# ROCKET
Rational Oligonucleotide design Calculated with Kinetic parameter for Enhanced in vitro Transcription.
ROCKET designed DNA oligonucleotide for in vitro transcription....

## Installation

### pypy3
```python
pypy3 -m pip ROCKET
```



## Usage

### Python
```python
from ROCKET import ROCKET

### just return Forward and Reverse primer. ###
ROCKET("GCCGGGGUGGUGUAGCCUGGUUAGCACAGGGGACUGUGGAUCCCCUAGCCCGGGUUCAAAUCCCGGCCCCGGCCCCA")

### if you need to add -1G, CCA end and precursor.(g, c and p correspond to -1G, CCA end and precursor respectively.) ###
ROCKET("GCCGGGGUGGUGUAGCCUGGUUAGCACAGGGGACUGUGGAUCCCCUAGCCCGGGUUCAAAUCCCGGCCCCGGCCCCA",g=True, c=True, p=True )

```

### CLI
```
usage:usage: ROCKET [-h] [-s --sequence] [-d --dir_path] [-M --maximum_primer] [-m --minimum_primer] [-p --precursor] [-g --g_addition] [-c --cca_addition]

Rational Oligonucleotide design Calculated with Kinetic parameter for Enhanced in vitro Transcription

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
### example
```python
