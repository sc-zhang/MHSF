## Introduction

MHSPF(**M**icro**H**omologous **S**equences **P**airs **F**inder) is a simple tool for searching 3-tuples MHS pairs from genome.

## Dependencies

### Software

* Python (>=3.7)

### Python Modules
* pathos

## Installation

```bash
git clone https://github.com/sc-zhang/MHSPF.git
chmod +x mhspf.py
# Optional, add follow line to your .bash_profile or .bashrc
export PATH=/path/to/MHSPF:$PATH
```

## Usage

```bash
usage: mhspf.py [-h] -g GENOME [--min_size MIN_SIZE] [--max_size MAX_SIZE] [-c COUNT] [--min_pair_dist MIN_PAIR_DIST] [--max_pair_dist MAX_PAIR_DIST] [--max_element_dist MAX_ELEMENT_DIST] -o OUTPUT [-t THREADS]

options:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Input genome file
  --min_size MIN_SIZE   Minium size of MHS, default=5
  --max_size MAX_SIZE   Maximum size of MHS, default=7
  -c COUNT, --count COUNT
                        Count of most frequency MHS for each length, default=20
  --min_pair_dist MIN_PAIR_DIST
                        Minimum distance between two MHSs tuples, default=500
  --max_pair_dist MAX_PAIR_DIST
                        Maximum distance between two MHSs tuples, default=3000
  --max_element_dist MAX_ELEMENT_DIST
                        Maximum distance between two MHSs, default=50
  -o OUTPUT, --output OUTPUT
                        Output directory
  -t THREADS, --threads THREADS
                        Threads count, default=10
```

> **Notice:** 
> 1. MHS should contain at least two kind of bases.
> 2. MHSs pairs on genome are like:  
>    MHS1---MHS2---MHS3---------------------MHS1---MHS2---MHS3  
>    (MHSx is a short nucleotide sequence, the threshold of distance between two MHS and between two 3-tuple MHSs can be
>    set by user)
> 3. The 3-tuple MHSs are the top 10% combinations of most frequency MHS
> 4. For each 3-tuple MHSs, we only retain one order of these MHSs randomly.


## Example

```bash
mhspf.py -g genome.fasta -o wrkdir -t 10
```


## Reesults
The output file is a table file seperated by tab, the table file is like below:

| ID | MHS1          | MHS2           | MHS3            |
|----|---------------|----------------|-----------------|
| 1  | AATCC,3923952 | AATGCA,3273232 | TGACCAG,2797233 |

> each record is sequence,frequency
>