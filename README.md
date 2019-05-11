# BA-Thesis
My BA Thesis concerning bioinformatical analysis of IDP (intrinsically disordered proteins)

# Programms:
## divider.py
This script returns n-meres sequences of proteins which were provided. 
It takes as an input .fasta file with many protein sequences and n which is the length of meres.
Additionaly this program counts the number of occurances of all meres.

As an example:
```bash
./divider.py -i text.fasta -n 5 
```
returns a file 
```bash
cat test_5-meres
AAAAA	5490
LLLLL	3391
HTGEK	3059
TGEKP	2937
GEKPY	2593
SSSSS	2589
QQQQQ	2407
CGKAF	2075
IHTGE	1731
GGGGG	1363
ECGKA	1358
...
```


## someStats.py
This program is for anylysis of short n-meres-peptides in proteins. It returns z-scores of n-meres-peptides based on the size of its permutation group. 
Input file should contain nmeres with its counts in dataset.

```bash
./someStats.py -i test_5-meres -n 5
```
Output file is all possible n-mers with its z-score.
```bash
cat test_5-meres_ZScores
HTGEK	496.95
TGEKP	415.91
GEKPY	407.37
CGKAF	395.1
IHTGE	312.81
ECGKA	269.92
RIHTG	252.2
KPYEC	223.42
...
```

