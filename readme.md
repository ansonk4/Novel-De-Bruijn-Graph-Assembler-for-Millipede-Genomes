# Novel De Bruijn Graph Assembler for Millipede Genomes

This repository is the code for the Novel De Bruijn Graph Assembler for Millipede Genomes project - CUHK Library Data Analytics Practice Opportunity 2022/23

## DBG without clustering
1. Clone this repository

```
git clone https://github.com/ansonk4/DBG.git
```

2. Put the target fasta file to the cloned repository

3. Run dbg.py 

```
cd dbg
python dbg.py
```

4. Input the target fasta file name, for example:

```
Input file name (e.g. file.fasta): T_c_100k.fasta
Input min kmer (e.g: 10): 20
Input max kmer (e.g. 70): 40
```

5.  The resultant genome info and dbg is stored in `dbg_output`

