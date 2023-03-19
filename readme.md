# Novel De Bruijn Graph Assembler for Millipede Genomes

This repository is developed for [Novel De Bruijn Graph Assembler for Millipede Genomes - CUHK Library Data Analytics Practice Opportunity 2022/23](https://dsprojects.lib.cuhk.edu.hk/en/projects/dbg-genome/home/?edit&language=en)

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

4. Input the info, for example:

```
Input file name (e.g. T_c_100k.fasta): T_c_100k.fasta
Input min kmer (e.g: 10): 10
Input max kmer (e.g. 40): 40
Input the number of genome to read (e.g. 1000): 1000
```

The above input means the program performs genome assembly on the first 1000 genome of TC 100k.fasta, and tests the kmer size from 10 to 40.

5.  The resultant genome and dbg is stored in `dbg_output/`

## DBG with clustering
1. Follow the instruction here to obtain the `*.clstr` file

2. Clone this repository

```
git clone https://github.com/ansonk4/DBG.git
```

3. Run cluster.py to generate `clustred_fasta.fasta`

```
cd dbg
python cluster.py
```

4. Input the info, for example:

```
Input the clstr file (e.g. T_c_1000_filter_CDH_80%.clstr): T_c_1000_filter_CDH_80%.clstr
Input the fasta file (e.g. T_c_100k.fasta): T_c_100k.fasta  
```

5. Use `clustred_fasta.fasta` as the input to run dbg.py

```
python cluster.py
```

6. Input the info, for example:

```
Input file name (e.g. T_c_100k.fasta): clustred_fasta.fasta
Input min kmer (e.g: 10): 10
Input max kmer (e.g. 40): 40
Input the number of genome to read (e.g. 1000): 1000
```

7. The resultant genome and dbg is stored in `dbg_output/`

## Test case
`T_c_100k.fasta` and `T_c_1000_filter_CDH_80%.clstr` are provided for testing purposes. They are obtained from `Trigoniulus_corallinus_genomic.fastq`. For more info, please go to [Novel De Bruijn Graph Assembler for Millipede Genomes - CUHK Library Data Analytics Practice Opportunity 2022/23](https://dsprojects.lib.cuhk.edu.hk/en/projects/dbg-genome/home/?edit&language=en)