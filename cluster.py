from Bio import SeqIO
from itertools import groupby
from dbg import test_kmer, visualize_debruijn
import time
'''This code:
    1. read and extract the id from each cluster
    2. each cluster do DBG,
    3. output fasta for all contigs'''

def read_clstr():
    with open("T_c_1000_filter_CDH_80%.clstr", "r") as fcd:
        lines = fcd.readlines()

    clusters = {}
    read_list = []
    clstr_id = None
    for line in lines:
        if line[0] == ">":
            if len(read_list) > 0:
                clusters[clstr_id] = read_list
            read_list = []
            clstr_id = line.strip()[1:] 
        else: 
            read_id = line.strip().split(",")[1]
            read_id = read_id.split(" ")[1].replace(">", "").replace("...", "") 
            read_list.append(read_id)

    return clusters 
        

def run_dbg_for_cluster(all_dict: dict):
    fasta_file = "T_c_100k.fasta"
    index = SeqIO.index(fasta_file, "fasta")

    clstr_dict = {}
    for i in all_dict:
        id_list = all_dict[i] 
        seq_list = []
        for id in id_list:
            if id in index:
                sequences = index[id].seq
                seq_list.append(str(sequences))            
        clstr_dict[i] = seq_list

    no_clusters = len(clstr_dict)
    with open(f"sequences_all_contigs_.fasta", "w") as file: 
        for i, c in enumerate(clstr_dict):
            reads = clstr_dict[c]
            print(f'Cluster {i+1} / {no_clusters}')
            if (len(reads) == 1):
                contigs = reads[0]
            else:
                contigs, v, e, k  = test_kmer(10, 40, reads)
                # visualize_debruijn(v, e, False, f'cluster{i}', 'cluster_dbg')
            file.write(f">Final_contigs_cluster_{i} Lenght={len(contigs)}\n") 
            file.write(contigs + "\n")
            print()

if __name__ == '__main__':
    start = time.time()
    clusters = read_clstr()
    run_dbg_for_cluster(clusters)
    print(f'Time Taken: {time.time() - start}s')