from de_bruijn import *
import time

def read_fasta(file_name: str):
    reads = []
    with open(file_name, encoding='utf-8') as file:
        i = 0
        seq = ""
        for line in file.readlines():
            if not line.startswith('+') and not line.startswith('@') and not line.startswith('#') and not line.startswith('>'):
                seq += line[:-1]
            else:
                reads.append(seq)
                i += 1
                seq = ""
    
    reads.pop(0)

    return reads

def assembly(reads, k):
    """this function produces a list of assembly results 
        from a set of given sequencing reads

    Args:
        reads (List(String)): a list of strings of sequencing reads
        k (int): the string size for the de bruijn graph nodes

    Returns:
        final_contigs: a dictionary mapping the actual contigs 
        and their corresponding reads
    """
    start_time = time.time()
    vertices, edges = construct_graph(reads, k=k)
    contigs = output_contigs(vertices, edges)
    final_contigs = find_contigs(contigs, reads)
    
    print(f'Time Take: {time.time() - start_time}')
    return final_contigs, vertices, edges

reads = read_fasta("T_c_100k.fasta")
contigs, vertices, edges = assembly(reads, k=10)
print(list(contigs.keys())[0])
print(len(contigs))