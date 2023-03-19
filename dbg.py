import math # math.inf
import time
import graphviz 
from multiprocessing import Pool #starmap
from itertools import repeat
from copy import deepcopy 

def build_k_mer(reads: list[str], k: int) -> tuple[list, list]:
    k_mers, starts = [], []
    for read in reads:
        k_mers += [read[i:i+k] for i in range(len(read)-k+1)]
        starts.append(read[0:k-1])

    return k_mers, starts

def debruijnize(kmers: list) -> tuple[dict, dict]:
    nodes, edges = dict(), dict() 
    for kmer in kmers:
        kmer1 = kmer[:-1]
        kmer2 = kmer[1:]
        if kmer1 not in nodes:
            nodes[kmer1] = kmer1
        if kmer2 not in nodes:
            nodes[kmer2] = kmer2

        if kmer1 not in edges:
            edges[kmer1] = []
        edges[kmer1].append(kmer2)

    return nodes, edges

def construct_eulerian_path(edges: dict[str: list[str]], starts: list[str]) -> str:
    best_contig = ""
    best_length = 0

    for start in starts:
        E = deepcopy(edges)
        contig = start
        current = start

        while current in E:
            next = E[current][0]
            contig += next[-1]
            del E[current][0]
            if not E[current]:
                del E[current]
            current = next

        if len(contig) > best_length:
            best_contig = contig
            best_length = len(contig)

    return best_contig 

def assembly(reads: list, k: int) -> tuple[str, dict, dict, int]:
    k_mers, starts = build_k_mer(reads,k)
    nodes, edges = debruijnize(k_mers)
    trail = construct_eulerian_path(edges, starts)

    return trail, nodes, edges, k 


def single_kmer_testing(minK: int, maxK: int, step: int, reads:list):
    test_range = list(range(minK, maxK, step))
    with Pool() as pool:
        results = pool.starmap(assembly, zip(repeat(reads), test_range))

    max_length = -math.inf
    best_id = -1
    for i, result in enumerate(results):
        print(f"K: {result[3]} length: {len(result[0])}")
        if len(result[0]) >= max_length:
            best_id = i
            max_length = len(result[0])

    return results[best_id] 

def test_kmer(minK: int, maxK: int, reads:list):
    """testing kmer size in range [minK, minK]""" 
    start_time = time.time()
    
    print("First test")
    contigs, vertices, edges, bestK = single_kmer_testing(minK, maxK+1, 5, reads)

    print("Second test")
    minK = max(minK, bestK-5)
    maxK = min(maxK, bestK+5)
    contigs, vertices, edges, bestK = single_kmer_testing(minK, maxK, 1, reads)

    print(f"Best Kmer size: {bestK}")
    print(f'Time Taken: {time.time() - start_time}s')

    return contigs, vertices, edges, bestK

def visualize_debruijn(Veritces, Edges, view: bool, file_name = 'DBG', dir = 'dbg_output'):
    dot = graphviz.Digraph(file_name, comment='DBG', engine='sfdp', node_attr={'color': 'lightblue2', 'style': 'filled', 'width':'0.5', 'shape':'point'})
    dot.attr(overlap='false', beautify='true', bgcolor='black')
    dot.edge_attr.update(arrowsize='1.5', color='white')
    for vertice in Veritces:
        dot.node(vertice, '')
    for vertice, neighbours in Edges.items():
        for neighbour in neighbours:
            dot.edge(vertice, neighbour)

    dot.render(directory=dir, view=view)  
    return dot

def read_fasta(file_name: str, num_seq: int):
    print(f'Reading {file_name}')
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
                if i > num_seq:
                    break
    
    reads.pop(0)
    print('Finish reading fasta file')
    return reads

if __name__ == '__main__':
    file = input("Input file name (e.g. file.fasta): ")
    kmin = int(input("Input min kmer (e.g: 10): "))
    kmax = int(input("Input max kmer (e.g. 70): "))
    reads = read_fasta(file, 100)
    contigs, vertices, edges, bestK = test_kmer(kmin, kmax, reads)
    visualize_debruijn(vertices, edges, True)
    with open('dbg_output/result.txt', "w+") as file:
        file.write(f'Best k_mer size: {bestK}\n')
        file.write(f'Length of the resultance sequence: {len(contigs)}s\n')
        file.write(f'Sequence: {contigs}')

