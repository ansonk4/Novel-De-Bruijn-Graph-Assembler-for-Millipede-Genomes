import streamlit as st
import graphviz
from io import StringIO

import dbg 
import dbg_ult

st.title('DBG')

fasta = st.file_uploader('Fasta')


if not fasta:
    st.stop()

minK = st.slider('minK', 5, 50, 10)
maxK = st.slider('maxK', 5, 100, 40)

stringio = StringIO(fasta.getvalue().decode("utf-8"))
reads = dbg_ult.st_read_fasta(stringio)
contig, v, e = dbg_ult.old_testing_kmer_size(minK, maxK, reads)

st.write('Output')
st.write(contig)
st.write('Lenght:')
st.write(len(contig))

    
