![linux](https://github.com/gsiekaniec/ORI/workflows/linux/badge.svg)
![macOS](https://github.com/gsiekaniec/ORI/workflows/macOS/badge.svg)
# <img src="img/ORI.png" alt="ORI" width="3000"/>

ORI (Oxford nanopore Reads Identification) is a software using long nanopore reads to identify bacteria present in a sample at the strain level. 

There are two sub-parts in the ORI program: (1) the creation of the index containing the reference genomes of the interest species and (2) the query of this index with long reads from Nanopore sequencing in order to identify the strain(s). 

The index is based on kmindex [1]



<sub>1. Lemane, Téo, et al. "Indexing and real-time user-friendly queries in terabyte-sized complex genomic datasets with kmindex and ORA." Nature Computational Science 4.2 (2024): 104-109.

----

## Installation
TODO

----

### I) First step: create your own index
```bash
# anticipate ORI needs
python ../ORI.py length -g /projects/holopig/these/mock_pantax/ORI_genomes/ -o lengths.txt

# BF sizes: 
sort -k 2  -n lengths.txt | tail -n 1 
NC_021658.1     14782125        1
# Using https://sebllns.github.io/bf_size_calc/
# bf size is 482443264
# total index size: 9.44 GB


# Create fof
for acc in `ls /projects/holopig/these/mock_pantax/ORI_genomes/*.fna`; do  
	name=`basename $acc | tr '.' '_'`;
	echo "$name: $acc" >> fof.txt;  
done

# build index
kmindex build -i index_km -f fof.txt -d dir_kmindex -r test_pantax -k 15 --hard-min 1 --nb-partitions 32 --bloom-size 482443264
```

### II) Query reads and retreive best set of strains explaining the query resutls


```bash
for reads in `ls /projects/holopig/these/mock_pantax/MADRe_reads/by_species/*.fastq`;
do
	out_ori_name=res_ori_`basename $reads | cut -d '.' -f 1`.txt 
	echo "Creating "$out_ori_name
	head -n 16000 $reads > /tmp/q.fq
	rm -rf output ; time kmindex query -i index_km -q /tmp/q.fq -a -f matrix -r 0.5 -t 40 -a
	time timeout 10m python ../ORI.py identification -m output/test_pantax.tsv -le lengths.txt -c /projects/holopig/rvicedom/env_ori/bin/clingo --output $out_ori_name
done > logs_query_kmindex_ori.txt 2>&1
```
