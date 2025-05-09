conda activate cutadapt

mkdir PLATE_reads/

for i in $(cat primer_combinations.csv | sed "s/,/_/g"); do 
mkdir PLATE_reads/${i};
done

cutadapt -j 1 -e 2 -m 1 -Q 30,30 --action=trim -g file:golay_linker_barcodes.fasta -G file:golay_linker_barcodes.fasta -o PLATE_reads/{name1}_{name2}.fwd.fq.gz -p PLATE_reads/{name1}_{name2}.rev.fq.gz 01_raw_reads/PLATE_R1_001.fastq.gz 01_raw_reads/PLATE_R2_001.fastq.gz

find PLATE_reads/* -type f -size 20c -exec rm {} \;

for i in $(cat primer_combinations.csv); do 
mv PLATE_reads/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,).fwd.fq.gz PLATE_reads/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,)/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,).fwd.fq.gz;
seqkit seq --reverse PLATE_reads/$(echo $i | cut -f 2 -d,)_$(echo $i | cut -f 1 -d,).fwd.fq.gz | gzip >> PLATE_reads/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,)/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,).fwd.fq.gz;
rm PLATE_reads/$(echo $i | cut -f 2 -d,)_$(echo $i | cut -f 1 -d,).fwd.fq.gz;

#mv PLATE_reads/$(echo $i | cut -f 2 -d,)_$(echo $i | cut -f 1 -d,).fwd.fq.gz PLATE_reads/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,)/$(echo $i | cut -f 2 -d,)_$(echo $i | cut -f 1 -d,).fwd.fq.gz;
mv PLATE_reads/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,).rev.fq.gz PLATE_reads/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,)/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,).rev.fq.gz;
seqkit seq --reverse PLATE_reads/$(echo $i | cut -f 2 -d,)_$(echo $i | cut -f 1 -d,).rev.fq.gz | gzip >> PLATE_reads/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,)/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,).rev.fq.gz;
rm PLATE_reads/$(echo $i | cut -f 2 -d,)_$(echo $i | cut -f 1 -d,).rev.fq.gz;

#mv PLATE_reads/$(echo $i | cut -f 2 -d,)_$(echo $i | cut -f 1 -d,).rev.fq.gz PLATE_reads/$(echo $i | cut -f 1 -d,)_$(echo $i | cut -f 2 -d,)/$(echo $i | cut -f 2 -d,)_$(echo $i | cut -f 1 -d,).rev.fq.gz;
done

mkdir PLATE_reads/unknown; mv PLATE_reads/*unknown* PLATE_reads/unknown/
mkdir PLATE_reads/index_hopped; mv PLATE_reads/*.fq.gz PLATE_reads/index_hopped/

