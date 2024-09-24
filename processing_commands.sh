#Quality control and assembly
java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

megahit -1 output_forward_paired.fq -2 output_reverse_paired.fq -o megait_out

# viral contigs detection
genomad end-to-end contigs.fna genomad_output_dir genomad_db

# Clustring viral contigs into vOTUs
mmseqs easy-cluster viral_contigs.fna vOTU tmp --min-seq-id 0.9 -c 0.9 --cov-mode 1

# Call genes
prodigal -p meta -a vOTU_prot.faa -m -d vOTU_gene.fna -o genes.gff -f gff -s poteintial.stat -i vOTU.fna

# Relative abundance estimation
## mapping with bowtie2
bowtie2 -build vOTU.fna vOTU_db 
bowtie2 -x vOTU_db -S abun.sam -1 forward.fq -2 reverse.fq --no-unal --very-sensitive
## abundance estimation using checkM
samtools view -bS abun.sam > abun.bam
samtools sort abun.bam -o abun.sorted.bam
samtools index abun.sorted.bam
checkm coverage -x fna vOTU_dir/ vOTU_abun.out abun.sorted.bam

# RdRP extraction
psiblast -query contigs_6frame_translation.faa -db rdrp_ref_v3 -o match -num_iterations 3
hmmsearch -o mot.XX.match mot.XX.afa.hmm rdrp_draft.faa
# Libraries of motif A\B\C\D were searched separately and RdRP drafts cover at least 3 of them were considered as complete core RdRP.

# Tree buiding
muscle -align core_rdrp.faa -output core_rdrp.afa
hhconsensus -M 99 -i core_rdrp.afa -oa3m core_rdrp_with_consensus.a3m
fasttree core_rdrp_merged_by_consensus.afa > tree.nwk
