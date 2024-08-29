1. Identification of editing sites

# Identify all editing sites
motif_scraper --motif sequence --outputFile cas_species_sites.csv --search_strand=both Species.fasta

# Identify unique editing sites
awk -F "," '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' cas_species_sites.csv > cas_species_sites_chuli.txt
sed -i '1d' cas_species_sites_chuli.txt
awk '{print $6}' cas_species_sites_chuli.txt | sed 's/...$//' > cas_species_sites_sgRNA.txt
paste cas_species_sites_chuli.txt cas_species_sites_sgRNA.txt | sort -k 7 | uniq -f 6 -u > cas_species_sites_uniq.txt

# Identify PAM type
awk '{print $6}' cas9_species_sites_chuli.txt | cut -c 21-23 | sort |uniq -c > cas9_species_sites_kind.txt # Take Cas9 as an example


2. Distribution of editing sites

# Build a bed file of editing sites
awk '{if ($4 == "+") {print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4} else {print $1"\t"$3"\t"$2"\t"$5"\t"$6"\t"$4}}' cas_species_sites_chuli.txt > cas_species_sites.bed
bedtools sort -i cas_species_sites.bed > cas_species_sites.sort.bed

# Build bed files of exon, intron, intergenic and promotorpromoter regions
awk '{if ($3 == "transcript") {print $0}}' Species.gtf > Species.transcript.gtf
sed -i 's/"//g' Species.transcript.gtf
awk 'BEGIN{FS="\t| |;";OFS="\t"}{print $1,$4,$5,$10,$22,$7}' Species.transcript.gtf > Species.transcript.bed
awk '{if ($3 == "exon") {print $0}}' Species.gtf > Species.exon.gtf
sed -i 's/"//g' Species.exon.gtf
awk 'BEGIN{FS="\t| |;";OFS="\t"}{print $1,$4,$5,$10,$22,$7}' Species.exon.gtf > Species.exon.bed
bedtools sort -i Species.transcript.bed > Species.transcript.sorted.bed
bedtools sort -i Species.exon.bed > Species.exon.sorted.bed
bedtools subtract -s -a Species.transcript.sorted.merged.bed -b Species.exon.sorted.merged.bed > Species.intron.bed
bedtools subtract -a Species.genome.bed -b Species.transcript.sorted.merged.
bed > Species.intergenic.bed
awk '{if ($6 == "+") {print $1"\t"$2-1000"\t"$2-1"\t"$4"\t"$5"\t"$6} else {print $1"\t"$3+1"\t"$3+1000"\t"$4"\t"$5"\t"$6}}' Species.transcript.sorted.bed > Species.promotor.bed
awk '{if ($2<1) {print $1"\t"1"\t"$3"\t"$4"\t"$5"\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' Species.promotor.bed > Species.promotor.chuli.bed
bedtools sort -i Species.promotor.chuli.bed > Species.promotor.sorted.bed
bedtools merge -i Species.promotor.sorted.bed -c 4,5,6 -o distinct > Species.promotor.sorted.merged.bed

# Identify the distribution of editing sites
intersectBed -a cas_species_sites.sort.bed -b Species.exon.sorted.merged.bed -wa > cas_species_sites_exon.bed
intersectBed -a cas_species_sites.sort.bed -b Species.intron.bed -wa > cas_species_sites.intron.bed
intersectBed -a cas_species_sites.sort.bed -b Species.intergenic.bed -wa > cas_species_sites.intergenic.bed
intersectBed -a cas_species_sites.sort.bed -b Species.promotor.sorted.merged.bed -wa > cas_species_sites_promotor.bed

# Judge if any editing sites were in exon or promoter
python judge_gene_editing_site.py Species.transcript.bed cas_species_sites_exon_gene.bed cas_species_sites_exon_gene.txt
python judge_gene_editing_site.py Species.transcript.bed cas_species_sites_promotor_gene.bed cas_species_sites_promotor_gene.txt

3. Variation sites analysis

# Build bed files of variations (take 1000 Genomes Project [human] as an example)
zcat ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz | awk '{if (!/#/) {print $1"\t"$2"\t"$2}}' > human_variation_sites.bed
sort -u human_variation_sites.bed > human_variation_sites_chuli.bed
bedtools sort -i human_variation_sites_chuli.bed > human_variation_sites.sorted.bed

# Identify the variations in editing sites
intersectBed -a human_variation_sites.sorted.bed -b cas_human_sites.bed -wa > cas_human_variation_sites_in_editing_sites.bed
intersectBed -a human_variation_sites.sorted.bed -b cas_human_sites.bed -wb > cas_human_editing_sites_with_variation_sites.bed
intersectBed -a human_variation_sites.sorted.bed -b cas_human_sites_exon.bed -wa > cas_human_variation_sites_in_editing_sites_exon.bed
intersectBed -a human_variation_sites.sorted.bed -b cas_human_sites_exon.bed -wb > cas_human_editing_sites_with_variation_sites_exon.bed
intersectBed -a human_variation_sites.sorted.bed -b cas_human_sites_promotor.bed -wa > cas_human_variation_sites_in_editing_sites_promotor.bed
intersectBed -a human_variation_sites.sorted.bed -b cas_human_sites_promotor.bed -wb > cas_human_editing_sites_with_variation_sites_promotor.bed

# Calculate the number of variation sites in editing sites
python get_variation_site_in_editing_sites.py cas_human_editing_sites_with_variation_sites.bed cas_human_editing_sites_with_variation_sites_site.txt
python get_variation_site_in_editing_sites.py cas_human_editing_sites_with_variation_sites_exon.bed cas_human_editing_sites_with_variation_sites_exon_site.txt
python get_variation_site_in_editing_sites.py cas_human_editing_sites_with_variation_sites_promotor.bed cas_human_editing_sites_with_variation_sites_promotor_site.txt















