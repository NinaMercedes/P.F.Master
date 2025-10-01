# _P. falciparum_ Dataset: Curation and Analysis :dna:
All code used in this repository is dedicated to the curation and population genetics analysis of a global _P.falciparum_ dataset. This can help us with the creation and upkeep of big data in the future! If you have any reccomendations, please feel free to raise a pull request or drop me an email or Teams message. :smiley_cat:

## Downloading the samples from ENA :arrow_down:
The first step is the addition and curation of new samples from the dataset and automatically running the fastq2matrix pipeline developed by **Jody Phelan**. This code was developed by **Joseph Thorpe** and is useful for speeding up the curation process. Using the Advanced Search of the ENA browser (Search -> Advanced Search), select Data type: raw reads and either search for a project accession or another 'wgs' query with taxa id. At the bottom, you should see a 'Copy Curl Request'. Copy this once you are happy with the search terms. Note: I would suggest reading through studies associated with ENA data to ensure this is data you would like to include into the dataset. Such studies may also include useful metadata that you can download and manually curate. Once you have copied the curl request replace the existing curl sequest in the f2m/run_get_new_samples.txt. You will also want to change any paths to where you want to store the files (Guyana directory in the example provided) and paths to your 'f2m' folder containing Joe's super handy code! 

```
# Some installation (fastq2matrix)
cd ~
conda create -n fastq2matrix
conda install python=3.7 bwa samtools bcftools parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree plink
git clone https://github.com/LSHTMPathogenSeqLab/fastq2matrix.git
cd fastq2matrix
python setup.py install
cd ~
conda deactivate
```

Now to download the samples and generate some cram and g.vcf.gz files... These will be aligned to out Pf3D7 reference sequence (using bwa) and variants are called using GATK software. Additional bits and bobs such as trimming (trimmomatic) and qc are performed in this pipeline. We use bqsr, using *P. falciparum* crosses to optimise the variant calls. Any changes you wish to make can be done by editing the run_fastq2matix.sh file. I have also provided an example for in-house sequences (uploaded to ENA once published). :heavy_exclamation_mark: For these **in-house new sample sequences**, please put them in a new directory in new_samples_09_24_v2 and fill out the excel spreadsheet found in the directory. 
```
# Download and process raw read data
conda activate fastq2matrix
cd Pf_09_24_v2/f2m
bash run_get_new_samples.sh
# You can also run_new_sample.sh with a single accession to make this work for one individual sample.
# For in-house files, the following code can be run in a dedicated directory (recommend using xargs for parallelisation)
cat samples.txt| xargs -I {} -P 10 sh -c "fastq2vcf.py all -1 {}_1.fastq.gz -2 {}_2.fastq.gz --ref /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta --p {} --threads 10 --bqsr-vcf
/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/3d7_hb3.combined.final.vcf.gz,/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/7g8_gb4.combined.final.vcf.gz,/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/hb3_dd2.combined.final.vcf.gz	--cram --bam-qc"
# Note the combination of --cram and --bam-qc may not work, could run:
cat samples.txt| xargs -I {} -P 10 sh -c "fastq2vcf.py all -1 {}_1.fastq.gz -2 {}_2.fastq.gz --ref /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta --p {} --threads 10 --bqsr-vcf
/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/3d7_hb3.combined.final.vcf.gz,/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/7g8_gb4.combined.final.vcf.gz,/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/hb3_dd2.combined.final.vcf.gz	--bam-qc" #or just with --cram for larger datasets

# If you are interested in a specific accession list you can also run thr following code, remember again to change the file paths
cat samples.txt| xargs -I {} -P 10 sh -c "bash run_new_sample.sh {}"
conda deactivate
```

## Get Multi-sample BAM from Pf7k and Indonesia 📁
Some samples have multiple runs that we have merged. To double check compare results to those mergeed in the Pf7k project. Also downloaded bams from Indonesia because no fastq are available. Ran fastq2 matrix on these with --bam option.
```
cd ~/ Pf_09_24_v2
mkdir multi_sample
cat multi_samples.txt| xargs -I {} -P 10 sh -c "wget https://pf7_release.cog.sanger.ac.uk/bam/{}.bam"
cat Indonesia.txt| xargs -I {} -P 10 sh -c "wget https://pf7_release.cog.sanger.ac.uk/bam/{}.bam"
#Index
cat multi_samples.txt| xargs -I {} -P 10 sh -c "wget https://pf7_release.cog.sanger.ac.uk/bam/{}.bam.bai"
cat Indonesia.txt| xargs -I {} -P 10 sh -c "wget https://pf7_release.cog.sanger.ac.uk/bam/{}.bam.bai"
#fastq2matrix
cat multi_samples.txt| xargs -I {} -P 10 sh -c "fastq2vcf.py gatk --bam {}.bam --ref /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta --p {} --threads 10"
```

## Coverage :five:
Once we have the g.vcf.gz and cram files, we want to assess the genome coverage. For the dataset we will aim for 60% of genome coverage greater than 5. We will assess this using mosdepth and the following code. First install mosdepth- I reccomend installing in its own conda just to separate out the steps. Code based off previous analysis by **Emilia Manko**. 
```
# Some Installation
conda create -n mosdepth
conda activate mosdepth
conda install bioconda::mosdepth
conda install conda-forge::r-base
# Open R and install data.table, dplyr, ggplot2, readr
```
Now it is time to assess coverage using 500bp windows... You can use 1000 windows instead, you will just need to alter '500' in the parse results script to reflect this change. Instead of a number you can also supply a specific bed/ tabulated file for coverage estimates over a specific region. Also you can change the coverage threshold/ percentage. 
```
# Run Mosdepth, example in Zambia directory
cd ~/Projects/July_24_Pf/Zambia
cat samples.txt| xargs -I {} -P 10 sh -c "mosdepth -n --fast-mode --by 500 {} {}.bqsr.cram --fasta /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta"
cat samples.txt| xargs -I {} -P 10 sh -c "mosdepth -n --fast-mode --by "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Core_genome_Pf3D7_v3_ext.bed" {} {}.bqsr.cram --fasta /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta"
# Change sample file paths etc and run for a few samples, with chr graphs
Rscript code/coverage/parse_mosdepth.R
# or the following for batch/ large-scale analysis
Rscript code/coverage/parse_mosdepth_coverage.R
# Once we have the samples we need, move them over to the big cram and vcf directories!
conda deactivate
```


## Merging VCFs :heavy_plus_sign:
We can perform joint variant calling when merging our vcfs to optimise our variant calls. Code below- the first part is using our dummy dataset to check everything works. We don't want to start merging only to realise we don't have all the correct packages installed: check you have GATK, samtools, bedtools, bcftools, plink and R. Make sure to run in screen and linux encoded txt file!
```
conda activate fastq2matrix
# to test using smaller dataset in 'dummy' directory
cd ~Pf_09_24/dummy
"/mnt/storage13/nbillows/fastq2matrix/scripts/merge_vcfs.py" all --sample-file wgs_id.txt --prefix dummy_data --ref /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta --vcf-dir "/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/vcfs/" --threads 10
conda deactivate

```

## Filtering multi-sample VCF
Once we have successfully merged the VCF (which can take an age for large datasets), we can now begin the filtering process. The parameters/ options used here are standard for our *P. falciparum* pipeline, but feel free to change as you wish. 
```
# Again using the dummy file as an example
cd ~Pf_09_24/dummy
conda activate fastq2matrix
filter_merged_vcf.py --merged-file "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.genotyped.vcf.gz" --prefix dummy_data.2024_10_29.filt --ref "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta" --bqsr-vcf	/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/3d7_hb3.combined.final.vcf.gz,/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/7g8_gb4.combined.final.vcf.gz,/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/hb3_dd2.combined.final.vcf.gz	 --include-region /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Core_genome_Pf3D7_v3_ext.bed --vqslod 0 --missing-sample-cutoff 0.4 --cutoff-mix-GT 0.8 --gff-file /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.modified.new.gff3 --threads 8
```


```
python "/mnt/storage13/nbillows/Pf_09_24/dummy/DR_haplotypes.py" --vcf "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.GT.miss0.4.vqslod.filt.snps.vcf.gz" --prefix dummey_test --reference "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta" --gff "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.modified.new.gff3" --DR_bed "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/DR.tab"
```
