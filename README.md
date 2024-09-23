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
cat samples.txt| xargs -I {} -P 10 sh -c "fastq2vcf.py all -1 {}_2.fastq.gz -2 {}_2.fastq.gz --ref /mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta --p {} --threads 10 --bqsr-vcf	/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/3d7_hb3.combined.final.vcf.gz,/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/7g8_gb4.combined.final.vcf.gz,/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/hb3_dd2.combined.final.vcf.gz	--cram"
conda deactivate
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
Now it is time to assess coverage...
```
# Run Mosdepth, example in Zambia directory
cd ~/Projects/July_24_Pf/Zambia
cat samples.txt| xargs -I {} -P 10 sh -c "mosdepth -n --fast-mode --by 500 {} {}.cram"

```

