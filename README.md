# _P. falciparum_ Dataset: Curation and Analysis
All code used in this repository is dedicated to the curation and population genetics analysis of a global _P.falciparum_ dataset. This can help us with the creation and upkeep of big data in the future! If you have any reccomendations, please feel free to raise a pull request or drop me an email. :smiley_cat:

## Downloading the samples from ENA
The first step is the addition and curation of new samples from the dataset and automatically running the fastq2matrix pipeline developed by **Jody Phelan**. This code was developed by **Joseph Thorpe** and is useful for speeding up the curation process. Using the Advanced Search of the ENA browser (Search -> Advanced Search), select Data type: raw reads and either search for a project accession or another 'wgs' query with taxa id. At the bottom, you should see a 'Copy Curl Request'. Copy this once you are happy with the search terms. Note: I would suggest reading through studies associated with ENA data to ensure this is data you would like to include into the dataset. Such studies may also include useful metadata that you can download and manually curate. Once you have copied the curl request replace the existing curl sequest in the f2m/run_get_new_samples.txt. You will also want to change any paths to where you want to store the files (Guyana directory in the example provided) and paths to your 'f2m' folder containing Joe's super handy code! 

```
# Some installation (fastq2matrix)

```
