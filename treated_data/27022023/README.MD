# Treatment of samples and metadata

Some database samples downloaded from the MG-RAST contains incorrect metadata information. 
To address this issue, we have manually reevaluated the sample metadata for accuracy. 
Mistakes include contradictions between sequencing type and sample name, 
incorrect coordinates, and discrepancies between sample categories metadata and 
project descriptions. Our goal is to provide more accurate data information.

---

## Procedures for metadata reavaliation

1. We used a python script to download all samples metadata from the MG-RAST api.
Then we have selected the important metadata columns to reavaliate. 
The considereted columns were the sample name,feature, biome, material, "env_package", 
"metagenome_taxonomy", coordinates and country. 

2. For the manual repair of the contradictions we have analyzed first the individual
column of the sample name to remove samples with contigs, scaffolds and assemblies.
Futhermore we have filtered samples with the name "RNA" and "cDNA".

3. MUST WRITE THE METHODOLGY OF THE COMPARISON BETWEEN THE SAMPLES CATEGORIES METADATA 
(FEATURE, BIOME, MATERIAL, ENV_PACKAGE, SAMPLE_NAME) FOR ANALYZING MISTAKES IN THE SAMPLE CLASSIFICATION.
(EDU)

4. To detect incorrect coordinates we reviewed the project description in combination 
to metadata categories such as feature, biome, material, metagenome taxonomy, 
"env_package" and if necessary the country.


---

## Repository organization

1. The base_table.csv file is a table of the current state of our metagenomic 
metadata. This table will be should be divided into biomes for a manual repair
of incongruences.

2.  life-style, habitat, biosphere and ecosystem are thee problematic 
categorizations wwwe should review. In addition, wwe should evaluate all details
such as accuracy of geographic coordinates

3. the folder biome_tables divide the data by biomes got from MG-RAST original
metadata.
