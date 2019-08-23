##  Proteomics AMR detection ## 

### The project contains scripts for the manuscript on AMR inference from proteomics.<br />
1. AMR_00_main.r is the main script that sources scripts for different components in the manuscript.<br />
2. AMR_01_RGI.r summarizes results for isolates from RGI (Resistance Gene Identifier).<br />
3. AMR_02_ResFinder.r summarizes results for isolates from ResFinder.<br />
4. AMR_03_functions.r and AMR_04_proteomics.r process and analyze labelled LC-MS/MS proteomics data from MaxQuant output of searches on the CARD (The Comprehensive Antibiotic Resistance Database) and SwissProt databases.<br />

### The four folders contain data processed by each tool. <br />
Both genomics and proteomics datasets on four Campylobacter jejuni isolates were previously published in Clart et al. PLoS One 2018 https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190836<br />
1. contigsFromSpades: Fastq files were downloaded from NCBI and reanalyzed with SPAdes 3.11.1 tool to obtain contigs with no minimum length cut-off. Contigs from the isolates are stored here.<br />
2. RGI: Genomic AMR detection results from RGI	v4.1.0 tool <br />
3. ResFinder: Genomic AMR detection results from ResFinder 3.1 tool <br />
4. MaxQuant: Results from reanalyzing proteomics data on MassIVE MSV000081410 https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=0efcd5aed02c4abc98e93571dbac2b6d<br />
    - ExperimentalSetup.txt: The corresponding iTRAQ labels of the isolates within each data<br />
    - proteinGroups.txt: Protein abundance reported by MaxQuant <br />
    - evidence.zip: Abundance at peptide level reported by MaxQuant <br />
    

## Legal ##


Copyright Government of Canada 2018

Written by: Julie Chih-yu Chen @ National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

## Contact ##


**Julie Chih-yu Chen**: chih-yu.chen@canada.ca
 
