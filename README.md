# Disease-Assoicated Microbiome base Blood Screening in python (DAMBSpy)
![linux-64](https://img.shields.io/badge/Linux-FCC624?style=flat&logo=linux&logoColor=black)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/)
![snakemake](https://img.shields.io/badge/snakemake-7.32.4-green.svg)
[![Documentation Status](https://readthedocs.org/projects/abstar/badge/?version=latest)](https://abstar.readthedocs.io/en/latest/?badge=latest)
![](https://img.shields.io/badge/license-MIT-blue.svg)


# ✨ Introduction
Recent studies have demonstrated that microbial DNA, which is typically absent from the blood, can be detected in patients with non-bacterial infectious diseases such as cancer. These findings indicate a potential association between the presence of specific microbes and certain diseases ([_Poore_ et al., Nature 2020](https://doi.org/10.1038/s41586-020-2095-1) ). The study revealed that low levels of foreign microbial DNA in the blood could serve as predictors of disease in patients. This ability to predict diseases based on the presence of specific microbes suggests that conditions might be diagnosed through simple liquid biopsies, circumventing the need for invasive tissue biopsies.

To detect these low levels of DNA in blood, several advanced technologies have been employed. Among these, droplet digital PCR (ddPCR) is particularly notable for its superior sensitivity compared to other methods. ddPCR's high sensitivity is attributed to its approach of partitioning the sample into thousands of droplets, which allows for the amplification and detection of individual DNA molecules within each droplet. This partitioning minimizes the impact of inhibitors and reduces background noise, resulting in highly accurate and quantitative detection of low-abundance DNA sequences. Consequently, ddPCR is exceptionally well-suited for identifying rare microbial DNA in blood samples, making it a valuable tool for early disease diagnosis.

This project aims to standardize the analysis pipeline using Python and various open-source tools. Additionally, it includes the design of ddPCR primers to verify the detection of blood microbes. Ultimately, the goal is to establish a foundational pipeline for disease diagnosis based on blood microbial markers.


# 🚩 Key Features
- **Blood WGS Data Analysis**: This component involves the comprehensive analysis of whole genome sequencing (WGS) data obtained from patient blood samples. The process begins by mapping the sequencing data to the human genome to identify sequences that align with human DNA. Subsequently, the unmapped reads—those sequences that do not match the human genome—are analyzed. These unmapped reads are scrutinized to detect the presence of non-human microbial DNA, which could be indicative of microbial infections or associations with specific diseases. This step is crucial for identifying rare microbial DNA signatures that may be linked to disease states.

- **Automatic ddPCR Primer Design for Bacterial Targets**: This component focuses on enhancing the sensitivity of microbial DNA detection through the design of droplet digital PCR (ddPCR) primers. Once the microbial DNA has been identified from the WGS data, specific primers for ddPCR are designed automatically. These primers are tailored to target the microbial DNA sequences detected in the blood samples. The use of ddPCR allows for highly sensitive and precise quantification of these microbial DNA sequences, facilitating accurate and reliable detection even at low abundance levels. This step ensures that the detection method is both robust and scalable for diagnostic applications.

- **Standardization of Analysis Pipeline**: This tool aims to improve and standardize the entire analysis pipeline. The goal is to develop a streamlined, efficient, and reproducible workflow that can be consistently applied across different samples and settings. This involves refining existing methods to eliminate inefficiencies and complexities, thereby making the process more accessible and user-friendly. By standardizing the pipeline, the project ensures that the analysis can be reliably replicated and that results are consistent and comparable across different studies and applications. This standardization is essential for translating the pipeline into a practical diagnostic tool for clinical use, enabling the routine detection of microbial markers in blood for disease diagnosis.


# ✅ Installation
**DAMBSpy is only available for linux. Windows user can still execute DAMBSpy with a Windows Subsystem for Linux (WSL) such as Ubuntu.**

Installing DAMBSpy currently **REQUIRES 1) conda AND 2) pip.** Make sure both are installed prior to following the installation instructions below. 

## Exported environment (x86-64 & linux-64)
You can set up the conda environment to run DAMBSpy using the YAML definition found in this repository:
<https://haje01.github.io/2020/04/21/snakemake-tutorial.html>
```bash
# Clone DAMBSpy GitHub directory
git clone https://github.com/Tropini-lab/DAMBSpy.git

# Change directory
cd DAMBSpy

# Create and set up conda environment
conda deactivate
conda env create -f dambs_env.yml
conda activate dambs
```

# ⚠️ Important: Before You Start !!
- Processing WGS data is **time-consuming** and requires **substantial RAM**. Additionally, the **huge storage size** needed for bacterial DNA analysis is significant.
    - [Hardware specifications for human data processing](https://github.com/kaist-ina/BWA-MEME?tab=readme-ov-file#notes)
    - [Hardware specifications for bacterial data processing](https://ccb.jhu.edu/software/kraken/MANUAL.html#system-requirements)
- Bacterial DNA is present in very low quantities and **may not be detectable in low-depth WGS data**.
- PCR can be unpredictable, and while primers may appear perfect in silico, **it is strongly recommended to confirm their specificity in vitro before use**.

# 🏷️ How it works
DAMBSpy takes any number of bacterial CDS files as input. ~~~~~~~~~~~~~~~~~



# ✒️ Usage
**Always activate your conda environment before running DAMBSpy.**

```sh
conda activate dambs
```
**Check that you are using the most recent PUPpy version and update if needed**

```sh
conda list dambs # Check version
conda update dambs # Update puppy package
```

DAMBSpy operates in 2 main steps:

1) ``dambs-align`` - performs local pairwise sequence alignment of all the input genes against each other, and
2) ``dambs-primers`` - designs taxon-specific primers based on user-determined parameters.

**IMPORTANT:** The syntax to run PUPpy varies depending on the installation mode.

If installed with **conda**, you do not need to specify the path to the scripts. For example:

```sh
dambs-align -h
```

PUPpy can be executed either from **command-line**.

## Command-line execution

Detailed usage information, including all the primer design parameters, can be seen by running ``-h`` or ``--help`` at each step.

```
dambs-align -h
dambs-primers -h
```

Remember to specify the **path** to the script if you installed PUPpy from the exported environment! 

### 1. Gene alignment (dambs-align)

The alignment step must always be run first for any **new** defined bacterial community.

```
puppy-align -pr <PATH>/test/INPUT_primerTarget -nt <PATH>/test/INPUT_nonTarget -o <PATH>/test/OUTPUT_puppy-align
```

This command creates the output file ``<PATH>/test/OUTPUT_puppy-align/ResultDB.tsv`` which can be used as input for the primer design command (step 2). The command `puppy-primers` can be run as many times as desired without having to rerun `puppy-align` again, as long as the bacterial community remains unchanged.

### 2. Primer design (dambs-primers)

The second step consists in designing taxon-specific primers unique to individual members or shared by groups in the bacterial community.

``puppy-primers`` **requires** 2 arguments as input:



# ⚙️ DAMBSpy parameters
Parameters for **``dambs-bacteria``**
```
usage: dambs-bacteria [-h] <채우기>

```

Parameters for **``dambs-primers``**
```
usage: dambs-primers [-h] <채우기>

```


# 📓 Reference
- Poore, G. D., et al. (2020). Microbiome analyses of blood and tissues suggest cancer diagnostic approach. Nature.
- Zmrzljak, U. P., et al. (2021). Detection of Somatic Mutations with ddPCR from Liquid Biopsy of Colorectal Cancer Patients. Genes.
- Ghezzi, H., et al. (2023). PUPpy: a primer design pipeline for substrain-level microbial detection and absolute quantification. doi:http://dx.doi.org/10.14288/1.0438913