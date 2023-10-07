# Genetic Data Analysis of BAM file

The primary outcomes are a line plot showcasing the average coverage across the genome and key statistics, including the number of total reads, the number of reads per chromosome, and the GC percentage. Script produces a report in PDF, containing the visualization and the calculated statistics.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)

## Installation

Clone the repository: 

```bash
git clone https://github.com/hrovatkatarina/GenomeAnalysis
```

List of Dependencies:
- pysam
- matplotlib
- pandas
- numpy
- fpdf

You can install them using the following command:

```bash
pip install pysam matplotlib pandas numpy fpdf
```

## Usage

To run the script, use the following command, where as bamfilename type name of your bam file:

```bash
python GenomeAnalysis.py bamfilename.bam
```
Note: You should have you bam file and index file associated with your bam file (.bam.bai) in the same directory.

## Output

- Report.pdf
- Coverage.png 

Outputs of the sample bam file are provided in Results folder.



