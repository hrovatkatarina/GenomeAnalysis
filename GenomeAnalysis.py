import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from fpdf import FPDF
import sys

if len(sys.argv) != 2:
    print("Usage: python GenomeAnalysis.py <filename.bam>")
    sys.exit()

filename = sys.argv[1]

samfile = pysam.AlignmentFile(filename, "rb")

# TASK 1: Average coverage across genome

# Retrieve a list oh chromosomes from the header (First 22 chromosomes and X and Y)
chromosome = []
chromosomes_header = samfile.header['SQ']
for ch in chromosomes_header:
    chromosome.append(ch['SN'])
    chromosome = chromosome[:24]

# Calculate the average coverage for each chromosome in the genome
genome_coverage = []
for chrom in chromosome:
    coverage = samfile.count_coverage(contig=chrom)
    arr_coverage = np.array(coverage)
    tot_coverage = arr_coverage.sum(axis=0)
    average_coverage = np.mean(tot_coverage)
    genome_coverage.append(average_coverage)

# Plot average coverage across genome
plt.figure(figsize=(10, 6))
plt.plot(chromosome, genome_coverage, marker='o')
plt.title('Average coverage across the genome')
plt.xlabel('Chromosome')
plt.ylabel('Coverage')
plt.savefig("Coverage.png")

samfile.close()

# TASK 2 - Total number of reads

# Total number of reads (including those which were not mapped to human chromosomes)
samfile = pysam.AlignmentFile(filename, "rb")
num_read = 0
for read in samfile:
    num_read += 1
#print("Total number of reads:", num_read)

samfile.close()

# The total count of reads that were mapped to human chromosomes
samfile = pysam.AlignmentFile(filename, "rb")
total_reads = 0
for chrom in chromosome:
    nm_reads = samfile.count(contig=chrom)
    total_reads += nm_reads
#print("Total number of reads that were mapped to human chromosomes:", total_reads)
samfile.close()

# TASK 3: Number of reads for each chromosome
samfile = pysam.AlignmentFile(filename, "rb")
reads_list = []

for chrom in chromosome:
    nm_reads = samfile.count(chrom)
    reads_list.append((chrom, nm_reads))
    #print(f"Number of reads for chromosome {chrom}: {nm_reads}")

# Store number of reads per chromosome in dataframe
df = pd.DataFrame(reads_list, columns=['Chromosome', 'Number of reads'])
samfile.close()

# TASK 4: GC content across entire genome
samfile = pysam.AlignmentFile("alignment-bam.bam", "rb")
total_gc = 0
total_b = 0

for read in samfile:
    sequence = read.query_sequence
    gc_count = sequence.count('G') + sequence.count('C')
    total_gc += gc_count
    total_b += len(sequence)

gc_perc = (total_gc / total_b) * 100
gc_percentage = round(gc_perc, 2)
#print("GC percentage of all reads:", gc_percentage)

samfile.close()

# Create PDF report
pdf = FPDF()
pdf.add_page()
pdf.set_font("Helvetica", style = 'B', size = 13)
pdf.cell(200, 10, txt = "Report: Genetic data Analysis", ln = 1, align = 'C')
pdf.set_font("Helvetica", size = 11)
pdf.cell(200, 10, txt = "1. Average coverage across the genome", ln = 1, align = 'L')
pdf.image("Coverage.png", w=200, h=120)
pdf.cell(200, 10, txt = "2.1 Total number of reads: " + str(num_read), ln = 1, align = 'L')
pdf.cell(200, 10, txt = "2.2 Total number of reads that were mapped to human chromosomes: " + str(total_reads), ln = 1, align = 'L')
pdf.cell(200, 10, txt = "3. Number of reads for each chromosome: ", ln = 1, align = 'L')
pdf.cell(200, 10, txt = "Chromosome: Number of reads", ln = 1, align = 'L')
for index, row in df.iterrows():
    pdf.cell(200,5, txt = "{}: {}".format(row["Chromosome"],row["Number of reads"]), ln = 1, align = 'L')
pdf.cell(200, 10, txt = "4. GC percentage across the entire genome : " + str(gc_percentage) + "%", ln = 1, align = 'L')
pdf.output("Report.pdf")
print("Report created.")
