
# GuitarPlotFast.R 
Draw a metagene plot based on PARE/Degradome dataset.

## Installation

Download GuitarPlotFast.R in your $HOME/bin direction which should be added to $PATH variable
```
wget -O ~/bin/GuitarPlotFast.R https://raw.githubusercontent.com/LabHMChenABRC/PARE-analysis/main/Metagene-plot/GuitarPlotFast.R
```
Make GuitarPlotFast.R executable 
```
chmod +x ~/bin/GuitarPlotFast.R
```
## Dependency
R packages: Guitar, data.table, R.utils, ggplot2 and cowplot

## Usage
### Make density files:
``` shell
GuitarPlotFast.R -b <Directory_of_bed.gz> -g <gtf> -o <Output_folder>

-b <Directory_to_bed> contains "bed.gz" files created by bedtools based on genomic bam files which are produced by STAR.
                      bedtools bamtobed -i <filename.bam/cram> | bgzip -@ <number of core  >directory_to_bed/filename.bed.gz
-g <gtf>              gtf file, using gffread can convert gff to gtf format:
                      gffread -T -o <.gtf> <.gff>
-o <Output_folder>    program will save .mrna.density files in this direction
```

 This command will extract coding genes with 5'UTR, CDS, 3'UTR, upstream sequence, and downstream sequence lengths of at least 100bp each. And count the 5'Ends of unique alignments (MAPQ=255) on the forward strand of selected genes. The density distribution of the upstream sequence, 5'UTR, CDS, 3'UTR, and downstream sequence is scaled to 1:2:4:2:1 ratio.  
 * Calculate the density of 20M library for less than 15 min.  
 * STAR assigns unique alignments with MAPQ 255.  
 * To reduce RAM requirements and improve performance, some functions of R package Guitar are modified, and their name tail with .fast.
  
### Create metagene plot:
``` shell
GuitarPlotFast.R -m -p <output_file_prefix> -d <metagenePlot_info_file> -o <output_folder>

-m                          switch to plot metagene figure
-p <output_file_prefix>     output file prefix
-d <metagenePlot_info_file> a tabular file with a header of 'file' and 'plotgroup' provide this:
                            <Path of density file> <Group1[,Group2,...]>
-o <output_folder>          output direction.
```
`GuitarPlotFast.R -m` will produces <prefix>.mRNA-metaplot.pdf in <output_folder>
* The below example of <metagenePlot_info_file> will used to create a plot with group1 (WT, xrn4-6) and group2 (WT, fry1-6).

  | file                                 | plotgroup |
  | ------------------------------------ | --------- |
  | Path-of-density/WT.mrna.density      | 1,2       |
  | Path-of-density/xrn4-6.mrna.density  | 1         |
  | Path-of-density/fry1-6.mrna.density  | 2         |

<img src="https://github.com/LabHMChenABRC/PARE-analysis/assets/58059039/71238714-8b56-4bb6-93e4-d19d4ac1e0b6" width=50% height=50%>

