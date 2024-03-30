#!/usr/bin/env bash
wd=$PWD
density_info_file=$wd/Density.info.txt
GuitarPlotFast.R -m -p PARE -d $density_info_file -o $wd/PARE/metagenePlot
