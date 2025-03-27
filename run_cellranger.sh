# Trong-Hieu Nguyen, last modified on 24.08.2022.

current_dir=$(pwd)

path_to_cellranger="/home/hieunguyen/CRC1382/src_2023/CellRanger/cellranger-7.1.0"
export PATH=${path_to_cellranger}:$PATH



##############################################################################
# M A I N - C E L L R A N G E R - C O M M A N D S
for sampleid in K1 K2 P1 PS;do \
    cellranger multi --id=${sampleid} --csv=multi_config_${sampleid}.csv --localcores=40;
    done
# EOF


