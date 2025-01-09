# 1000WGS
The main file for generation intervals is 00_intervel_according2fasta.dict.r
Delete all the file under 00intervals and rerun the R script , the interval range can be changed within the R context

split -l 1 -d -a 3  --additional-suffix=.bed ../windows.bed hg38interval.
