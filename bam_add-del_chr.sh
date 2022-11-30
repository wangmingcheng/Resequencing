#! /usr/bin/bash

# 为bam文件染色体信息增加或者删除"chr"字符
#add chr
#samtools view -H $1 | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $1 > $2

#del chr
samtools view -H $1 | sed -e 's/SN:chr\([0-9XY]\)/SN:\1/' -e 's/SN:chrM/SN:MT/' | samtools reheader - $1 > $2

echo "Finished!"

