# Paramecium_MNase
Analyzing the reads mapping of Paramecium MNase experiment

bowtie2-build ()

bowtie2 -p (#threads number) -x (#index file produced by bowtie2-build command) -1 (#reads file 1) -2 (#reads file 2)  -S (#output SAM file)

samtools view -bF 4 | samtools sort - | samtools mpileup - > OUTPUTFILE (#samtools is used for dealing with sam file from bowtie2)

sam2pro

analysis profile files
