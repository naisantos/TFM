######################
#Filtering + Mapping
#####################

#Primero filtramos los archivos bam para quedarnos solo con los mappeos unicos y no repetidos
#PicardCommandLine MarkDuplicates -> Premove duplicates from BAM files with the following command -> Identifies duplicate reads.
# -i One or more input SAM or BAM files to analyze.
#-o The output file to write marked records to  Required.
# -M File to write duplication metrics to  Required.

#Despues, tomamos solamente los que tienne una calidad específica
#Tenemos que coger datos en los que la largura de estos supere a 30, ya que normalmente estan entre 35-50 -> estos valores los logramos de los files de TrimmedQC
#-q INT   only include reads with mapping quality >= INT [0]

#Finalmente, volvemos a hacer un index, como con el genoma


#########
#LOOP
#########

cd /media/sequentia/synology_office/Naiara/H3K4me3/bam_files

for f in /media/sequentia/synology_office/Naiara/H3K4me3/rawdata/*fastq.gz

do
    name=$(basename $f | sed 's/.fastq.gz//');
    echo "Procesing sample: $name"; 
    
    PicardCommandLine MarkDuplicates I=${name}.sort.bam O=./Filtered/${name}.sort.rmdup.bam M=${name}.sort.rmdup.txt
    samtools view -q 30 -b ${name}.sort.rmdup.bam > ${name}.sort.rmdup.MQ30.bam 
    samtools index ${name}.sort.rmdup.MQ30.bam 
done





