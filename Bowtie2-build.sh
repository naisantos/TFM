##Se hace el mappeo de los datos trimmeados. Para ello, seguimos 4 pasos necesarios que numeraremos y explicaremos, y despues los juntaremos todos en un loop para ejecutarlo con nuestros archivos. 

##1- Hacer un indice del genoma. Para ello, utilizamos dos comandos que son bowtie2-build y bowtie2-inspect. 
#bowtie2-build builds a Bowtie index from a set of DNA sequences. bowtie2-build outputs a set of 6 files with suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2. 
#The bowtie2-build indexer

nohup bowtie2-build GRCm38.primary_assembly.genome.fa GRCm38.primary_assembly.genome.idx & 

## bowtie2-inspect extrae información de un índice Bowtie sobre qué tipo de índice es y qué secuencias de referencia se utilizaron para construirlo.

bowtie2-inspect GRCm38.primary_assembly.genome.idx &

##2 - Alineando lecturas al genoma con Bowtie2
#Usage:
  ##bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i>} [-S <sam>]

  #-x: /path/to/genome_indices_directory
  #-U: /path/to/FASTQ_file
  #-S: /path/to/output/SAM_file


bowtie2 -p 10 -q --local \
-x ~/Naiara/Referencia/GRCm38.primary_assembly.genome.idx \
-U ~/Naiara/H3K4me3/Trimmed_output/Trimmed_files/CF1_k9ac_01286AAD_GATCAG_R1_001.fastq.gz \
-S ~/Naiara/H3K4me3/Bowtie2/CF1_k9ac_01286AAD_GATCAG_R1_001.sam &

## 3- Pasar de formato SAM a BAM, que es un formato binario 
#The command we will use is samtools view with the following parameters:

#-h: include header in output
#-S: input is in SAM format
#-b: output BAM format
#-o: /path/to/output/file

 samtools view -h -S -b \
-o CF1_k9ac_01286AAD_GATCAG_R1_001.sam_unsorted.bam \
CF1_k9ac_01286AAD_GATCAG_R1_001.sam_unsorted.sam

##4 - Ordenar los archivos BAM usando las coordenadas del genoma
#The command we will use is sambamba sort with the following parameters:

#-t: number of threads / cores
#-o: /path/to/output/file

 sambamba sort -t 2 \
-o CF1_k9ac_01286AAD_GATCAG_R1_001.sam_sorted.bam \
CF1_k9ac_01286AAD_GATCAG_R1_001.sam_unsorted.bam

#################################################
#Un loop para todos los archivos 
#################################################

cd /media/sequentia/synology_office/Naiara/ 

for f in /media/sequentia/synology_office/Naiara/H3K4me3/Trimmed_output/Trimmed_files/*.tr.fq.gz 
do
    name=$(basename $f | sed 's/.tr.fq.gz//'); 
    echo "Procesing sample: ${name}" ;  ç

    bowtie2 -p 10 -q --local \ -x ./Referencia/GRCm38.primary_assembly.genome.idx \ -U ./H3K4me3/Trimmed_output/Trimmed_files/${name}.tr.fq.gz | samtools view -Sb - > ${name}.bam; 
    sambamba sort -o ${name}.sort.bam ${name}.bam

     
done


# Este script toma un archivo fastq de datos ChIP-seq, ejecuta FastQC y produce un archivo BAM listo para la llamada de picos. Bowtie2 es el alineador utilizado, y el archivo BAM resultante se ordena por coordenadas genómicas y se eliminan las lecturas duplicadas utilizando sambamba.
