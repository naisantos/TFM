###############################
#Statistics for mapping
###############################

#Realiza un recorrido completo por el archivo de entrada para calcular e imprimir estadísticas en la salida estándar

cd /media/sequentia/synology_office/Naiara/H3K4me3/bam_files

for f in /media/sequentia/synology_office/Naiara/H3K4me3/rawdata/*fastq.gz
do
    name=$(basename $f | sed 's/.fastq.gz//');
    echo "Procesing sample: ${name}";
    samtools flagstat ${name}.sort.bam 
    
done



