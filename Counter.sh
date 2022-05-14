#Hay que mirar si coinciden el numero de lecturas con las teóricos que nos mandan. 

# ##############################################
# Count
# ##############################################

cd /media/sequentia/synology_office/Naiara/H3K4me3/

for f in /media/sequentia/synology_office/Naiara/H3K4me3/rawdata/*.fastq.gz
do
    name=$(basename $f | sed 's/.fastq.gz//');
    echo "Procesing sample: $name" ;
    zcat rawdata/$name.fastq.gz | wc -l ##El comando zcat se usa para contar cuantas lineas hay; al obtener el numero, hay que dividirlo entre 4 para ver el número real de secuencias y compararlas 
    # wc stands for word count 
    # -l, list compressed file contents
done

# ##############################################
# Count for Trimmed files
# ##############################################

cd /media/sequentia/synology_office/Naiara/H3K4me3/Trimmed_output

for f in /media/sequentia/synology_office/Naiara/H3K4me3/Trimmed_output/Trimmed_files/*.tr.fq.gz
do
    name=$(basename $f | sed 's/.tr.fq.gz//');
    echo "Procesing sample: $name" ;
    zcat Trimmed_files/$name.tr.fq.gz | wc -l 
done





