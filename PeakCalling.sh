# ##############################################
# PEAK CALLING
# ##############################################


#El peak calling se utiliza para ver que genes estan activos. Para saber que argumentos utilizar con las funciones, hay que mirar en articulos cuales han utilizado y que resultados les han dado, para asi poder compararlos con los nuestros

cd /media/sequentia/synology_office/Naiara/H3K4me3/bam_files/Filtered/


#peak calling - (H3K4me3)
ls /media/sequentia/synology_office/Naiara/H3K4me3/rawdata/*fastq.gz | egrep -v "input" | 
while read line; do
    name=$(basename $line | sed 's/.fastq.gz//');
    echo "Procesing sample: $name";
    macs2 callpeak -t ${name}.sort.rmdup.MQ30.bam -g mm -q 0.05 -n ./peaks/${name}_narrow
    macs2 callpeak -t ${name}.sort.rmdup.MQ30.bam -g mm -m 5 50 -p 1e-5 --broad -n ./peaks/${name}_broad
    
done











