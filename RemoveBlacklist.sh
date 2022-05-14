# ##############################################
# REMOVE REGIONS IN THE BLACKLIST FROM PEAKS
# ##############################################

#Se denomina Blacklist al archivo que contiene coordenadas que deben ser filtradas, por lo cual hay que eliminar los picos que caen en esas regiones. En nuestro caso, el archivo de Blacklist se nos daba ya de pasados experimentos. 

# nohup bash ChIPSeq_BlackList.sh > ChIPSeq_BlackList.out 2>&1

#Usage:   bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>
#-a y -b son los files que se intersecan
#-v      Only report those entries in A that have _no overlaps_ with B. - Similar to "grep -v" (an homage).

cd /media/sequentia/synology_office/Naiara/H3K4me3/bam_files/Filtered/peaks/


#remove blacklist - narrow (H3K4me3)
ls /media/sequentia/synology_office/Naiara/H3K4me3/bam_files/Filtered/peaks/*_narrow_peaks.narrowPeak  | egrep -v "input" | while read line; do
    name=$(basename $line | sed 's/.fastq.gz//');
    echo "Procesing sample: $name";
    bedtools intersect -a ${name} -b /media/sequentia/synology_office/projects/482-Salva/references/mus_musculus_release104/blacklist/mm39-blacklist.v2.sorted.bed -v > ./blacklist/${name}.noBlackLists.narrowPeak -nonamecheck
done




