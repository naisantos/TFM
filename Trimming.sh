#Hacer el recorte de los datos 

## python2.7 /scripts/trimming.py -h -> La ayuda sobre el comando. En mi caso, los argumentos utilizados y sus respectivas descripciones son las siguientes:


# -k == Minimum match against the adapters file to recognize part of a sequence as an adapter
# -g == Trimming algorithm to be used. Options: 'bbduk' and 'trimmomatic' [default bbduk]
# -- single == Input files are considered single-end
# -t == THREADS (ponemos 10 siempre)
# -c == --compression (por defecto es low)
# -o == Output directory to store the results
#-i DIRECTORY, --input DIRECTORY Directory where the input .fq/.fq.gz/.fastq/fastq.gz are in.



nohup python2.7 /scripts/trimming.py -i ./rawdata -k 11 -g trimmomatic --platform illumina --single -t 10 -c high  -o ./Trimmed_output &

