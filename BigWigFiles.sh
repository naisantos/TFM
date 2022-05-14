###########
#Create BigWig files
##########

#Para empezar con los peak calling, tenemos que normalizar los datos y cambiarlos de formato a BigWig. 


#This tool takes an alignment of reads or fragments as input (BAM file) and generates a coverage track (bigWig or bedGraph) as output. The coverage is calculated as the number of reads per bin, where bins are short consecutive counting windows of a defined size. It is possible to extended the length of the reads to better reflect the actual fragment length. *bamCoverage* offers normalization by scaling factor, Reads Per Kilobase per Million mapped reads (RPKM), counts per million (CPM), bins per million mapped reads (BPM) and 1x depth (reads per genome coverage, RPGC).

# --bam BAM file, -b BAM file
#-o FILENAME  -> Output file name. (default: None)
#  --normalizeUsing {RPKM,CPM,BPM,RPGC,None} -> Use one of the entered methods to normalize the number of reads per bin.
#-binSize INT bp, -bs INT bp -> Size of the bins, in bases, for the output of the bigwig/bedgraph file. (Default: 50)



cd /media/sequentia/synology_office/Naiara/H3K4me3/bam_files/Filtered

#create bigwig files
for f in /media/sequentia/synology_office/Naiara/H3K4me3/rawdata/*fastq.gz
do
    name=$(basename $f | sed 's/.fastq.gz//');
    echo "Procesing sample: $name";
    bamCoverage -p 10 --binSize 50 --normalizeUsing BPM -b ${name}.sort.rmdup.MQ30.bam -o ${name}.bpm.bw
done




df = pd.....(ensemble:symbols)

dictenssym = {}
for index,row in df.iterrows():
    dictenssym[row["ensembl"]] = row["symbol"]

    





