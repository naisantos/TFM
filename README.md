# TFM

En este documento aparecen los archivos de los diferentes programas y un breve resumen de cada uno de ellos. El orden de los programas es el mismo orden que el de su utilización durante el trabajo. 



* H3K4me3_nRC_download.sh -> Para descargar los archivos
* Counter.sh -> Para mirar si coinciden el número de lecturas con las teóricas que nos mandan. En el archivo aparecen los comandos a utilizar tanto con los datos originales como con los recortados.
* Trimming.sh -> Para poder hacer el recorte de los datos
* fastqc.sh -> Llevar a cabo el control de calidad de los archivos
* Bowtie2_build.sh -> Para hacer el mapreo de los datos
* Filtering.sh -> Filtrar los archivos BAM para quedarnos con los mapeos únicos
* StatisticsMapping.sh -> Para obtener estadísticas básicas sobre los datos mapeados, entre ellos el número de lecturas
* BigWig.sh -> Normalizar los datos y cambiarlos de formato a BigWig
* PeakCalling.sh -> Identificar los dos picos de enriquecimiento: picos estrechos y amplios
* RemoveBlacklist.sh -> Remover los duplicados de regiones genómicas espefíficas
* HomerMergedPeaks_BiologicalReplicated_v2.py -> Quedarnos con los picos solapantes entre las tres repicas y sus respectivos diagramas de Venn
* HomerMergedPeaks2Condition.py -> Quedarnos con los picos solapantes entre control y tratado y sus respectivos diagramas de Venn
* VennDaigram_HomerMerge.R -> El código necesario para hacer los gráficos de los diagramas de Venn
* ChIPseq_DiffBind.R -> Quedarnos con las coordenadas de los picos comunes después de haber hecho el análisis diferencial 
* PeakAnnotation.R -> Quedarnos solamente con los genes que no estén en común, ya sean en control o en tratado 
* Filter.R -> Filtrar los picos que van unidos a regiones reguladoras. En este script también aparece cómo lograr los diagramas de Venn y obtener varios subconjuntos para poder visualizar los datos en Enrichr y en IGV
* Ens2Gb.R -> Utilizando la base de datos de los picos comunes, crear 2 archivos para hacer el análisis de GO y KEGG
* Run GO, KEGG from DiffBind Result.sh -> Introduciendo los datos de los 2 nuevos archivos creados, obtenemos gráficos que resumen las funciones biológicas y las metabólicas que regulan los genes
* GOEA.R -> En uno de los archivos de GO y KEGG obtenidos gracias al anterior script, para que al lado del simibolo de los genes apareciesen también sus nombres
* symbolgenes.oy -> Hacer librerias para acceder más fácil a los nombres de los genes
* RNAChIP.R -> Integrar los datos de ChIP con RNA-seq, e intersecar los genes que se co-regulan y quedarnos con ellos. Además, creamos dos nuevos archivos para hacer el análisis de GO y KEGG. 



Asimismo, también está subida la carpeta donde aparecen las funciones biológicas de los genes de los diferentes análisis. Entre ellos están la ChIP de adultos y neonatales, la RNA-seq de los neonatalaes, la intersección entre ChIP y RNA-seq de adultos y neonatales y superposición entre neonatales y adultos de los análisis de ChIP y RNA-seq. La carpeta donde están estos archivos se llama 'Ontologías génicas'. 




