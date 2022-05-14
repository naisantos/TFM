#Introduciendo los datos de los archivos creados con Ens2Gb.R, llevamos a cabo los ánalisis de Go y KEGG, utilizando dos scripts que ya estaban hechos

###Script to run GO from DiffBind Results: DiffBindGO.sh
bash /media/sequentia/sharing/src/012-GOEA+statistics_script/Script_GOEA_Alfonso/goea-script.sh --FC_TSV /media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/Male_FC.tsv --g /media/sequentia/synology_office/Naiara/Referencia/go2genome_mouse.txt --output /media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/GOEA/
  

  
  ###Script to run kegg from DiffBind Results: DiffBindkegg.sh
bash /media/sequentia/sharing/src/015-KEGG_enrichment/KEGG-ea-script.sh --FC_TSV /media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/Male_FC.tsv --g /media/sequentia/synology_office/Naiara/Referencia/genome2kegg.tsv --ens2gb /media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/Male_entrez.tsv --specie mmu --GeneID kegg --output /media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/KEGG/
  
  
  










