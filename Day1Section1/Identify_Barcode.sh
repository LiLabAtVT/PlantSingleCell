
#Run this command in your terminal. Make sure the gem_classification.csv file is located in your current working directory.
awk -F',' 'NR > 1 && $4 == "Arabidopsis_thaliana" { print $1 }' gem_classification_Pos24hpi_1.csv > arabidopsis_barcodes.txt
awk -F',' 'NR > 1 && $4 == "Pcap" { print $1 }' gem_classification_Pos24hpi_1.csv > pcap_barcodes.txt
awk -F','  'NR > 1 && $4 == "Multiplet" { print $1 }' gem_classification_Pos24hpi_1.csv > multiplet_barcodes.txt