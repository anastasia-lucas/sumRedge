for file in $(ls simulated_data/2pop/*bed | sed 's/.bed//'); do ./plink --bfile $file --hwe 0.0000000001 --make-bed --out $file"-hwe" ; done

for file in $(ls simulated_data/2pop/*-hwe.bed | sed 's/.bed//'); do ./plink --bfile $file --recodeA --out $file ; done
