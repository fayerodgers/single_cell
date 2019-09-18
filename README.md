# Single cell

## Retrieving data


Need a manifest with Sanger sample IDs and meta data

```
cut -f 1 samples.txt | while read -r sample; do
 CRAMS=($(imeta qu -z seq -d sample = ${sample} and target = 1 and type = cram | grep -o [0-9_#]*.cram))
 for cram in ${CRAMS[@]}; do
  echo ${sample}$'\t'${cram}$'\t''CRAM' >> data_locations.txt
 done
 COLLECTIONS=($(imeta qu -z seq -d sample = ${sample} and target = 1 and type = cram | grep 'collection' | sort | uniq | grep -o '\/seq.*'))
#check CRAMs were all sequenced in the same flow cell (all CRAMS are in the same collection)
 if [[ "${#COLLECTIONS[@]}" = 1 ]] ; then
  #get CellRanger reports by sample ID
  PATHS=($(ils ${COLLECTIONS[0]}/cellranger | grep ${sample}| grep -o '\/seq.*'))
  #!!!! NB 2 samples are named differently- need to add paths manually for first CR run of 4672STDY6814755 and 4672STDY6814756!!! 
 #CellRanger reports for samples sequenced in >1 flow cell end up in /seq/illumina 
 else
  PATHS=($(ils /seq/illumina/cellranger | grep ${sample}| grep -o '\/seq.*'))
 fi
 for path in ${PATHS[@]}; do 
  echo ${sample}$'\t'${path}$'\t''CELLRANGER' >> data_locations.txt
 done
done
```

CRAMs are named like "flowcell_lane#index.cram". Can later reconstruct paths on irods like this: /seq/flowcell/flowcell_lane#index.cram.

Retrieve all available cellranger directories for these samples:

