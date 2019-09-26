# Single cell

## Retrieving data


Need a manifest with Sanger sample IDs and meta data (samples.txt)

```
kinit #initialise irods

cut -f 1 samples.txt | while read -r sample; do

 CRAMS=($(imeta qu -z seq -d sample = ${sample} and target = 1 and type = cram | grep -o [0-9_#]*.cram))
 
 for cram in ${CRAMS[@]}; do
  echo ${sample}$'\t'${cram}$'\t''CRAM' >> data_locations.txt    #in case we need them later
 done
 
 #Retrieve CellRanger reports 
 #Check CRAMs were all sequenced in the same flow cell (all CRAMS are in the same collection)
 #CellRanger reports for samples sequenced in >1 flow cell end up in /seq/illumina 
 
 COLLECTIONS=($(imeta qu -z seq -d sample = ${sample} and target = 1 and type = cram | grep 'collection' | sort | uniq | grep -o '\/seq.*'))
 
 if [[ "${#COLLECTIONS[@]}" = 1 ]] ; then   
  PATHS=($(ils ${COLLECTIONS[0]}/cellranger | grep ${sample}| grep -o '\/seq.*'))
 else
  PATHS=($(ils /seq/illumina/cellranger | grep ${sample}| grep -o '\/seq.*')) 
 fi
 for path in ${PATHS[@]}; do 
  echo ${sample}$'\t'${path}$'\t''CELLRANGER' >> data_locations.txt
 done
 
done

#2 samples are named differently- need to add paths manually for first CellRangerv1 runs of 4672STDY6814755 and 4672STDY6814756
```

CRAMs are named like "flowcell_lane#index.cram". Can later reconstruct paths on irods like this, if needed: /seq/flowcell/flowcell_lane#index.cram.

Retrieve cellranger metrics for these samples. 

```
regex="cellranger$"
grep CELLRANGER data_locations.txt | cut -f 1,2 | while read -r sample path; do 

 version=$(echo $path | grep -o 'cellranger[0-9]*_count' | grep -o 'cellranger[0-9]*')
 annotation=$(echo $path | grep -o mm10\.*$)
 
 if [[ $version =~ $regex ]]; then
  version="cellranger131"
 fi
 
  #pull down the HTML files
 
 if [[ ! -d websummaries/${version} ]]; then   
  mkdir websummaries/${version}
 fi
 
 iget ${path}/web_summary.html websummaries/${version}/${sample}.html   
 
 #get the same data in CSV format
 
 if [[ ! -e ${version}.txt ]]; then      #because we want the headers first time round
  echo -n "sample_id,transcriptome," > ${version}.txt
  iget ${path}/metrics_summary.csv - | head -n 1 >> ${version}.txt 
 fi 
 
 echo -n $sample","$annotation"," >> ${version}.txt
 iget ${path}/metrics_summary.csv - | tail -n -1 >> ${version}.txt
 
 #get the counts matrices
 
 if [[ ! -d count_matrices/${version} ]]; then   
  mkdir count_matrices/${version}
 fi 
 
 counts_dir=$(ils -r ${path} | grep 'filtered' | grep 'C' | sort | tail -n 1 | sed -e 's/C-//' | sed -e 's/\s//g')

 iget -r ${counts_dir} count_matrices/${version}/${sample}
 
done
```

Metrics summaries are combined into a master file, cellranger_metrics.tsv.

```cellranger_metrics.R``` compares the samples that have been run with both cellranger v2.1.1 and cellranger v3.0.2. Plots in QC_plots.

## Seurat 

To run a standard analysis (eg):

```./basic_seurat.R --samples list.txt --metadata samples.txt --cellranger_version cellranger302 --results_directory dir```

First run it like this with the default options, then look at the plots to decide on filtering and clustering parameters.


