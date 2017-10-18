#Reformat GEMMA outputs to be on the log scale. Filter files down to meta-analysis ready files. Meta analyse, make QQ and Manhattan plots, do forest plots, LDSC

#Set working directoy=ry
workingdir=$(pwd)

ls | grep pcs$ > study_files.txt

#Convert all data into the log scale
for file in $(cat study_files.txt)
do
 sname=$(echo $file | awk 'BEGIN{FS="_"}{print $1}')
 
 #Phenotype files must have the same name as the study!
 phenofile="$sname".gemma.pheno
#Take average value of phenotype, which is the prevalence
 #Assumes that phenotype is 1/2 coded!
 phi=$(grep -v NA $phenofile | awk '{print $1 - 1}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')

 qsub -lwalltime=00:25:00 beta_to_logscale.pbs -e errandout/ -o errandout -d $workingdir -F "-s beta_to_logscale.rscript -i $file -v $phi -b beta -m af -p p_wald -o "$sname".logscale"
 
 #Get extra annotation info (cases/controls counts)
 ncon=$(grep 1 $phenofile | wc -l | awk '{print $1}')
 ncase=$(grep 2 $phenofile | wc -l | awk '{print $1}')
  
 echo $sname $ncase $ncon >> case_control.counts
 
done


beta_to_logscale.pbs -s beta_to_logscale.rscript -i $file -p $phi -b beta -m af -p p_wald -o "$sname".logscale

#Make the sequence of Manhattan plots (all gemma studies, then add in qimr, then add in all external ones)

Rscript make_mi_files.R

#Leave one out analysis (or run multiple data) 

 ls metal_scripts | grep LOOM  > metal_jobs12.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

 cat metal_jobs114.txt metal_jobs115.txt | grep -v lat > metal_jobs116.txt

 ncommands=$(wc -l metal_jobs14.txt | awk '{print $1}')
 nodesize=16
 nodeuse=$(($nodesize ))
 totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1-$totjobs -l walltime=00:05:00 run_meta_v2_loo.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m metal_jobs14.txt -p $metalpath -n $nodeuse"

 cat results/GTMALeur_39_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/GTMALeur_39.results
 cat results/GTFEMeur_47_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/GTFEMeur_47.results

 ls results_cat | grep -E eur | grep SOBP > metal_outputsA.txt

 ncommands=$(wc -l metal_outputsD.txt | awk '{print $1}')
 nodesize=16
 nodeuse=$(($nodesize ))
 totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1-$totjobs -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsA.txt -c green -p 0.05"
 
 
 awk '{print $6}'  results_cat/lat_38.results | grep -v NA > results_cat/lat_38.results.p

 ls results_cat | grep .p$ > qq_plot.files
 qsub -l walltime=00:20:00 qq_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-q qq_plot.R -e 1 -p qq_plot.files "

