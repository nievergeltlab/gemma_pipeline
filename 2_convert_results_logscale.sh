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
 #head $phenofile
 #echo $sname $ncase $ncon >> case_control.counts
 
done


ls *.logscale > logscale_studies.txt
echo "pris.logscale
safr.logscale" > logscale2.txt

#If you missed a study, create a list like this:
#echo daner_comc_aam_analysis_run3.gz  daner_meg2_af3_analysis_run5.gz | sed 's/ /\n/g' > studylist.txt

njobs=$(wc -l logscale_studies.txt | awk '{print $1}')

qsub -t1-$njobs -l walltime=00:05:00 split_files_gemma.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f logscale_studies.txt -s  $workingdir -o metal_inputs -m 0.01 -i 0.6"
qsub -t1-2 -l walltime=00:10:00 split_files_gemma.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f logscale2.txt -s  $workingdir -o metal_inputs -m 0.01 -i 0.6"

#have to clean up queensland, dont know exactly what the format difference is with the other data
for study in qimr.logscale # $(cat studylist_related.txt )
do
 for chr in {1..22}
 do
  awk -v chr=$chr -v maf=$maf  '{if  (NR == 1) print "SNP","A1","A2","MAF","OR","SE","P" ;  if (($2 ==chr) && ($7 > maf) && ($7 < 1-maf) && ($15 > -5) && ($15 < 5))  print $1,$5,$6,$7,exp($15),$16,$12}'  $study > metal_inputs/"$study"_"$chr"
 done
# rm -f  metal_inputs/"$study".unzip
done



#Split the non gemma data
njobs=$(wc -l externalsites.studies | awk '{print $1}')

qsub -t1-$njobs -l walltime=00:10:00 split_files.qsub -d $workingdir -e errandout/ -o errandout/ -F "-f externalsites.studies -s  $workingdir -o metal_inputs -m 0.01 -i 0.6"



#Make the sequence of Manhattan plots (all gemma studies, then add in qimr, then add in all external ones)

Rscript make_mi_files.R

#Leave one out analysis (or run multiple data) 

metalpath=/home/cnieverg/gemma_gwas/results/generic-metal/metal

 ls metal_scripts > metal_jobs.txt #Put in a grep for the last row of the study input data if you ONLY want the total data meta analysis

 ncommands=$(wc -l metal_jobs.txt | awk '{print $1}')
 nodesize=16
 nodeuse=$(($nodesize ))
 totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1-$totjobs -l walltime=00:15:00 run_meta_v2_loo_v2.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m metal_jobs.txt -p $metalpath -n $nodeuse"

 
 
 cat results/all_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/all.results
 cat results/ext_*.tbl | awk '{if (NR == 1 || $1 != "MarkerName") print }' > results_cat/ext.results
 
 grep -v ?????? results_cat/all.results > results_cat/all.results.fix
 grep -v ?????? results_cat/ext.results > results_cat/ext.results.fix

 #Files for correlation between results
 awk '{print $1,$6}' results_cat/all.results.fix > results_cat/all.results.fix.snpp
 awk '{print $1,$6}' results_cat/ext.results.fix > results_cat/ext.results.fix.snpp
  

 ls results_cat | grep fix  > metal_outputsA.txt

 ncommands=$(wc -l metal_outputsA.txt | awk '{print $1}')
 nodesize=16
 nodeuse=$(($nodesize ))
 totjobs=$(( ($ncommands + $nodeuse - 1 ) / $nodeuse ))

 qsub -t1 -l walltime=00:20:00 mh_plot_rapid.qsub -d $workingdir -e errandout/ -o errandout/ -F "-m mh_plot_pgc.R -s ManhattanPlotterFunction_colorfixed_max10ylim2_pgc.R -l phase3.locations2 -d metal_outputsA.txt -c blue -p 0.05"
 
 
 awk '{print $6}' results_cat/all.results.fix | grep -v NA > results_cat/all.results.fix.p
 awk '{print $6}' results_cat/ext.results.fix | grep -v NA > results_cat/ext.results.fix.p
 

 ls results_cat | grep .p$ > qq_plot.files
 qsub -l walltime=00:20:00 qq_plot.qsub -d $workingdir -e errandout/ -o errandout/ -F "-q qq_plot.R -e 1 -p qq_plot.files "

 #LDSC prepare
 awk '{if (NR==1) print "SNP","A1","A2","BETA","SE","P","N_Case","N_Control"; else print $1,$2,$3,$4,$5,$6,"21667","58465"}' results_cat/ext.results.fix > results_cat/ext.results.fix.ldsc
 awk '{if (NR==1) print "SNP","A1","A2","BETA","SE","P","N_Case","N_Control"; else print $1,$2,$3,$4,$5,$6,"18990","52417"}' results_cat/all.results.fix > results_cat/all.results.fix.ldsc
 cd results_cat
 zip ext.results.fix.ldsc.zip ext.results.fix.ldsc
 zip all.results.fix.ldsc.zip all.results.fix.ldsc
 
awk '{print $1}' /home/cnieverg/gemma_gwas/results/results_cat/ext.results.fix.ldsc | LC_ALL=C sort -k1b,1  > /home/cnieverg/gemma_gwas/results/results_cat/ext.results.fix.ldsc.snplist
  
 munge_sumstats.py --sumstats results_cat/ext.results.fix.ldsc --out ext.results.fix.ldsc.munge

ldsc.py \
--h2  ext.results.fix.ldsc.munge.sumstats.gz \
--ref-ld-chr /home/cnieverg/report_data/testrun/results_cat/eur_w_ld_chr/ \
--w-ld-chr  /home/cnieverg/report_data/testrun/results_cat/eur_w_ld_chr/ \
--out ext.results.fix.ldsc

