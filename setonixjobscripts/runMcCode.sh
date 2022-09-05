#!/bin/bash
#Script: run.sh
#Purpose: script for running a production version of the monte carlo code on
#         setonix. Sets up necessary directory structure and submits jobs to run
#         the code. Submits a final job to automatically collect output data and save it to
#         Acacia if all other jobs run successfully.
#Author: Reese Horton


#Create list of incident energies. 
#Fine grid creates around 275 jobs to submit
vals=$(seq 0.0 0.25 50.0)
morevals=$(seq 52.0 2.0 100.0)
evenmorevals=$(seq 110.0 10.0 120.0)
if [ -e tempens.txt ]; then
	 rm -rf tempens.txt
fi
echo "${vals[@]}" >> tempens.txt
echo "${morevals[@]}" >> tempens.txt
echo "${evenmorevals[@]}" >> tempens.txt
readarray -t enList < tempens.txt
rm -rf tempens.txt

#Set up required directories on scratch
if [ ! -d "$MYSCRATCH"/MCProduction ]; then
   mkdir "$MYSCRATCH"/MCProduction
fi

cd "$MYSOFTWARE"/MCSimulation
cp ./main "$MYSCRATCH"/MCProduction
cp ./data.in "$MYSCRATCH"/MCProduction
cd "$MYSCRATCH"/MCProduction

#Write energies to a file for data extraction scripts to use
rm energies.txt
for en in "${enList[@]}" 
do
   echo "$en" >> energies.txt
done

#Change the incident energy in the input file 'data.in' to a placeholder used in 
#in the submission scripts 'yyy'
lineNum=`(grep -n 'Incident particle energy (eV)' data.in | cut -c1-1)`
currEn=`(grep -n 'Incident particle energy (eV)' data.in | cut -c3-7)`
sed -i "${lineNum}s/${currEn}/yyy/g" data.in 
#Input file called "input" in submission scripts, rename accordingly
mv data.in input

if [ ! -d "$MYSCRATCH"/MCProduction/data-1 ]; then
   echo "Copying input data from Acacia"
   mc cp --recursive rhorton/data-1 "$MYSCRATCH/MCProduction"
fi

declare -a idArray
idInd=0
for en in "${enList[@]}"
do
   if [ ! -d "$MYSCRATCH"/MCProduction/EN"${en}" ]; then
      mkdir "$MYSCRATCH"/MCProduction/EN"${en}"
   fi
   cd "$MYSCRATCH"/MCProduction/EN"${en}"

   cp ../input ./data.in
   sed -i "s/yyy/${en}/" data.in
   #Submit job and get the jobId
   jobId=$(sbatch --partition=work -J MC${en}eV --error=MC${en}eV.err $MYSOFTWARE/workflow/runJob.script | cut -f 4 -d ' ')
   #jobId=$(sbatch --partition=work -J TEST  $MYSOFTWARE/workflow/test.sh | cut -f 4 -d ' ')  #Used for debugging, keep for convenience

   #Save ID of submitted job to array
   idArray[${idInd}]="$jobId"
   idInd=$((idInd+1))

   cd ../
done

#Setup list of depencencies for final cleanup script, slurm format is jobId:jobId:jobId:etc for all required jobs
dependencies=""
for id in "${idArray[@]}"
do
   dependencies+=":${id}"
done


#Wait for all jobs to finish successfully(afterok option) and run cleanup script
sbatch --dependency=afterok"${dependencies}"  --partition=work -J MCCleanUp $MYSOFTWARE/workflow/runDataCollect.script
#sbatch --dependency=afterok"${dependencies}"  --partition=work -J MCCleanUp $MYSOFTWARE/workflow/test.sh


