#!/bin/bash
#Script: collectData.sh
#Purpose: collects data from the output files 'output.txt and 
#         'diagnostics.txt' produced by the monte carlo code, then
#         writes the data to a files 'results.txt' for plotting 
#         as a function of incident energy
#Author: Reese Horton


#Function: getVal
#Gets the value from the given file with a line containing the given string
#Usage example:  getVal "Mean energy per ion pair" output.txt
getVal() {
   lineNum=`(grep -n "$1" "$2" | cut -c1-4)`
   lineNum=`( echo "$lineNum" | cut -d ":" -f1 )`
   #Piping sequence: reverse order of input string characters, split on spaces and extract first value, reverse result
   #Gets final value in the input line without needing to know the line length or number of spaces
   val=`(sed -n "$lineNum"p "$2" | sed -e 's/[[:space:]]*$//' | rev | cut -d " " -f1 | rev)`
   echo "$val"
}

#Function: collectPoint
#Collects a point of data from the specified files and saves the size of the
#data point to the input array
#Uses bash 4.3 feature (declare -n) to pass arrays by reference
#Usage: collectPoint "Mean energy per ion pair" output.txt widths line
collectPoint() {
   declare -n inWidths=$3
   declare -n inLine=$4

   val=`(getVal "$1" "$2")`
   width=${#val} 
   diff=`( echo "25 - $width +3" | bc )` 
   inLine+="$val" 
   #Pad the line with spaces between data points
   for (( jj=0; jj<${diff}; jj++));
   do
      inLine+=" "
   done

   inWidths+=("$width")
}

cd "$MYSCRATCH"/MCProduction

#Remove old version of results file
rm results.txt

declare -a widths
widths[0]="7"
declare -a headers
headers[0]="EIn(eV)"

readarray -t enArray < energies.txt
for en in ${enArray[@]}
do
   echo "$en"
   cd EN"$en"
   line="$en"
   enLen=${#en} 
   numSpaces=`( echo "21 - $enLen" | bc)`
   for (( ii=0; ii<${numSpaces}; ii++ ))
   do
      line+=" "
   done

   headInd=1
 
   #Ionisation and electronic state excitation data
   collectPoint "reciprocal of mean energy per ion pair:" output.txt widths line
   headers[${headInd}]="1/w  (eV^-1)"
   headInd=$((headInd+1))
   collectPoint "average B1Su excitations per ion pair:" output.txt widths line
   headers[${headInd}]="B1SuExcIonPair"
   headInd=$((headInd+1))
   collectPoint "average C1Pu excitations per ion pair:" output.txt widths line 
   headers[${headInd}]="C1PuExcIonPair "
   headInd=$((headInd+1))
   collectPoint "average singlet excitations per ion pair:" output.txt widths line
   headers[${headInd}]="SingletIonPair"
   headInd=$((headInd+1))
   collectPoint "average triplet excitations per ion pair:" output.txt widths line
   headers[${headInd}]="TripletIonPair"
   headInd=$((headInd+1))

   #Dissociation and vibrational state excitation data
   collectPoint "average X1Sg(v=1) excitations:" output.txt widths line
   headers[${headInd}]="X1Sg(v=1) Exc"
   headInd=$((headInd+1))
   collectPoint "average X1Sg(v=2) excitations:" output.txt widths line 
   headers[${headInd}]="X1Sg(v=2) Exc"
   headInd=$((headInd+1))
   collectPoint "average ratio of X1Sg(v=2) to (v=1) excitations" output.txt widths line
   headers[${headInd}]="X1Sg(v=2,v=1 ratio)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations per ion pair (total):" output.txt widths line
   headers[${headInd}]="DissIonPair(Tot)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations per ion pair (primary):" output.txt widths line
   headers[${headInd}]="DissIonPair(Primary)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations per ion pair (secondary):" output.txt widths line
   headers[${headInd}]="DissIonPair(Secondary)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations per ion pair (singlet):" output.txt widths line
   headers[${headInd}]="DissIonPair(Singlet)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations per ion pair (triplet):" output.txt widths line
   headers[${headInd}]="DissIonPair(Triplet)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations per ion pair (b3Su):" output.txt widths line
   headers[${headInd}]="DissIonPair(b3Su)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations per ion pair (C1Pu):" output.txt widths line
   headers[${headInd}]="DissIonPair(C1Pu)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations per ion pair (B1Su):" output.txt widths line
   headers[${headInd}]="DissIonPair(B1Su)"
   headInd=$((headInd+1))
   collectPoint "ave kinetic energy released through dissociation (eV):" output.txt widths line
   headers[${headInd}]="DissEKRelease"
   headInd=$((headInd+1))
   collectPoint "ave fraction of primary energy deposited as heat through dissociation:" output.txt widths line
   headers[${headInd}]="DissHeatPrimary"
   headInd=$((headInd+1))
   collectPoint "ave fraction of primary energy lost through vibrational excitation:" output.txt widths line 
   headers[${headInd}]="PrimaryVibEn"
   headInd=$((headInd+1))
	 collectPoint "ave dissociations (tot)" output.txt widths line 
	 headers[${headInd}]="Dissociations(Tot)"
   headInd=$((headInd+1))
	 collectPoint "ave dissociations (primary)" output.txt widths line 
	 headers[${headInd}]="Dissociations(primary)"
   headInd=$((headInd+1))
	 collectPoint "ave dissociations (secondary)" output.txt widths line 
	 headers[${headInd}]="Dissociations(secondary)"
   headInd=$((headInd+1))
   collectPoint "ave total number of bound excitations in vcs mode:" output.txt widths line 
	 headers[${headInd}]="BoundExcsInVcs (tot)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations (singlet):" output.txt widths line 
	 headers[${headInd}]="Dissociations (singlet)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations (triplet):" output.txt widths line 
	 headers[${headInd}]="Dissociations (triplet)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations (b3Su):" output.txt widths line 
	 headers[${headInd}]="Dissociations (b3Su)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations (B1Su):" output.txt widths line 
	 headers[${headInd}]="Dissociations (B1Su)"
   headInd=$((headInd+1))
   collectPoint "ave dissociations (C1Pu):" output.txt widths line 
	 headers[${headInd}]="Dissociations (C1Pu)"
   headInd=$((headInd+1))

   #Write the created line to file
   echo "$line" >> ../results.txt
   cd ../
done

#Produce the header for the file based on the widths of the data read
header="${headers[0]}              "
length=${#headers[@]}
for (( ii=1; ii<${length}; ii++));
do
   header+="${headers[${ii}]}"
   headLen=${#headers[${ii}]}
   #diff=`( echo "${widths[${ii}]} - $headLen +3" | bc )` 
   diff=`( echo "25 - $headLen +3" | bc )` 
   for (( jj=0; jj<${diff}; jj++));
   do
      header+=" "
   done
done

#Convoluted method to prepend header to results.txt
echo "$header" | cat - results.txt > temp && mv temp results.txt

#Once file is created, copy it to Acacia for storage
#IMPORTANT: add the date of creation to the object name to avoid overwriting older results
today=$(date '+%Y-%m-%d') 
mc cp results.txt rhorton/simresults/mcresults"${today}".txt 


