#!/bin/bash

#define FileList as the first input value
FileList=$1
InputDir=$2
N_Parts=$3
N_Part=$4
OutputFile=$5

#set the number of lines in the file
N_Lines=`wc -l $FileList | awk '{print $1}'`

#set the number of lines per part
N_Lines_Per_Part=`echo "scale=0; $N_Lines / $N_Parts" | bc`

#set the starting line
Start_Line=`echo "scale=0; ($N_Part - 1) * $N_Lines_Per_Part + 1" | bc`

#set the ending line
End_Line=`echo "scale=0; $N_Part * $N_Lines_Per_Part" | bc`

#output variables above into one line
echo "FileList = $FileList"
echo "InputDir = $InputDir"
echo "N_Parts = $N_Parts"
echo "N_Part = $N_Part"
echo "OutputFile = $OutputFile"
echo "N_Lines = $N_Lines"
echo "N_Lines_Per_Part = $N_Lines_Per_Part"
echo "Start_Line = $Start_Line"
echo "End_Line = $End_Line"

#read file in FileList from Start_Line to End_Line, check the file in the InputDir using IsRootFileGood(cpp return value:0 good, 1 bad), generate the command to merge all good files into OutputFile_N_Part.root. print the command to the screen and execute it.
Command_Merge="hadd -f $OutputFile\_$N_Part.root"
Command_Check="IsRootFileGood"
for File in `sed -n "$Start_Line,$End_Line p" $FileList`
do
    $Command_Check $InputDir/$File
    if [ $? -eq 0 ]; then
        Command_Merge="$Command_Merge $InputDir/$File"
    fi
done

echo "================================================="
echo $Command_Merge