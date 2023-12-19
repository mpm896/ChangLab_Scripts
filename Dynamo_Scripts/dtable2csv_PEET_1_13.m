%{
Matthew Martinez
9/17/2020
Y-Wei Chang Lab

This file is to turn a Dynamo particle table into an IMOD-compatible .csv
file. This is originally created to turn a refined table into a .csv that
IMOD can recognize to map particles back into their original tomograms with
their final orientations

1. All you need for input is the name of your alignment project of interest

2. Creates a directory in the folder where it's being a run called
"nameOfProject_AlignedParticleTables"

3. 3 types of files will be generated per tomogram: _EulerAngles.csv
(contains Euler angles in ZXZ format), _SlicerAngles.csv (contains slicer
angles in ZYX format after MOTL2Slicer), and _Summary.csv, which is used in
the clonemodel function

4. In line 43, change binning factor to proper number. For example, if
average is unbined and you want to put into bin 4 tomograms, 'f' = (source
bin) / (destination bin) = 1/4 = 0.25

%}

run = char(input('Name of alignment project: ','s'));
runResults = strcat(run,'/results');
dirInfo = dir(runResults); %Takes note of what's in the results directory of your alignment project
tf = ismember({dirInfo.name},{'.','..'}); %looks at file and directory names from input pathway. If '.' or '..', puts into variable tf
dirInfo(tf) = []; %sets rows from tf to empty, basically erasing them
numdir = length(dirInfo);
iteLast = num2str(numdir-1);

%Read in the refined table as variable rt
if numdir-1 < 10
    path2rt = strcat(runResults,'/ite_000',iteLast,'/averages/refined_table_ref_001_ite_000',iteLast,'.tbl');
else
    path2rt = strcat(runResults,'/ite_00',iteLast,'/averages/refined_table_ref_001_ite_00',iteLast,'.tbl');
end

old_rt = dread(path2rt);
rt = dtrescale(old_rt,'f',0.25);

%Create an empty table of the right dimensions for the csv file
rows = size(rt,1);
initCSV = zeros(rows,20);

%This block counts how many particles there are per tomogram and puts them
%in a verticle array called numTables
count = 0;
numTables = [];
for i = 1:rows
    if i == rows
        diff = rt(i,20)-rt(i-1,20);
        if diff ~= 0
            count = 1; 
            numTables = [numTables; count];
            break
        else
            count = count+1;
            numTables = [numTables; count];
        end
    else 
        diff = rt(i+1,20)-rt(i,20);
        if diff ~= 0
            count = count+1;
            numTables = [numTables; count];
            count = 0;
        else
            count = count+1;
        end
    end
end

%This block will create individual tables from the refined table for each set of particles per
%tomogram
rowNT = size(numTables,1);
count = 0;
tomoRefinedTables = cell(1,rowNT);
for i = 1:rowNT
    tomoRefinedTables{1,i} = rt(count+1:count+numTables(i,1),:);
    count = count+numTables(i,1);
end
    
%This block creates a table of zeros for each tomogram with the correct
%particle number of rows and 20 columns
rowNT = size(numTables,1);
tomoTables = cell(1,rowNT);
for i = 1:rowNT
    tomoTables{1,i} = zeros(numTables(i,1),20);
end

%This block makes a table for each tomogram containing 3 columns, each
%pertaining to the Euler angles Z(1), X(2), Z(3). Z(1) and Z(3) are swapped for now between the dynamo table and angle file. It will then save to an
%angles.csv file, be input into MOTL2Slicer, and saved as a new angles csv
%file in XYZ format and read back into the workspace cell slicerTable
anglesTable = cell(1,rowNT);
for i = 1:rowNT
    temp = tomoRefinedTables{1,i};
    anglesTable{1,i} = zeros(size(temp,1),3);
    tempAngles = anglesTable{1,i};
    
    for j = 1:size(temp,1)
        tempAngles(j,1) = -temp(j,9);
        tempAngles(j,2) = -temp(j,8);
        tempAngles(j,3) = -temp(j,7);
    end
    
    anglesTable{1,i} = tempAngles;
end

newDir = strcat(run,'_AlignedParticleTables');
mkdir(newDir)
cd(newDir)

for i = 1:rowNT
    tomoNum = num2str(tomoRefinedTables{1,i}(1,20));
    nameEulers = strcat('tomo_',tomoNum,'_EulerAngles.csv');
    writematrix(anglesTable{1,i},nameEulers);
end

!ls *EulerAngles.csv | awk -F "Euler" '{print "MOTL2Slicer "$1"EulerAngles.csv "$1"SlicerAngles.csv"}' | sh

slicerTable = cell(1,rowNT);
for i = 1:rowNT
    tomoNum = num2str(tomoRefinedTables{1,i}(1,20));
    slicerName = strcat('tomo_',tomoNum,'_SlicerAngles.csv');
    slicerTable{1,i} = csvread(slicerName);
end

%This block ceates a table of zeros for each tomogram with the correct
%particle number of rows and 7 columns, like an IMOD Summary.csv
headerTable = ["#contour","X","Y","Z","xAngle","yAngle","zAngle"];
csvSummaryTable = cell(1,rowNT);
for i = 1:rowNT
    tempRT = tomoRefinedTables{1,i};
    tempTable = zeros(size(tempRT,1),7);
    tempSlicer = slicerTable{1,i};
    
    for j = 1:size(tempRT,1) %loop to set final positions of XYZ by subtracting offset (columns 4-6 in dynamo table) from original positions (columns 24-26 in dynamo table) and placing in columns 2-4 of IMOD csv
        tempTable(j,1) = j;
        
        tempTable(j,2) = tempRT(j,24)+tempRT(j,4);
        tempTable(j,3) = tempRT(j,25)+tempRT(j,5);
        tempTable(j,4) = tempRT(j,26)+tempRT(j,6);
        
        tempTable(j,5) = tempSlicer(j,1);
        tempTable(j,6) = tempSlicer(j,2);
        tempTable(j,7) = tempSlicer(j,3);
    end
    
    tempTable = [headerTable; tempTable];
    csvSummaryTable{1,i} = tempTable;
end

for i = 1:rowNT
    tomoNum = num2str(tomoRefinedTables{1,i}(1,20));
    nameSummary = strcat('tomo_',tomoNum,'_Summary.csv');
    writematrix(csvSummaryTable{1,i},nameSummary);
end
