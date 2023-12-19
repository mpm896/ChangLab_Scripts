%{
Matthew Martinez
2/23/2021
Yi-Wei Chang Lab

THIS VERSION IS COMPATIBLE WITH PEET 1.13.0. IF USING A LATER VERSION,
CHECK FOR DIFFERENT SCRIPT.

This script will create aligned particles to view in IMOD, similar to what
PEET does with createAlignedParticles.

Line 41 and 42 will both read in the refined table. Only one of these
should be uncommented. Uncomment line X1 if you want to keep the same
binning factor. Uncomment line 42 if you want to change the binning factor.
If you want to change the binning factor, change the number after 'f' as
follows:
    'f' = (source binning) / (destination binning)
    For example, if my particles are unbinned but I want to view them as
    bin 4, then 'f' = 1 / 4 = 0.25



%}

run = char(input('Name of alignment project: ','s'));
runResults = strcat(run,'/results');
dirInfo = dir(runResults); %Takes note of what's in the results directory of your alignment project
tf = ismember({dirInfo.name},{'.','..'}); %looks at file and directory names from input pathway. If '.' or '..', puts into variable tf
dirInfo(tf) = []; %sets rows from tf to empty, basically erasing them
numdir = length(dirInfo);
iteLast = num2str(numdir-1);
particleDir = input('Enter the directory containing your particles: ','s');
newDir = strcat(run,'_alignedParticles');
mkdir(newDir)

    
%Read in the refined table as variable rt
if numdir-1 < 10
    path2rt = strcat(runResults,'/ite_000',iteLast,'/averages/refined_table_ref_001_ite_000',iteLast,'.tbl');
else
    path2rt = strcat(runResults,'/ite_00',iteLast,'/averages/refined_table_ref_001_ite_00',iteLast,'.tbl');
end

rt = dread(path2rt);
%rt = dtrescale(rt,'f',0.25);

keepBoxSize = input('Keep particles the same size? (y/n): ','s');

%Bin particles if desires
binParticles = input('Bin the particles? (y/n): ','s');

if binParticles == 'y'
    cd(particleDir)
    binSize = input('Factor to bin? (enter an integer) ')
    command = strcat('ls particle*.em | awk -F "." ''{print "binvol -b ',{' '},string(binSize),{' '},'"$1".em "$1"_bin4.em"}'' | sh ')
    status = system(command)
end

%Move particles to the aligned particles directory, either with the same
%size or a new crop size
if keepBoxSize == 'y'
    if binParticles == 'n'
        %command = strcat('cp particle_*.em ../',newDir)
        command = strcat ('ls particle_*.em | awk -F "." ''{print "cp "$1".em ../',newDir,'/"$1".mrc"}'' | sh ')
        status = system(command);
        cd('..');
        cd(newDir);
    else
        command = strcat ('ls particle*bin',string(binSize),'.em | awk -F "." ''{print "mv ./"$1".em ../',newDir,'/"$1".mrc"}'' | sh ')
        status = system(command)
        cd('..')
        cd(newDir)
    end
else
    sizeToCrop = char(input('Box size? ','s'));
    tableMapFile = input('Enter table map file: ','s');
    dtcrop(tableMapFile,rt,newDir,sizeToCrop);
    cd(newDir);
    command = strcat ('ls particle_*.em | awk -F "." ''{print "cp ./"$1".em ./"$1".mrc"}'' | sh ')
    status = system(command);
end

%This block makes a table containing 3 columns, each
%pertaining to the Euler angles Z(1), X(2), Z(3). Z(1) and Z(3) are swapped for now between the dynamo table and angle file. It will then save to an
%angles.csv file, be input into MOTL2Slicer, and saved as a new angles csv
%file in XYZ format and read back into the workspace cell slicerTable

EulerAngles = zeros(size(rt,1),3);


for i = 1:size(rt,1)
        EulerAngles(i,1) = rt(i,7);
        EulerAngles(i,2) = rt(i,8);
        EulerAngles(i,3) = rt(i,9);
end

writematrix(EulerAngles,'refinedEulerAngles.csv');

!MOTL2Slicer refinedEulerAngles.csv refinedSlicerAngles.csv
slicerTable = csvread('refinedSlicerAngles.csv');

%Create a table with 6 columns: X,Y,Z,Xangle,Yangle,Zangle
csvSummary = zeros(size(rt,1),6);
for i = 1:size(rt,1)
    if binParticles == 'y'
        csvSummary(i,1) = (rt(i,24)+rt(i,4))/binSize;
        csvSummary(i,2) = (rt(i,25)+rt(i,5))/binSize;
        csvSummary(i,3) = (rt(i,26)+rt(i,6))/binSize;
            
        csvSummary(i,4) = slicerTable(i,1);
        csvSummary(i,5) = slicerTable(i,2);
        csvSummary(i,6) = slicerTable(i,3);
    else
        csvSummary(i,1) = rt(i,24)+rt(i,4);
        csvSummary(i,2) = rt(i,25)+rt(i,5);
        csvSummary(i,3) = rt(i,26)+rt(i,6);
            
        csvSummary(i,4) = slicerTable(i,1);
        csvSummary(i,5) = slicerTable(i,2);
        csvSummary(i,6) = slicerTable(i,3);
    end
end

%Apply rotatevol to each particle
for i = 1:size(rt,1)
    if i < 10 
        num = strcat('00000',string(rt(i,1)));
    elseif i >= 10 && i < 100
        num = strcat('0000',string(rt(i,1)));
    elseif i >= 100 && i < 1000
        num = strcat('000',string(rt(i,1)));
    elseif i >= 1000 && i < 10000
        num = strcat('00',string(rt(i,1)));
    end
    
    if binParticles == 'y'
        command = strcat('rotatevol -a',{' '},string(csvSummary(i,6)),',',string(csvSummary(i,5)),',',string(csvSummary(i,4)),{' '},'particle_',num,'_bin',string(binSize),'.mrc particle_',num,'_bin',string(binSize),'_aligned.mrc')
    else
        command = strcat('rotatevol -a',{' '},string(csvSummary(i,6)),',',string(csvSummary(i,5)),',',string(csvSummary(i,4)),{' '},'particle_',num,'.mrc particle_',num,'_aligned.mrc')
    end
    system(command)
end