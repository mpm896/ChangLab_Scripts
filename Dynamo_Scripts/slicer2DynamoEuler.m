%{
11/1/2020
Matt Martinez
Yi-Wei Chang Lab

Convert slicer angles into Dynamo-compatible angles and insert into
particle crop table

*** PREREQUISITES ***

Before this, run makeMOTL.sh or slicer2MOTL to convert slicer angles into
Euler angles, and put into a PEET-style motive list (so that angles are in
the proper columns in PEET convention)

Put the PEET-style motive lists (in the case of makeMOTL.sh, the
initMOTL.csv files) into a subdirectory of the Dynamo project. Each motive
list must be named differently (hint: add the tomogramID to the front of
each name) so that they can all be added to the same directory

In line 41, input directory name in 'dir('')' of the directory containing the .csv files
In line 52, input directory name after 'cd' of the directory containing the .csv files
In line 77, change the input in 'dread' to the directory containing your paticles + '/crop.tbl'
In line 86, change the input in 'dread' to the directory containing your paticles + '/crop.tbl'
In line 88, change the input in 'daverage' to the directory containing your paticle

If you have your tomograms sorted by number, uncomment line 45 and the
block of lines from line 56-64, and comment out the block of lines from line 68-71. 
If you do not have your tomograms sorted by number (but rather sorted however
they showup with the 'ls' command in the terminal), uncomment lines 68-71
and comment out line 45 and lines 56-64.

If you want your tomograms sorted by number, you need:
natsort.m
natsortfiles.m

%}

% !mkdir ParticleAngles
% !ls /ChangLab1-hd2/matt/Plasmodium/Volumes/ManualRecon/*/imod/txtModel*/*initMOTL.csv | awk '{print "cp "$1" ./ParticleAngles/"}' | sh

dirInfo = dir('ParticleAngles_Nov2020_Set2'); %Takes note of what's in the results directory of your alignment project
tf = ismember({dirInfo.name},{'.','..'}); %looks at file and directory names from input pathway. If '.' or '..', puts into variable tf
dirInfo(tf) = []; %sets rows from tf to empty, basically erasing them
numdir = length(dirInfo);
%csvSorted = natsortfiles({dirInfo.name});

csvTable = cell(1,numdir);

MOTL = [];


cd ParticleAngles_Nov2020_Set2

%read in csv files and concatenate them into one large csv file in order.
%COMMENT OUT IF LEAVING UNSORTED
%for i = 1:numdir
 %   for j = 1:numdir
  %      filename = convertCharsToStrings(dirInfo(j).name)
   %     if csvSorted{1,i} == filename
    %        csvTable{1,i} = readmatrix(dirInfo(j).name);
     %       MOTL = [MOTL; csvTable{1,i}];
      %  end
    %end
%end

%read in csv files and concatenate them into one large csv file in order.
%COMMENT OUT IF YOU WANT TO SORT IN ORDER BY NUMBER
for i = 1:numdir
    csvTable{1,i} = readmatrix(dirInfo(i).name);
    MOTL = [MOTL; csvTable{1,i}];
end



cd ../

cropTable = dread('Particles_Nov2020_Set2/crop.tbl');

%move angles from csv file to particle crop table
for i = 1:size(MOTL,1)
    cropTable(i,7) = -MOTL(i,18);
    cropTable(i,8) = -MOTL(i,19);
    cropTable(i,9) = -MOTL(i,17);
end

dwrite(cropTable,'Particles_Nov2020_Set2/crop.tbl'); %write the crop.tbl file with the new angles

osb = daverage('Particles_Nov2020_Set2','t',cropTable); %create a zero search initial average

%dwrite(osbManual.average,'AverageZeroSearch_bin4.mrc');

