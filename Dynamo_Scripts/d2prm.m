%{
03/22/22
Matthew Martinez
Yi-Wei Chang Lab

From a subtomogram alignment project, create a PEET-compatible .prm file with 
all the necessary input parameters to load into a PEET project and quickly 
run a zero search average. This will make it easier to perform various PEET
functions, such as calcFSC, varianceMap, etc.

The tilt range is assumed -60 -> +60. If needed, chance in PEET after
running this program or change in lines 184 + 186

Only necessary input is name of the alignment run and the pixel size of
your tomgorams.
%}

run = char(input('Name of alignment project: ','s'));
pixSize = str2num(input('Enter pixel size in A: ','s'));
runResults = strcat(run,'/results');
dirInfo = dir(runResults); %Takes note of what's in the results directory of your alignment project
tf = ismember({dirInfo.name},{'.','..'}); %looks at file and directory names from input pathway. If '.' or '..', puts into variable tf
dirInfo(tf) = []; %sets rows from tf to empty, basically erasing them
numdir = length(dirInfo);
iteLast = num2str(numdir-1);

newDir = strcat(run,'_prm');
mkdir(newDir);

%Get catalog name
fileCtlg = dir('*.ctlg'); 
ctlg = erase(fileCtlg.name,'.ctlg');

%Get table map name
fileDoc = dir('*.doc');
tblMap = fileDoc.name;

%Open table map file and read tomogram names into cell array (will be in
%cell {1,2}
fid = fopen(tblMap,'r');
tomograms_all = textscan(fid,'%f %s');
fclose(fid)

%Read in the final refined table and average as variables rt and avg
if numdir-1 < 10
    path2rt = strcat(runResults,'/ite_000',iteLast,'/averages/refined_table_ref_001_ite_000',iteLast,'.tbl');
    path2avg = strcat(runResults,'/ite_000',iteLast,'/averages/average_ref_001_ite_000',iteLast,'.em');
else
    path2rt = strcat(runResults,'/ite_00',iteLast,'/averages/refined_table_ref_001_ite_00',iteLast,'.tbl');
    path2avg = strcat(runResults,'/ite_00',iteLast,'/averages/average_ref_001_ite_00',iteLast,'.em');
end
rt = dread(path2rt);
avg = dread(path2avg);

%Get tomogram numbers
tomoNums = unique(rt(:,20));

tomograms = cell(size(tomoNums,1),1);
tiltAngles = cell(size(tomoNums,1),1);
%Find which tomograms were used in the averaging, in case there are extras
%in the table map
for i = 1:size(tomoNums,1)
    if any(tomoNums(:) == i)
        index = find(tomoNums == i);
        tomograms{i,1} = tomograms_all{2}{i};
        tiltAngles{i,1} = '[-60, 60]';
    else
        continue
    end
end

%Create xyz coordinates and Euler angles to create model files and initial
%motive lists
motlCell = cell(size(tomoNums,1),1);
coordCell = cell(size(tomoNums,1),1);
motlHeader = ["CCC","reserve","reserved1","pIndex","wedgeWT","NA","NA1","NA2","NA3","NA4","xOffset","yOffset","zOffset","NA5","NA6","reserved2","EulerZ1","EulerZ3","EulerX2","reserved3"];
motl = [];
modelPoints = [];
for i = 1:size(tomoNums,1)
    coordTable = [];
    angleTable = [];
    count = 1;
    for k = 1:size(rt,1)
        if rt(k,20) == tomoNums(i)
            coord = [rt(k,24)+rt(k,4), rt(k,25)+rt(k,5), rt(k,26)+rt(k,6)];
            coordTable = [coordTable; coord];
            
            angles = [1,0,0,count,0,0,0,0,0,0,0,0,0,0,0,0,-rt(k,9),-rt(k,7),-rt(k,8),0];
            angleTable = [angleTable; angles];
            
            count = count + 1;
        end
    end
    
    coordCell{i} = coordTable;
    angleTable = [motlHeader; angleTable];
    motlCell{i} = angleTable;
    
    nameCoord = strcat('tomo_',string(tomoNums(i)),'_modelPoints.txt');
    writematrix(coordCell{i},strcat(newDir,'/',nameCoord),'Delimiter','space');
    
    nameMOTL = strcat('tomo_',string(tomoNums(i)),'_initMOTL.csv');
    writematrix(motlCell{i},strcat(newDir,'/',nameMOTL));
    
end

path = pwd;
cd(newDir)

%Convert model text files to .mod files
!ls *.txt | awk -F "." '{print "point2model "$1".txt "$1".mod"}' | sh

newDirInfo = dir('.'); %Takes note of what's in the results directory of your alignment project
tf = ismember({newDirInfo.name},{'.','..'}); %looks at file and directory names from input pathway. If '.' or '..', puts into variable tf
newDirInfo(tf) = []; %sets rows from tf to empty, basically erasing them
tfBackup = contains({newDirInfo.name},{'~'}); %looks for backup files to remove
newDirInfo(tfBackup) = [];

mods = contains({newDirInfo.name},{'.mod'}); %Finds model files
motls = contains({newDirInfo.name},{'.csv'}); %Finds initial motive lists

modFiles = []
motlFiles = []
for i = 1:size(mods,2)
    if mods(i) == 1
        modFiles = [modFiles; string(newDirInfo(i).name)];
    end
    if motls(i) == 1
        motlFiles = [motlFiles; string(newDirInfo(i).name)];
    end
end

%Get the reference (average from previous iteration) and mask files and
%correct pixel size
if numdir-1 < 10
    card = dread(strcat('../',run,'/cards/ite_000',iteLast,'/card_iteref_ref_001_ite_000',iteLast,'.card'));
else    
    card = dread(strcat('../',run,'/cards/ite_00',iteLast,'/card_iteref_ref_001_ite_00',iteLast,'.card'));
end

if card.file_mask(1) == '/' %Determines if it is a user input mask with absolute path or if an editor mask
    mask = string(card.file_mask);
else
    mask = strcat(path,'/',card.file_mask);
end 

mask = dread(mask);
dwrite(mask,'mask.mrc');
mask = strcat(path,'/',newDir,'/mask.mrc');

if numdir-2 < 10
    ref = strcat(path,'/',runResults,'/ite_000',iteLast-1,'/averages/average_ref_001_ite_000',iteLast-1,'.em');
else
    ref = strcat(path,'/',runResults,'/ite_00',iteLast-1,'/averages/average_ref_001_ite_00',iteLast-1,'.em');
end
ref = dread(ref);
dwrite(ref,'reference.mrc');
ref = strcat(path,'/',newDir,'/reference.mrc');

cmdMask = sprintf(strcat('alterheader -del %0.2f,%0.2f,%0.2f'," ",mask), pixSize,pixSize,pixSize)
cmdRef = sprintf(strcat('alterheader -del %0.2f,%0.2f,%0.2f'," ",ref), pixSize,pixSize,pixSize)

system(cmdMask)
system(cmdRef)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%CONSTRUCT THE PRM FILE%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prmConstruct = cell(44,1);

prmConstruct{1} = "fnVolume = {'";
prmConstruct{2} = "fnModParticle = {'";
prmConstruct{3} = "initMOTL = {'";
for i = 1:size(tomoNums,1)
    if i ~= size(tomoNums,1)
        prmConstruct{1} = strcat(prmConstruct{1},string(tomograms{i}),''',''');
        prmConstruct{2} = strcat(prmConstruct{2},path,'/',newDir,'/',modFiles(i),''',''');
        prmConstruct{3} = strcat(prmConstruct{3},path,'/',newDir,'/',motlFiles(i),''',''');
    else
        prmConstruct{1} = strcat(prmConstruct{1},string(tomograms{i}),'''}');
        prmConstruct{2} = strcat(prmConstruct{2},path,'/',newDir,'/',modFiles(i),'''}');
        prmConstruct{3} = strcat(prmConstruct{3},path,'/',newDir,'/',motlFiles(i),'''}');
    end
end

prmConstruct{4} = "tiltRange = {[";
for i = 1:size(tomoNums,1)
    if i ~= size(tomoNums,1)
        prmConstruct{4} = strcat(prmConstruct{4},'-60, 60], [');
    else
        prmConstruct{4} = strcat(prmConstruct{4},'-60, 60]}');
    end
end

prmConstruct{5} = "dPhi = {0:0}";
prmConstruct{6} = "dTheta = {0:0}";
prmConstruct{7} = "dPsi = {0:0}";
prmConstruct{8} = "searchRadius = {[0]}";
prmConstruct{9} = "lowCutoff = {[0, 0.05]}";
prmConstruct{10} = "hiCutoff = {[0.15, 0.05]}";
prmConstruct{11} = strcat('refThreshold = {',string(size(rt,1)),'}');
prmConstruct{12} = "duplicateShiftTolerance = [NaN]";
prmConstruct{13} = "duplicateAngularTolerance = [NaN]";
prmConstruct{14} = strcat('reference = ''',ref,'''');
prmConstruct{15} = strcat('fnOutput = ''',ctlg,'''');
prmConstruct{16} = strcat('szVol = [',string(card.feature_data_sidelength),', ',string(card.feature_data_sidelength),', ',string(card.feature_data_sidelength),']');
prmConstruct{17} = "alignedBaseName = ''''";
prmConstruct{18} = "debugLevel = 3";
prmConstruct{19} = strcat('lstThresholds = [',string(size(rt,1)),']');
prmConstruct{20} = "refFlagAllTom = 1";
prmConstruct{21} = "lstFlagAllTom = 1";
prmConstruct{22} = "particlePerCPU = 20";
prmConstruct{23} = "yaxisType = 0";
prmConstruct{24} = "yaxisObjectNum = NaN";
prmConstruct{25} = "yaxisContourNum = NaN";
prmConstruct{26} = "flgWedgeWeight = 1";
prmConstruct{27} = "sampleSphere = 'none'";
prmConstruct{28} = "sampleInterval = NaN";
prmConstruct{29} = strcat('maskType = ''',mask,'''');
prmConstruct{30} = "maskModelPts = []";
prmConstruct{31} = "insideMaskRadius = 0";
prmConstruct{32} = "outsideMaskRadius = NaN";
prmConstruct{33} = "nWeightGroup = 8";
prmConstruct{34} = "flgRemoveDuplicates = 0";
prmConstruct{35} = "flgAlignAverages = 0";
prmConstruct{36} = "flgFairReference = 0";
prmConstruct{37} = "flgAbsValue = 1";
prmConstruct{38} = "flgStrictSearchLimits = 1";
prmConstruct{39} = "flgNoReferenceRefinement = 0";
prmConstruct{40} = "flgRandomize = 0";
prmConstruct{41} = "cylinderHeight = NaN";
prmConstruct{42} = "maskBlurStdDev = 0";
prmConstruct{43} = "flgVolNamesAreTemplates = 0";
prmConstruct{44} = "edgeShift = 1";

prm_final = []
for i = 1:44
    prm_final = [prm_final; prmConstruct{i}; nan];
end

fileName = strcat(ctlg,'_',run,'.prm');
writematrix(prm_final,fileName,'FileType','text','QuoteStrings',0);

