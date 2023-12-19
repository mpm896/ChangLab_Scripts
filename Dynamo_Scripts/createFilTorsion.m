%{

Matthew Martinez
Created 6/2/2020
Script for particles along a filament (filament models)

Read IMOD model text files into workspace variables
Create a filamentWithTorsion model for each IMOD model
Add and update crop points to Dynamo filament model from the IMOD model
Add to catalogue
Update crop points to fit the filament dimensions
Extract table for each filament model for one tomogram
Merge tables, crop all the points

**PREREQUISITES**
1. Convert all model files (.mod) to .txt files with the PEET command model2point. For example: 
    ls *.mod | awk -F "." '{print "model2point "$1".mod "$1".txt"}' | sh

2. Put all your .txt model files into a .csv file. For example, I used:
 ls */imod/ModelFile/txtModel/*.txt > imodModels.csv

As of now, this .csv containing all your .txt models should be the only
.csv in your parent directory

3. In your tableMap.doc file, tomograms must be listed in the SAME ORDER as
the model files when you use the 'ls' command in the terminal

4. If you want to crop from tomograms with different bin factor than model
update line 103 with proper scaling factor. Otherwise comment out.

5. In line 111, update with the proper tableMap.doc file for your project

***NOTE: IMPORTANT*****
Prior to running this, make sure to edit values for dz (pixel spacing) at line 90,
dphi (angular twist with respect to the previous particle) at line 89

%}

model_path = input('Enter the pathway to the IMOD model txt files (without final /): ','s'); % 's' doesn't check your input type
catal_path = char(input('Enter pathway to directory that contains your catalogue (without final /): ','s')); 
file = dir([model_path,'/*.csv']) %Finds the .csv file containing all your .txt models to be imported
filename = strcat(file.folder,'/',file.name);
dirinfo = dir(model_path); % info from the directory in model_path put into variable dirinfo

cd(catal_path) %change to the directory of input

fid = fopen(filename,'r'); %Opens the .txt or .csv file containing all your .txt model files
txtModels = textscan(fid,'%s'); %reads all lines into a cell
fclose(fid);

numfiles = length(txtModels{1,1}); %number of model .txt files

%Input catalogue name, volume number, source directory name for cropped
%particles, and crop size
catal = char(input('Name of catalogue: ','s'));
volnum = input('Enter volume number corresponding to these models (type an integer): ');
pName = char(input('Name of source directory for cropped particles: ','s'));
cropSize = input('Particle crop size: ');

allTables = cell(1,numfiles); %Cell to input all tables for merging

%Put all .txt model files with full pathway into a cell
allFiles = cell(1,numfiles);
for i = 1:numfiles
    allFiles{1,i} = strcat(model_path,'/',txtModels{1,1}{i});
end

%Placeholder cells for initial table (ri), table after first 2 columns are
%fixed (rm), and final table (rf) for inspection
ri = cell(1,numfiles);
rm = cell(1,numfiles);
rf = cell(1,numfiles);


for i = 1:length(allFiles)
    fprintf(allFiles{i});
    r = readtable(allFiles{i}); %read the IMOD txt file into r
    
    
    %print(' ');

    rf{i} = table2array(r); %places final table of model points into cell rf
    %rf{i} = r; %places final table of model points into cell rf
    
    %Generate filamentWithTorsion model, add points to it from IMOD, set
    %spacing to 15 pixels, create crop points
    m = dmodels.filamentWithTorsion();
    m.addPoint(rf{i});
    m.subunits_dphi = 0;
    m.subunits_dz = 14;
    m.updateCrop();
    
    %This block of lines creates the model name with the proper extension
    m_name = erase(allFiles(i),'txt');
    m_omd = char(strcat(m_name,'omd'));
    dest = char(strcat(model_path,'/',m_omd));
    
    %Add the model to your volume in your project catalogue
    dcm('-c',catal,'-i',volnum,'-am',m,'-modelname',m_omd)
    
    %create crop table for cropping points in model
    t = m.grepTable();
    %t_rescaled = dtrescale(t,'f',8); %Rescales your points to a different bin factor
    t(1,20) = volnum;
    allTables{i} = t;
    volnum = volnum + 1;
end

%Merge all the tables of models in this tomogram and crop particles
tGlobal = dtmerge({allTables{:}},'linear_tags',1);
oGlobal = dtcrop('APR1.doc',tGlobal,pName,cropSize);
osb = daverage(pName,'t',tGlobal);
