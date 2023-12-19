%{

Matthew Martinez
Created 9/13/2023
Script for cropping particles along a filament (filament models)

1. Read in saved models from catalogue
2. Update crop points to Dynamo filament model with dphi and dz
3. Save to catalogue (save with new name for now so we don't mess up original
files
3. Extract table for each filament model
4. Merge tables, crop all the points

**PREREQUISITES**
1. Models must already be made and saved into catalogue
2. Run from the directory that contains the catalog (.ctlg file)

%}

% Check if directory contains a catalogue, get catalogue name,
% and set working path
ctlg_files = dir('*.ctlg')
if size(ctlg_files, 1) == 1
    ctlg_path = pwd;
    ctlg = ctlg_files(1).name;
elseif size(ctlg_files, 1) > 1
    ctlg_path = pwd

    % Also get catalogue name
    ctlg = char(input('Name of catalogue: ','s'));
else
    fprintf("Must run from directory containing the Dynamo catalogue")
    return
end

% Get Number of tomograms
tomo_path = strcat(ctlg_path, '/', ctlg, '/tomograms');
tomograms = dir(tomo_path);
tf = ismember({tomograms.name},{'.','..'}); 
tomograms(tf) = []; % Erase '.' and '..'
numtomos = length(tomograms);


% Read in all models for each tomogram, update crop points with desired
% dphi and dz, and crop the particles
dcmodels(ctlg, 'ws', 'model_files')
dphi = input("Enter dphi: ");
dz = input("Enter dz: ");

ext = '_xr556'; % Extension for new files. Random chars

T_table = cell(length(model_files.files), 1);
for i = 1:length(model_files.files)
    % Skip file if it is a new file with ext
    file_check = strfind(model_files.files(i), ext);
    if ~isempty(file_check{1})
        continue
    end


    % Modify the model
    m = dread(model_files.files(i));
    m.subunits_dphi = dphi;
    m.subunits_dz = dz;
    m.updateCrop();
    
    % Get uupdated table
    t = m.grepTable();
    T_table{i} = t;

    volume = m.cvolume.file();
    
    % Get volume index
    k = strfind(model_files.files(i), 'volume_');
    index = model_files.files(i);
    index = index{1};
    index = index(k{1}+7);

    % Adjust file name and write to the original location
    fileName = model_files.files(i);
    fileName = erase(fileName,'.omd');
    newFileName = char(strcat(fileName,ext,'.omd'));

    dwrite(m, newFileName);

end

% Combine tables into one and crop
doc_file = dir('*.doc');
if size(doc_file, 1) == 1
    doc_table = doc_file(1).name
else
    doc_table = char(input('Name of .doc file: ','s'));
end

newDir = char(input('Name of new directory for cropping: ','s'));
cropSize = input('Particle size for cropping: ');

tGlobal = dtmerge({T_table{:}},'linear_tags',1);
oGlobal = dtcrop(doc_table, tGlobal, newDir, cropSize);

    






    