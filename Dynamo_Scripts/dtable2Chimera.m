%{
5/4/2022
Matthew Martinez
Yi-Wei Chang Lab

Script to automate the visualization of a subtomogram average in UCSF 
Chimera/ChimeraX at the refined table positions 

1. It is assumed that the average is unbinned and that the binning for
depiction is 4x binned (***Adjust in lines 57-62).

2. The isosurface level is set, by default, to 6. If the average for
visualization has a bad thresholding, change this value (***Line 77)

4. Pixels below a certain threshold are masked out. If it looks like too
many pixels were masked out, change this value (***Line 71)

5. *** MUST PROVIDE A CUSTOM MASK IN ORDER TO CREATE A VOLUME ENVELOPE THAT
CONTAINS JUST ONE SUBUNIT ***

Particles are saved as child models within a parent model in ChimeraX
%}
PIX_SIZE = 10.617 %In angstroms. This is for depiction purposes in ChimeraX
homeDir = string(pwd);

%Get the table of interest
run = char(input('Name of alignment project or specific table file (ending in .tbl): ','s'));

%Determine if alignment project or specific table of interest and grab
%table and average of interest
x = strfind(run,'.tbl');
if ~isempty(x)
    t = dread(run);
    run = char(input('Name of alignment project with desired average: ','s'));
end

runResults = strcat(run,'/results');
dirInfo = dir(runResults);
tf = ismember({dirInfo.name},{'.','..'});
dirInfo(tf) = [];
numdir = length(dirInfo);
iteLast = num2str(numdir-1);

avgPath = ddb([run ':a']);
avg = dread(avgPath);

if isempty(x)
    tPath = ddb([run ':rt']);
    t = dread(tPath);
end

%WORK ON THIS. SO FAR IT SEEMS LIKE THIS PUTS THE AVERAGE INTO THE PROPER
%SPOT WITH PROPER SCALING, BUT BECAUSE THE AVERAGE IS NOT AN .MRC FILE IT
%DOES NOT HAVE THE PIXEL INFO AND THUS APPEARS VERY SMALL IN THE TOMOGRAM
%{
%Multiply coordinates by pixel size to match scaling with tomogram and any
%IMOD models that will also be depicted
t(:,4:6) = t(:,4:6).*PIX_SIZE;
t(:,24:26) = t(:,24:26).*PIX_SIZE;
%}
%Ask for custom mask and read it into workspace
maskFile = char(input('Name of mask file: ','s'));
ms = dread(maskFile);

%Make new directory where chimera .stl files will go
newDir = strcat(run,'_Chimera');
mkdir(newDir);

%Bin the volumes and the table from *unbinned* to *bin4*
ms = dynamo_bin(ms, 2);
avg = -avg; 
avg = dynamo_bin(avg, 2);
t = dtrescale(t, 'f', 0.25);

%%%
%Create envelope of the subtomogram average
%%%

%Mask and normalize the volume and eliminate voxels with denisities below
%the desired threshold
avgNorm = dynamo_normalize_roi(avg);
mask = zeros(size(avgNorm));
thresh = -0.1;
includedIndices = find(avgNorm>thresh);
mask(includedIndices) = avgNorm(includedIndices);
avgNorm = avgNorm.*mask;

%Create the envelope from the normalized/masked average
isolevel = 6;
dt = dynamo_isosurface(avgNorm.*ms,'isolevel',isolevel,'real_isolevel',true,'-show',false);

%Smooth the surface
tSmooth = mbgeom.triangulation.smoothingSubdivision(dt.tr,'mr',1);

%Remove the peaks created by the smoothing algorithm
tClean = mbgeom.triangulation.removePeaks(tSmooth,'std',2);

%Remove dust (unconnected pieces)
tNoDust = mbgeom.triangulation.removeDust(tClean);

%Impose a consistent orientation for the normals of the triangles
tFinal = mbgeom.triangulation.consistentNormals(tNoDust,'alignTo','out');

%Visualize the 4 envelopes that were just created
figure; f = gcf; f.Name = 'Subtomogram Average Envelope';
triangulationsToPlot = {dt.tr,tSmooth,tClean,tNoDust};
titles = {'Direct triangulation', 'Smoothed', 'Peaks removed', 'Dust removed'};
for i = 1:4
    ax(i) = subplot(2,2,i);
    h(i) = trisurf(triangulationsToPlot{i});
    axis('equal');
    h(i).LineStyle = 'none';
    shading(ax(i), 'interp');
    lightangle(-45,30);
    axis(ax(i),'off');
    ht(i) = title(ax(i),titles{i});
    ht(i).Color = 'w';
    mbgraph.cursors.setMouse3d(ax(i));
end

f.Color = 'k';

%Split up the table by tomogram and map the average into the tomograms
nTomo = unique(t(:,20));
tCell = cell(size(nTomo,1),1);
tMeshCell = cell(size(nTomo,1),1); %Cell for the tables with meshes
rc = (size(avgNorm,1)/2) + 0.5; %Volume rotation center. Very important for proper placement

stat = system("chimerax --version")
for i = 1:size(nTomo,1)
    if i == 1
        disp(i)
    end
    
    %Split up the table by tomogram
    for k = 1:size(t,1)
        if t(k,20) == nTomo(i)
            tCell{i} = [tCell{i}; t(k,:)];
        end
    end
    
    %Map the average into the table positions
    tCell{i}(:,24:26) = round(tCell{i}(:,24:26));
    tMeshCell{i} = dpktbl.triangulation.place(tCell{i},tFinal,'rc',rc);
    
    %Name and write the file
    fileName = strcat('AverageMesh_tomo_',string(nTomo(i)),'.stl');
    path2file = strcat(newDir,'/',fileName);
    dwrite(tMeshCell{i},path2file);
    
    %Check if the ChimeraX executable exists. If not, skip this chunk of
    %code
    if ~stat
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Create a file with the subtomogram average placed at the aligned positions
        % in ***ChimeraX***. This is not the surface file! Also this is
        % specific for ChimeraX. For regular UCSF Chimera this process is
        % easier
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %First run the command used for OG chimera to get the file with the
        %positions matrix, which will need to be reformatted for ChimeraX
        dtchimera(tCell{i},'-template',avgNorm.*ms,'-show',0)

        %Need to convert the OG Chimera position matrix file to one that is
        %compatible with ChimeraX
        posFile = dread("temp_dynamo_chimera_dir_transform.cmd");
        newPositions = [];

        %Reformat the position matrix to ChimeraX format
        for k = 1:size(posFile,2)
            if strfind(posFile{1,k},'Model')
                A = strsplit(posFile{1,k+1}); A(1) = [];
                B = strsplit(posFile{1,k+2}); B(1) = [];
                C = strsplit(posFile{1,k+3}); C(1) = [];

                line2write = strcat(string(posFile{1,k}),",",A{1},",",A{2},",",A(3),",",A(4),",",B{1},",",B{2},",",B{3},",",B{4},",",C{1},",",C(2),",",C{3},",",C{4});
                newPositions = [newPositions; line2write];
            end

            if k == size(posFile,2)
                filename2write = strcat('tomo_',string(nTomo(i)),'.positions');
                path2posFile = strcat(homeDir,'/',newDir,'/',filename2write);
                writematrix(newPositions,path2posFile,'FileType','text','QuoteStrings',0);
            end
        end

        %Construct the command file
        commandFile = strcat(newDir,'/tomo_',string(i),'_ChimeraX.cxc');
        volNum = size(tCell{i},1);
        line1 = "open";
        for k = 1:size(tCell{i},1)
            line1 = strcat(line1,{' '},homeDir,'/temp_chimera_table.em');    
        end
        line2 = "volume #* level 6";
        line3 = strcat('open',{' '},path2posFile,{' '},'model #1.*.* childModels true'); 
        line4 = "view orient";
        line5 = "color #* #30A6AD";
        line6 = "set bgColor white";

        lines2write = [line1; line2; line3; line4; line5; line6];
        writematrix(lines2write,commandFile,'FileType','text','QuoteStrings',0);
    end
    
    
end

%Launch ChimeraX on the first model
if ~stat
    cmd = strcat('chimerax',{' '},string(newDir),'/tomo_1_ChimeraX.cxc &');
    system(cmd);
end
