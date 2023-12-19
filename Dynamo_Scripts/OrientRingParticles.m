%{

Matthew Martinez
Created 9/3/2020
Script for orienting particles along a ring (model type General)
Requires a model that is already oriented orthogonally to a surface (so all
Z axes are aligned)

Change difZ (line 67) depending on bin of tomogram.
    for unbinned, difZ -> 300
    for bin4, difZ -> 75

%}

clear;
%Input catalogue name, volume number, source directory name for cropped
%particles, and crop size
catal = char(input('Name of catalogue: ','s'));
volnum = input('Enter volume number corresponding to these models (type an integer): ','s');
modelName = input('Enter (ring) model name: ','s');
modelPath = strcat(catal,'/tomograms/volume_',volnum,'/models/',modelName);
data = input('Enter name of data folder: ','s');
dataFolder = strcat(data,'/crop.tbl');

m = dread(modelPath);
t = m.grepTable();
tCrop = dread(dataFolder);
sz = size(tCrop,1);
szt = size(t,1);
x = 'yes';

%Adjust table to contain only the oriented particles of the model that's
%being oriented here
for i = 1:szt
    for j = 1:sz
        if (tCrop(j,24) > t(i,24)) && (t(i,24) > tCrop(j,24) - 1) && (tCrop(j,25) > t(i,25)) && (t(i,25) > tCrop(j,25) - 1) && (tCrop(j,26) > t(i,26)) && (t(i,26) > tCrop(j,26) - 1)
            t(i,7:9) = tCrop(j,7:9);
        end
    end
end

%Create a point in the center of the ring, oriented in Z with the particle
point = zeros(1,35);

Xmax = max(t(:,24)); Xmin = min(t(:,24)); Xavg = (Xmax+Xmin)/2;
Ymax = max(t(:,25)); Ymin = min(t(:,25)); Yavg = (Ymax+Ymin)/2;
Zmax = max(t(:,26)); Zmin = min(t(:,26)); Zavg = (Zmax+Zmin)/2;
point(1,24)=Xavg; point(1,25)=Yavg; point(1,26)=Zavg;

col7avg = mean(t(:,7)); col8avg = mean(t(:,8));
point(1,7)=col7avg; point(1,8)=col8avg; point(1,9)=0;

figure;
dtplot(t,'m','sketch','sketch_length',100,'sm',30);
dtplot(point,'m','sketch','sketch_length',100,'sm',30);

tOrdered = zeros(size(t,1),35);
tLengths = zeros(size(t,1),1);

%Find which plane the apical ring is oriented (XZ, XY, YZ) by calculating
%the distances between the min and max along each axis, and whichever one
%is the shortest determines the plane (i.e. if the distance along Y is the
%shortest, then the ring is closest to the XZ plane)
tDifference = zeros(1,3);
difX = max(t(:,24))-min(t(:,24)); tDifference(1,1)=difX;
difY = max(t(:,25))-min(t(:,25)); tDifference(1,2)=difY;
difZ = max(t(:,26))-min(t(:,26)); tDifference(1,3)=difZ;

if difZ > 75
    if tDifference(1,1) < tDifference(1,2)
        plane = 'YZ'
    elseif tDifference(1,2) < tDifference(1,1)
        plane = 'XZ' %This is the standard viewing plane
    end
else
    plane = 'XY' %This is a parasite whose tip is facing up or down
end

theta_table = zeros(size(t,1),1);

if plane == 'XZ'
    if col7avg > 5
        for i = 1:size(t,1)
            dx = t(i,24)-point(1,24); %Some geometry to calculate the angle needed to rotate the center point to be oriented towards the particle
            dz = t(i,26)-point(1,26);
            dl = sqrt(dx^2+dz^2);
            theta = acosd(dx/dl);
            theta_table(i,1)=theta;

            if dz<0
                t(i,9)=theta;
            else 
                t(i,9)=0-theta;
            end
        end
    else
        for i = 1:size(t,1)
            dx = t(i,24)-point(1,24); %Some geometry to calculate the angle needed to rotate the center point to be oriented towards the particle
            dz = t(i,26)-point(1,26);
            dl = sqrt(dx^2+dz^2);
            theta = acosd(dx/dl);
            theta_table(i,1)=theta;

            if dz<0
                t(i,9)=theta;
            else 
                t(i,9)=180-theta;
            end
        end
    end
    
    figure;
    dtplot(t,'m','sketch','sketch_length',100,'sm',30);
    dtplot(point,'m','sketch','sketch_length',100,'sm',30);
end

if plane == 'XY'
    for i = 1:size(t,1)
        dx = t(i,24)-point(1,24); %Some geometry to calculate the angle needed to rotate the center point to be oriented towards the particle
        dy = t(i,25)-point(1,25);
        dl = sqrt(dx^2+dy^2);
        theta = acosd(dx/dl);
        theta_table(i,1)=theta;
        
        if dy<0
            t(i,9)=theta;
        else 
            t(i,9)=180-theta;
        end
    end
    figure;
    dtplot(t,'m','sketch','sketch_length',100,'sm',30);
    dtplot(point,'m','sketch','sketch_length',100,'sm',30);
end

if plane == 'YZ'
    if col7avg < 0
        for i = 1:size(t,1)
            dy = t(i,25)-point(1,25); %Some geometry to calculate the angle needed to rotate the center point to be oriented towards the particle
            dz = t(i,26)-point(1,26);
            dl = sqrt(dy^2+dz^2);
            theta = acosd(dy/dl);
            theta_table(i,1)=theta;

            if dz<0
                t(i,9)=theta;
            else 
                t(i,9)=0-theta;
            end
        end
    else
        for i = 1:size(t,1)
            dy = t(i,25)-point(1,25); %Some geometry to calculate the angle needed to rotate the center point to be oriented towards the particle
            dz = t(i,26)-point(1,26);
            dl = sqrt(dy^2+dz^2);
            theta = acosd(dy/dl);
            theta_table(i,1)=theta;

            if dz>0
                t(i,9)=180+theta;
            else 
                t(i,9)=180-theta;
            end
        end
    end
    figure;
    dtplot(t,'m','sketch','sketch_length',100,'sm',30);
    dtplot(point,'m','sketch','sketch_length',100,'sm',30);
end

%{
    
    %Add the model to your volume in your project catalogue
    dcm('-c',catal,'-i',volnum,'-am',m,'-modelname',m_omd)
    
    %create crop table for cropping points in model
    t = m.grepTable();
    allTables{i} = t;
end

tomogramFile = m.cvolume.file();

%Merge all the tables of models in this tomogram and crop particles
o = dtcrop(tomogramFile,t,pName,cropSize);

%}