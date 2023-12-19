%{
Matthew Martinez
11/4/2020
Yi-Wei Chang Lab

This script can be used to apply symmetry around the rotational axis ("Z" axis or Azimuth) to your average. This is done similarly to PEET's symmetry operation.
 
In Dynamo, when you apply symmetry during an alignment project, it applies
symmetry to the reference to increase the chances of your particles
aligning in different positions, giving the appearance of increased
symmetry in the average.

In PEET, particles number are multiplied by the amound of symmetry you want
to apply and rotated by 360/symmetry # around the rotational axis. For
example, if you have 10 particles and want to apply 8-fold symmetry, you
will have 80 particles in the end. Every set of 10 particles will be
rotated by 360/8 (45) degrees from the last set. This is what this script
will do. 

%}

project = char(input('Which alignment project would you like to symmetrize? ','s'));
iteLast = char(input('Number of iterations: ','s'));
iteNum = str2num(iteLast);
sym = str2num(input('What symmetry would you like to apply (integer)? ','s'));
angle = 360/sym;
dataFolderName = char(input('Enter name of new data folder: ','s'));
tableMap = char(input('Enter name of table map file: ','s'));
cropSize = str2num(input('Enter particle crop size: ','s'));

if iteNum < 10
    refinedTable = strcat(project,'/results/ite_000',iteLast,'/averages/refined_table_ref_001_ite_000',iteLast,'.tbl');
    average = strcat(project,'/results/ite_000',iteLast,'/averages/average_ref_001_ite_000',iteLast,'.em');
else
    refinedTable = strcat(project,'/results/ite_00',iteLast,'/averages/refined_table_ref_001_ite_00',iteLast,'.tbl');
    average = strcat(project,'/results/ite_000',iteLast,'/averages/average_ref_001_ite_00',iteLast,'.em');
end

rt = dread(refinedTable);
particleNum = size(rt,1);

tableCopies = cell(1,sym);
rtNew = [];

for i = 1:sym
    tableCopies{1,i} = rt;
    
    for k = 1:particleNum
        tableCopies{1,i}(k,1) = k + (particleNum*(i-1));
        tableCopies{1,i}(k,9) = tableCopies{1,i}(k,9) + (angle*(i-1));
    end
    
rtNew = [rtNew; tableCopies{1,i}];
end

dtcrop(tableMap,rtNew,dataFolderName,cropSize);



