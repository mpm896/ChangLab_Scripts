%{

Matthew Martinez
Created 6/22/2020
Script for running FSC every n particles of a fimal average

From all particles, randomly select n particles
Split n particles into 2 groups of n/2 particles
Average each group of n/2 particles using the last refined table 
Run dynamo_FSC on the two resulting averages
Plot the FSC

Select another n particles, add it to the initial n particles (2n)
Split 2n particles into 2 groups of n particles
Average each group of n particles using the last refined table 
Run dynamo_FSC on the two resulting averages
Plot the FSC
Repeat until all particles are used

%}
clear;
%Enter project name and number of iterations
path2Avg = char(input('Name of project folder: ','s'));
iteLast = char(input('Number of iterations: ','s'));
iteNum = str2num(iteLast);
number = char(input('Number of particles to group: ','s'));
intNum = str2num(number); %Number of particles to group for each interval

%Creates variable with the name of the last refined table
if iteNum < 10
    refinedTable = strcat(path2Avg,'/results/ite_000',iteLast,'/averages/refined_table_ref_001_ite_000',iteLast,'.tbl');
else
    refinedTable = strcat(path2Avg,'/results/ite_00',iteLast,'/averages/refined_table_ref_001_ite_00',iteLast,'.tbl');
end

%Reads the last refined table into the ws variable, rt
rt = dread(refinedTable);
particleNum = size(rt,1); %number of rows of the refined table
rounds = ceil(particleNum/intNum);

data = input('Getting ready to average. Enter data folder: ','s');
amask = input('Enter filename for averaging mask: ','s');
fscmask = input('Enter filename for FSC mask: ','s');
pixSize = input('Enter angstroms per pixel: ','s');

o_table = cell(rounds,1);
rt_table = cell(rounds,1);

old_rt = rt;
for i = 1:rounds
    
    if particleNum < intNum
       rt_table{i,1} = rt;
       break
    end  
  
    X = randperm(particleNum,intNum); %generates 300x1 table with random integers between 1 and particleNum

    new_rt = old_rt(X,:); %Grabs the tags generated in X and puts them in this subset refined table
    rt_table{i,1} = new_rt; %Puts the new values of new_rt into a cell
    old_rt(X,:) = []; %Clears all the rows of tags that were generated randomly
    particleNum = size(old_rt,1); %Determines new # of particles remaining
end

%This block puts all the tables from the cell rt_table into a 3d table
%(except for the last table with remaining particles)
for i = 1:rounds-1
   rt_3d(:,:,i) = rt_table{i,1};
end
rt_last(:,:) = rt_table{rounds,1};

%This block concatenates tables in intervals so you have tables with n,
%2n, 3n, etc particles
rt_concat = cell(rounds,1);
A = [];
for i = 1:rounds
    if i == rounds
       rt_concat{i,1} = rt_last;
       break
    end
    
    A = [A; rt_3d(:,:,i)];
    rt_concat{i,1} = A;
end

%Averages each table, runs FSC, and outputs results into cell
for i = 1:rounds
    output = dynamo_average(data,'table', rt_concat{i,1},'mask',amask,'fsc',1,'fsc_mask',fscmask,'apix',pixSize);
    
    o_table{i,1} = output;
    clear output;
end

figure;
dynamo_fsc_plot(o_table(:),'xticks',10,'style',1);

