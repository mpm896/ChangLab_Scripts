%{

Matthew Martinez
Created 6/22/2020
Script for running FSC every n particles, ordered by cross correlation
score (column 10), of a final average

From all particles, sort by cc (column 10)
Average the first n particles, 2n particles, 3n particles, etc

%}

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
rtSorted = sortrows(rt,10,'descend'); %Sorts the refined table in descending order by cross correlation (row 10)
particleNum = size(rt,1); %number of rows of the refined table
rounds = ceil(particleNum/intNum);

data = input('Getting ready to average. Enter data folder: ','s');
amask = input('Enter filename for averaging mask: ','s');
fscmask = input('Enter filename for FSC mask: ','s');
pixSize = input('Enter angstroms per pixel: ','s');

o_table = cell(rounds,1);
rt_table = cell(rounds,1);

for i = 1:rounds
    
    if particleNum < intNum
       rt_table{i,1} = rtSorted;
       break
    end  
  
    X = intNum * i; %Number of tags to average for this round

    new_rt = rtSorted(1:X,:); %Grabs X tags from the sorted refined table and puts them in this subset refined table
    rt_table{i,1} = new_rt; %Puts the new values of new_rt into a cell
    particleNum = size(rtSorted,1)-X; %Determines new # of particles remaining
end


%Averages each table, runs FSC, and outputs results into cell
for i = 1:rounds
    output = dynamo_average(data,'table', rt_table{i,1},'fc',1,'mask',amask,'fsc',1,'fsc_mask',fscmask,'apix',pixSize);
    
    o_table{i,1} = output;
end

figure;
dynamo_fsc_plot(o_table(:),'xticks',10,'style',1);
