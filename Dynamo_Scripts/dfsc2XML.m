%{
4/3/2022
Matthew Martinez
Yi-Wei Chang Lab

Script to convert a .fsc file from Dynamo into an XML file that can be
uploaded during deposition of a map to the EMDB

Will grab the .fsc file from the iteration with the best calculated
resolution
%}

run = char(input('Name of alignment project: ','s'));
pixSize = str2num(input('Enter pixel size in A: ','s'));
runResults = strcat(run,'/results');

%Read in bandpass_resolution.txt and find which iteration gave the highest
%calculated resolution
fid = fopen(strcat(runResults,'/bandpass_resolution.txt'),'r');
res_vals = textscan(fid,'%s');

for i = 1:size(res_vals{1},1)
    if i == 1
        best = 1;
        continue
    end
    
    if char(res_vals{1}(i)) >= char(res_vals{1}(best))
        best = i;
    end
end
ite = best;

%Read in the card file and fsc file
if ite < 10
    card = dread(strcat(run,'/cards/ite_000',num2str(ite),'/card_ite_ite_000',num2str(ite),'.card'));
    fsc = dread(strcat(runResults,'/ite_000',num2str(ite),'/averages/bandpass_fsc_ite_000',num2str(ite),'.fsc'));
else    
    card = dread(strcat(run,'/cards/ite_00',num2str(ite),'/card_ite_ite_00',num2str(ite),'.card'));
    fsc = dread(strcat(runResults,'/ite_00',num2str(ite),'/averages/bandpass_fsc_ite_00',num2str(ite),'.fsc'));
end
boxSize = card.feature_data_sidelength;

%Make array of X values (1/Angs) to match those of the FSC

fsc_Y = fsc';
fsc_X = [];
for i = 1:size(fsc_Y,1)
    res = 1/((boxSize*pixSize)/i);
    fsc_X = [fsc_X; res];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%CONSTRUCT THE XML FILE%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


header = string('<fsc title="" xaxis="Resolution (A-1)" yaxis="Correlation Coefficient">');
footer = string('</fsc>');
line1 = sprintf("  <coordinate>");
line2 = sprintf("    <x>0.0</x>");
line3 = sprintf("    <y>1.0</y>");
line4 = sprintf("  </coordinate>");

xmlTable = [];
for i = 1:size(fsc_Y,1)
    str1 = sprintf("  <coordinate>");
    str2 = sprintf("    <x>%f</x>", fsc_X(i,1));
    str3 = sprintf("    <y>%f</y>", fsc_Y(i,1));
    str4 = sprintf("  </coordinate>");
    xmlTable = [xmlTable; str1; str2; str3; str4];
end
xmlTable = [header; line1; line2; line3; line4; xmlTable; footer];

%Write to file
filename = strcat(run,'_fsc.xml');
writematrix(xmlTable,filename,'FileType','text','QuoteStrings',0);





