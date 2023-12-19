%{
Created by Leon, inspired by Matt Martinez (Est. Graduation 2023)
Script to generate plots of alignment distributions per tomogram using the
last refined table of a dynamo run

First: Obtain the path to project/analysis

Second:Choose the number of tomograms to plot
If there is an error, double check you put a number less than or equal to
the number of tomograms available

Third:Choose which tomograms to plot
If there is an error, double check that you selected tomograms in the range
of total available. If there are 14 available, then each number should be
less than or equal to 14.

Fourth: Plot the selected tomograms
You may have to go here to edit the figures. Title, axis labels, makers,
etc.

%}

%Obtain the path the project/analysis run
j = 0;
while j == 0;
    project = char(input('Please enter the path to the project/analysis run: ', 's'));
    project_results = append(project, '/results') %Appends '/results' to the path
    if ~exist(project_results, 'dir');
        disp('Path to project results DOES NOT exist! Please input the path to the directory that contains results (Up 1 directory from results)');
    elseif exist(project_results, 'dir');
        disp('Path to project results exists! YAY!');
        j = 1;
    end;
end;

%Find the number of iterations in this project
all_files = dir(project_results); %Saves all the items of the directory
all_dir = all_files([all_files(:).isdir]); %Saves the items that are directories
num_dir = numel(all_dir) - 2; %Subtracts 2 because there are always '.' and '..' directories


%Obtain dimensions for the refined_cell based on the # of tomgorams
refined_table_str = append(project_results, '/ite_',sprintf('%04d',num_dir),'/starting_values/starting_table_ref_001_ite_',sprintf('%04d',num_dir),'.tbl'); %Saves the path to the refined table
refined_table = load (refined_table_str); %Loads the refined table
last_particle = size(refined_table,1); %Gets the number of the last particle in the refined_table
num_tomo = refined_table(last_particle,20); %Gets the tomogram number from the last particle. This should be the total number of tomograms

refined_cell = cell(1,num_tomo); %Generates a matrix with 1 (row) x number of tomograms

% Fills the refined_cell with values from the refined table to separate the
% orientations on a per tomogram basis
for i = 1:num_tomo; %Loop from 1 to number of tomomgrams
    for k = 1:size(refined_table,1); %For each tomogram loop through each particle of the refined table
        if refined_table(k,20) == i; %If the particle belongs to the tomogram currently being looped
            refined_cell{i} = [refined_cell{i};refined_table(k,:)]; %Add the entire row of that particle to the cell
        end;
    end;
end;


%Determine which tomograms to plot

j = 0; % Set variable to enter the while loop
while j == 0;
    num_tomo_to_plot = input(['Please enter the number of tomograms to plot (Note: You have a max of ' num2str(num_tomo) ' available): ']); %Enter the number of tomograms to plot
    tomo_table = zeros(1,num_tomo_to_plot); % Generate a matrix with the dimensions 1 row x number of tomgorams to plot
    j = 1; % Set to leave the while loop
    for i = 1:num_tomo_to_plot; % For the number of tomograms you would like to plot:
        tomo_number = input('Please enter the tomogram number you would like to plot: '); % Enter the tomogram one at a time in any order.
        tomo_table(i) = tomo_number; % Save the tomogram in the table
    end;
    disp(tomo_table); % Display which tomograms the use would like to use
    valid = input('Are these the tomograms you would like to plot (0 = no, 1 = yes): '); %Confirmation of selection
    while valid ~= 0 && valid ~= 1; % Loops for valid responses
        valid = input('Please enter a valid response (0 = no, 1 = yes): ');
    end;
    if valid == 0; % If these are not the tomograms you want to use
        j = 0;  % Then loop through and ask again
    elseif valid == 1; % If these are the tomograms you want to use
        j = 1; % Then leave the loop
    end;
end;


%Generating the figures
for i = 1:size(tomo_table,2); % Using the tomogram table of those you want to plot
    figure; %Generate a new figure for each tomogram
    dtplot(refined_cell{1,tomo_table(i)},'m','sketch','sm',30,'sketch_length',50); % Use this dynamo function to plot. Change 30 or 50 to change the size of the figure
    fig_title = 'Tomogram %d'; %Generates the title for each figure
    fig_title_complete = sprintf(fig_title, tomo_table(i));
    title(fig_title_complete);
    xlabel('Pixels in X');
    ylabel('Pixels in Y');
    zlabel('Pixels in Z');
end;


