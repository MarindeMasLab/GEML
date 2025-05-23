%{
Integrating Machine Learning with Metabolic Models for Precision Trauma Care: Relevant Reaction Analysis
Authors:

Igor Marin de Mas (Copenhagen University Hospital, Rigshospitalet)
Lincoln Moura (Universidade Federal do Ceará)
Fernando Luiz Marcelo Antunes (Universidade Federal do Ceará)
Josep Maria Guerrero (Aalborg University)
Pär Ingemar Johansson (Copenhagen University Hospital, Rigshospitalet)

Introduction

This script focuses on analyzing relevant metabolic reactions in genome-scale metabolic models. 
The goal is to identify key reactions, associated genes, and metabolites that play a significant 
role in patient-specific metabolic responses. This analysis supports the identification of potential 
metabolic targets for precision trauma care.

Workflow Overview

1. Load Data: Load the workspace containing preprocessed data for relevant reactions.
2. Define Relevant Reactions: Specify the subset of reactions to analyze based on predefined criteria.
3. Extract Reaction Properties: Retrieve reaction IDs, names, EC numbers, subsystems, and gene-protein-reaction (GPR) rules.
4. Gene Analysis: Identify genes associated with relevant reactions and map them to gene names using the NCBI API.
5. Metabolite Analysis: Extract substrates and products for each relevant reaction.
6. Subsystem Analysis: Summarize the distribution of relevant reactions across metabolic subsystems.
7. Summary Tables: Generate summary tables for reactions, genes, metabolites, substrates, and products.
8. Save Results: Save the analysis results to an Excel file for further exploration.

Libraries Used

- MATLAB: Used for implementing the workflow, performing data extraction, and generating summary tables.
- NCBI API: Accessed to map Entrez gene IDs to gene names for better interpretability.

Key Features

- Comprehensive analysis of relevant metabolic reactions, including associated genes and metabolites.
- Integration with the NCBI API for gene name mapping.
- Generation of detailed summary tables for reactions, subsystems, substrates, and products.
- Export of results to Excel for easy sharing and further analysis.
- Modular design for easy customization and reuse.
%}

% Load the data file
load("relevant_reactions_analysis_source_workspace.mat")

% Define the relevant reactions
% K-fold
%relevant_reactions_general = [229, 289, 1189, 1927, 893, 224, 2397, 823, 423, 2744, 1683, 2467, 2544, 343, 1560, 1315, 2592, 2249, 196, 2393, 329, 865, 2306, 2420, 2495, 1557, 1334, 1616, 885, 2305, 2394, 2882, 200, 480, 1185, 1782, 2450, 1556, 1321, 1958, 1163, 2276, 1951, 151, 2466, 311, 258, 2286, 2832, 1930, 1812, 2448, 2404, 457, 156, 1534, 2441, 169, 2474, 32, 2603, 1202, 1298, 590, 122, 1319, 2425, 2235, 84, 2568, 493, 2354, 960, 2481, 2388, 2378, 2266, 1630, 239, 2429, 2353, 1467, 988, 2166];
%relevant_reactions_female = [2451, 2163, 2410, 2453, 2802, 2499, 2440, 1201, 2472, 2007, 1161, 1334, 1782, 2424, 2447, 2249, 1143, 2465, 256, 2397, 2427, 2353, 1314, 2250, 250, 263, 2288, 2432, 1812, 156, 2394, 2393, 1163, 1537, 2486, 415, 1299, 1967, 258, 1764, 89, 151, 2436, 2235, 2461, 2467, 416, 1767, 1639, 1250, 2276, 2443, 375, 2484, 926, 96, 1665, 2459, 2584, 1298, 2167, 1595, 2444, 2454, 2430, 2277, 2166, 2420, 194, 364, 344, 2460, 1202, 2449, 2412, 2439, 457, 114, 1528];
%relevant_reactions_male = [598, 2882, 2390, 1315, 2744, 2249, 2397, 2453, 2235, 1927, 2354, 447, 1419, 2393, 2464, 2228, 423, 2443, 2305, 1185, 2388, 2300, 1314, 2425, 2467, 1827, 462, 1319, 2478, 2290, 2250, 2495, 2474, 1189, 2420, 1334, 2434, 232, 2446, 2030, 344, 2481, 1951, 2396, 107, 1143, 755, 426, 1163, 2288, 2193, 1747, 2466, 457, 1053, 885, 251, 2838, 32, 2272, 1764, 2276, 2084, 1683, 2394, 1594, 1639, 1782, 414, 2439, 2780, 2429, 865, 1338, 2306, 1930, 2197, 2015, 1556, 2452, 460, 1565];

% SHAP
relevant_reactions_shape = [1930, 1927, 1951, 1958, 2249, 2467, 1616, 2394, 2429, 2436, 2444, 2397, 480, 2227, 1163, 1234, 200, 1334, 2055, 2544, 1130, 1185, 2744, 2450, 2466, 32, 457, 2474, 1419, 2305, 2388, 89, 534, 2838, 2404, 2419, 2730, 2485, 2393, 2461, 462, 1665, 2497, 464];
relevant_reactions_shape_subset = [1616,1185, 2744, 2838, 1927, 2055, 1419];

%relevant_reactions = relevant_reactions_general;
%relevant_reactions = relevant_reactions_female;
%relevant_reactions = relevant_reactions_male;
relevant_reactions = relevant_reactions_shape_subset;


% Extract properties of the relevant reactions
relevant_rxn_id = sampleMetaOutC.rxns(relevant_reactions);
relevant_rxn_name = sampleMetaOutC.rxnNames(relevant_reactions);
relevant_rxn_EC = sampleMetaOutC.rxnECNumbers(relevant_reactions);
relevant_rxn_subsystem = sampleMetaOutC.subSystems(relevant_reactions);
relevant_rxn_gpr = sampleMetaOutC.grRules(relevant_reactions);

%%% Genes %%%
% Identify genes associated with the relevant reactions
relevant_reaction_genes = sampleMetaOutC.rules(relevant_reactions);

% Initialize a cell to store the extracted genes
extracted_genes = cell(size(relevant_reactions))';

% Iterate over the relevant reactions and extract the GPR rules
for g = 1:length(relevant_reactions)
    gth_gpr = sampleMetaOutC.rules{relevant_reactions(g)};
    % Use regular expressions to find all numbers
    genes = regexp(gth_gpr, '\d+', 'match');
    % Convert the found numbers to integers
    genes = str2double(genes);
    % Store the result in the cell
    extracted_genes{g} = genes;
end

% Convert Entrez IDs to gene names using the NCBI API
ncbi_gene_names = cell(size(extracted_genes));


for i = 1:length(extracted_genes)
    entrez_ids = extracted_genes{i};
    gene_names = {};
    for j = 1:length(entrez_ids)
        entrez_id = entrez_ids(j);
        url = sprintf('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%d&retmode=json', entrez_id);
        success = false;
        retry_count = 0; % Counter for retries
        max_retries = 5; % Set a maximum number of retries
        while ~success && retry_count < max_retries
            try
                options = weboptions('Timeout', 60); % Set the timeout to 60 seconds
                data = webread(url, options);
                % Access the dynamic field within data.result
                field_name = ['x', num2str(entrez_id)];
                if isfield(data.result, field_name) && isfield(data.result.(field_name), 'name')
                    gene_name = data.result.(field_name).name;
                    gene_names{end+1} = gene_name; % Add the gene name to the cell array
                end
                success = true; % Successfully retrieved the data
            catch ME
                if strcmp(ME.identifier, 'MATLAB:webservices:HTTP429StatusCodeError')
                    pause(1); % Wait 1 second before retrying for 429 error
                elseif strcmp(ME.identifier, 'MATLAB:webservices:Timeout')
                    pause(5); % Wait 5 seconds before retrying for timeout
                elseif contains(ME.message, '500')
                    retry_count = retry_count + 1; % Increment retry count for HTTP 500 error
                    pause(2 * retry_count); % Exponential backoff, wait before retrying
                    if retry_count >= max_retries
                        fprintf('HTTP 500 error for entrez_id: %d. Skipping after %d retries.\n', entrez_id, max_retries);
                        break; % Break out of the loop after max retries
                    end
                else
                    rethrow(ME); % Re-throw the error for unknown issues
                end
            end
        end
    end
    % Convert the cell array of gene names to a string vector
    ncbi_gene_names{i} = strjoin(gene_names, ', ');
end

relevant_reaction_entrez_genes = extracted_genes;
relevant_reaction_ncbi_genes = ncbi_gene_names;

%%% Metabolites %%%
relevant_substrates = cell(size(relevant_reactions))';
relevant_products = cell(size(relevant_reactions))';
relevant_M = full(sampleMetaOutC.S);
relevant_M = relevant_M(:,relevant_reactions);

for m = 1:length(relevant_reactions)
    mth_S = relevant_M(:,m);
    mth_index_substrate = find(mth_S == -1);
    mth_index_product = find(mth_S == 1);
    mth_relevant_substrates = sampleMetaOutC.metNames(mth_index_substrate);
    mth_relevant_products = sampleMetaOutC.metNames(mth_index_product);
    % Convert the cell of metabolite names to a string vector
    relevant_substrates{m} = strjoin(mth_relevant_substrates, '; ');
    relevant_products{m} = strjoin(mth_relevant_products, '; ');
end

%% Summary Subsystem
rxn_subsystems = cell2table(sampleMetaOutC.subSystems); % Extract the 'RxnSubsystem' column from sampleMetaOutC
[unique_subsystems, ~, idx] = unique(rxn_subsystems); % Get unique subsystems and their indices
filtered_idx = idx(relevant_reactions(:)); % Filter idx to get only the relevant indices

% Get the IDs of the relevant reactions
relevant_rxn_id = relevant_rxn_id(:); % Ensure relevant_rxn_id is a column vector
list_reactions_in_subsystem = accumarray(filtered_idx, (1:length(filtered_idx))', [], @(x) {strjoin(relevant_rxn_id(x), ', ')}); % List the reactions in each subsystem for the relevant reactions
% Remove subsystems without reactions
non_empty_idx = ~cellfun('isempty', list_reactions_in_subsystem);
list_reactions_in_subsystem = list_reactions_in_subsystem(non_empty_idx);
unique_subsystems = unique_subsystems(find(non_empty_idx),:);
num_reactions_in_subsystem = accumarray(idx(relevant_reactions(:)), 1);
num_reactions_in_subsystem = num_reactions_in_subsystem(find(non_empty_idx),:);
% Create the RxnSubsystems table
summary_subsystems = table(unique_subsystems, num_reactions_in_subsystem, list_reactions_in_subsystem, ...
    'VariableNames', {'SubsystemName', 'NumReactions', 'ListReactions'});

%% Summary Substrates
% Initialize cells to store unique substrates and their reactions
unique_substrates = {};
num_reactions_per_substrate = [];
list_reactions_per_substrate = {};

% Iterate over the relevant substrates
for i = 1:length(relevant_substrates)
    % Split the substrates by ';'
    substrates = strsplit(relevant_substrates{i}, '; ');
    for j = 1:length(substrates)
        substrate = substrates{j};
        % Find the index of the substrate in the list of unique substrates
        idx = find(strcmp(unique_substrates, substrate));
        if isempty(idx)
            % If the substrate is not in the list, add it
            unique_substrates{end+1} = substrate;
            num_reactions_per_substrate(end+1) = 1;
            list_reactions_per_substrate{end+1} = relevant_rxn_id{i};
        else
            % If the substrate is already in the list, update the counters
            num_reactions_per_substrate(idx) = num_reactions_per_substrate(idx) + 1;
            list_reactions_per_substrate{idx} = [list_reactions_per_substrate{idx}, ', ', relevant_rxn_id{i}];
        end
    end
end
% Create the Substrates table
summary_substrates = table(unique_substrates', num_reactions_per_substrate', list_reactions_per_substrate', ...
    'VariableNames', {'SubstrateName', 'NumReactions', 'ListReactions'});


%% Summary Products
% Initialize cells to store unique products and their reactions
unique_products = {};
num_reactions_per_product = [];
list_reactions_per_product = {};

% Iterate over the relevant products
for i = 1:length(relevant_products)
    % Split the products by ';'
    products = strsplit(relevant_products{i}, '; ');
    for j = 1:length(products)
        product = products{j};
        % Find the index of the product in the list of unique products
        idx = find(strcmp(unique_products, product));
        if isempty(idx)
            % If the product is not in the list, add it
            unique_products{end+1} = product;
            num_reactions_per_product(end+1) = 1;
            list_reactions_per_product{end+1} = relevant_rxn_id{i};
        else
            % If the product is already in the list, update the counters
            num_reactions_per_product(idx) = num_reactions_per_product(idx) + 1;
            list_reactions_per_product{idx} = [list_reactions_per_product{idx}, ', ', relevant_rxn_id{i}];
        end
    end
end
% Create the Products table
summary_products = table(unique_products', num_reactions_per_product', list_reactions_per_product', ...
    'VariableNames', {'productName', 'NumReactions', 'ListReactions'});


%% Summary Metabolites
% Initialize cells to store unique metabolites and their reactions
unique_metabolites = {};
num_reactions_per_metabolite = [];
list_reactions_per_metabolite = {};

% Use a container to avoid duplicates
metabolite_reaction_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Iterate over the relevant substrates and products
for i = 1:length(relevant_substrates)
    % Split the substrates by ';'
    substrates = strsplit(relevant_substrates{i}, '; ');
    for j = 1:length(substrates)
        metabolite = substrates{j};
        if isKey(metabolite_reaction_map, metabolite)
            % If the metabolite is already in the map, update the counters
            reactions_set = metabolite_reaction_map(metabolite);
            reactions_set = [reactions_set, relevant_rxn_id{i}];
            metabolite_reaction_map(metabolite) = unique(reactions_set);
        else
            % If the metabolite is not in the map, add it
            metabolite_reaction_map(metabolite) = {relevant_rxn_id{i}};
        end
    end
end
for i = 1:length(relevant_products)
    % Split the products by ';'
    products = strsplit(relevant_products{i}, '; ');
    for j = 1:length(products)
        metabolite = products{j};
        if isKey(metabolite_reaction_map, metabolite)
            % If the metabolite is already in the map, update the counters
            reactions_set = metabolite_reaction_map(metabolite);
            reactions_set = [reactions_set, relevant_rxn_id{i}];
            metabolite_reaction_map(metabolite) = unique(reactions_set);
        else
            % If the metabolite is not in the map, add it
            metabolite_reaction_map(metabolite) = {relevant_rxn_id{i}};
        end
    end
end
% Convert the map to cells to create the table
unique_metabolites = keys(metabolite_reaction_map);
num_reactions_per_metabolite = cellfun(@length, values(metabolite_reaction_map));
list_reactions_per_metabolite = cellfun(@(x) strjoin(x, ', '), values(metabolite_reaction_map), 'UniformOutput', false);
% Create the Metabolites table
summary_metabolites = table(unique_metabolites', num_reactions_per_metabolite', list_reactions_per_metabolite', ...
    'VariableNames', {'MetaboliteName', 'NumReactions', 'ListReactions'});


%% Summary Genes 
% Initialize cells to store unique genes and their reactions
unique_genes = {};
num_reactions_per_gene = [];
list_reactions_per_gene = {};

% Use a container to avoid duplicates
gene_reaction_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Iterate over the relevant genes
for i = 1:length(relevant_reaction_ncbi_genes)
    % Split the genes by ', '
    genes = strsplit(relevant_reaction_ncbi_genes{i}, ', ');
    for j = 1:length(genes)
        gene = genes{j};
        if isKey(gene_reaction_map, gene)
            % If the gene is already in the map, update the counters
            reactions_set = gene_reaction_map(gene);
            reactions_set = [reactions_set, relevant_rxn_id{i}];
            gene_reaction_map(gene) = unique(reactions_set);
        else
            % If the gene is not in the map, add it
            gene_reaction_map(gene) = {relevant_rxn_id{i}};
        end
    end
end

% Convert the map to cells to create the table
unique_genes = keys(gene_reaction_map);
num_reactions_per_gene = cellfun(@length, values(gene_reaction_map));
list_reactions_per_gene = cellfun(@(x) strjoin(x, ', '), values(gene_reaction_map), 'UniformOutput', false);

% Create the Genes table
summary_genes = table(unique_genes', num_reactions_per_gene', list_reactions_per_gene', ...
    'VariableNames', {'GeneName', 'NumReactions', 'ListReactions'});

% Sanitize and harmonize tables
relevant_reactions2 = table2cell(array2table(relevant_reactions'));
relevant_reaction_entrez_genes = relevant_reaction_entrez_genes';
relevant_reaction_ncbi_genes = relevant_reaction_ncbi_genes';
relevant_reactions2 = relevant_reactions2(:);
relevant_reactions2 = table2cell(array2table(relevant_reactions2'));
relevant_rxn_id = relevant_rxn_id(:);
relevant_rxn_name = relevant_rxn_name(:);
relevant_rxn_EC = relevant_rxn_EC(:);
relevant_rxn_subsystem = relevant_rxn_subsystem(:);
relevant_rxn_gpr = relevant_rxn_gpr(:);
relevant_reaction_entrez_genes = relevant_reaction_entrez_genes(:);
relevant_reaction_ncbi_genes = relevant_reaction_ncbi_genes(:);
relevant_substrates = relevant_substrates(:);
relevant_products = relevant_products(:);

%% Overall summary
summary_table = table(relevant_reactions2', relevant_rxn_id, relevant_rxn_name, relevant_rxn_EC, relevant_rxn_subsystem, relevant_rxn_gpr, relevant_reaction_entrez_genes, relevant_reaction_ncbi_genes, relevant_substrates, relevant_products, ...
    'VariableNames', {'Reactions', 'RxnID', 'RxnName', 'RxnEC', 'RxnSubsystem', 'RxnGPR', 'EntrezGenes', 'NCBIGenes', 'Substrates', 'Products'});

save('shap_subset_working_environment_2.mat', '-v7.3');

%% Write Excel
filename = 'MM-ML_results_summary_shap_subset.xlsx'; % Save the tables in an Excel file with different sheets
writetable(summary_table, filename, 'Sheet', 'Summary'); % Save the summary_table in the 'Summary' sheet
summary_subsystems_split = splitvars(summary_subsystems); % Split nested variables in summary_subsystems
writetable(summary_subsystems_split, filename, 'Sheet', 'Summary_pathways'); % Save the subsystems table
writetable(summary_genes, filename, 'Sheet', 'Summary_genes'); % Save the genes table
writetable(summary_metabolites, filename, 'Sheet', 'Summary_metabolites'); % Save the metabolites table
writetable(summary_substrates, filename, 'Sheet', 'Summary_substrates'); % Save the substrates table
writetable(summary_products, filename, 'Sheet', 'Summary_products'); % Save the products table

%% Check if the relevant reactions correlate with mortality rate in trauma group
folder = '/';

% List all files starting with "modelPatient_Sampled_"
files = dir(fullfile(folder, 'modelPatient_Sampled_*.mat'));
patient_data = table2array(readtable('Patient_Trauma_Groups.xls'));

female_label = 0;
male_label = 1;
label = female_label;

sex_idx = find(patient_data(:,4) == label);
patient_data = patient_data(sex_idx,:);
indices = []; % Initialize an empty vector to store the corresponding indices
% Loop through the vector of sample numbers
for i = 1:length(sex_idx)
    sample = sex_idx(i);    
    % Calculate the corresponding indices for this sample
    idx = (3*sample - 2):(3*sample);  % This gives [3n-2, 3n-1, 3n]    
    % Append the indices to the final vector
    indices = [indices, idx];
end
files = files(indices);

% Initialize cells to store file names and matrices
file_names = cell(length(files), 1);
points_matrices = cell(length(files), 8);

% Iterate over the files and load them
for i = 1:length(files)    
    file_path = fullfile(folder, files(i).name); % Construct the full path to the file  
    data = load(file_path); % Load the .mat file   
    points_matrix = data.sampleMetaOutC.points; % Extract the points matrix from sampleMetaOutC
    points_matrix_filtered = points_matrix(relevant_reactions,:);    
    clear points_matrix
    clear data
    % Extract parts of the file name
    pattern = 'modelPatient_Sampled_(\d+)_([a-zA-Z]+)';
    tokens = regexp(files(i).name, pattern, 'tokens');
    patient = str2num(tokens{1}{1});
    minmeanmax = tokens{1}{2};
    if strcmp(minmeanmax, 'Min')
        column = 6;
    elseif  strcmp(minmeanmax, 'Mean')
        column = 7;
    elseif  strcmp(minmeanmax, 'Max')
        column = 8;
    end
    patient_index_in_table_data = find(ismember(patient_data(:,1),patient));
    metabogroup = patient_data(patient_index_in_table_data,2);
    age = patient_data(patient_index_in_table_data,3);
    sex = patient_data(patient_index_in_table_data,4); 
    % Save the file name and points matrix in the cells
    points_matrices{patient,1} = patient;
    points_matrices{patient,2} = metabogroup;
    points_matrices{patient,3} = age;
    points_matrices{patient,4} = sex;
    points_matrices{patient,5} = files(i).name;
    points_matrices{patient,column} = points_matrix_filtered;
    disp(i)
end

points_matrices = points_matrices(~cellfun('isempty', points_matrices(:,1)), :);

for i = 1:length(points_matrices)
    points_matrices{i,9} = [points_matrices{i,6},points_matrices{i,7},points_matrices{i,8}];
end

% Graph representation and clustering
% Generate unique matrix with all solutions and the corresponding labels
total_points = size(points_matrices{1, 9}, 2);

% Calculate the sum of Euclidean distances for each column
num_labels = length(files)/3;
euclidean_sum = zeros(num_labels, total_points);
index_sum = euclidean_sum;
for i = 1:num_labels
    ith_point_matrices = points_matrices{i, 9};
    ith_euclidean_sum = zeros(1, total_points);    
    parfor j = 1:total_points        
        for k = 1:total_points
            if j ~= k
                ith_euclidean_sum(j) = euclidean_sum(j) + norm(ith_point_matrices(:, j) - ith_point_matrices(:, k));
            end
        end
    end
    ith_euclidean_sum = ith_euclidean_sum/max(ith_euclidean_sum);
    [ith_euclidean_sum, sortIdx] = sort(ith_euclidean_sum(1,:));
    euclidean_sum(i,:) = ith_euclidean_sum;
    index_sum(i,:) = sortIdx;
    clear ith_euclidean_sum 
    clear sortIdx
end

save('shap_subset_working_environment.mat', '-v7.3');


% Evaluate the results of the RESOS algorithm
indexes = 1:total_points;
euclidean_sum = sum(euclidean_sum);
to_plot = [indexes;euclidean_sum];

x = to_plot(1, :); % Extract the first row for x-axis values
y = to_plot(2, :); % Extract the second row for y-axis values
dy = diff(y) ./ diff(x); % Calculate the first derivative
d2y = diff(dy) ./ diff(x(1:end-1)); % Calculate the second derivative
[~, change_idx] = max(abs(d2y)); % Find the index of the maximum absolute change in the second derivative

% Plot the data
figure;
plot(x, y, '-o'); % Original data
hold on;
plot(x(change_idx + 1), y(change_idx + 1), 'r*', 'MarkerSize', 10); % Inflection point
title('Data Plot with Most Significant Inflection Point');
xlabel('X-axis');
ylabel('Y-axis');
grid on; % Adds a grid to the plot
legend('Data', 'Most Significant Inflection Point');

points_matrices_test = points_matrices;

for i = 1:num_labels
    points_matrices_test{i, 10} = points_matrices{i, 9}(:,index_sum(i,1:change_idx));   
end

result_matrix = [];
for i = 1:num_labels
   ith_label = strcat(string(points_matrices{i,1}),"-", string(points_matrices{i,2}));
   start_idx = (i-1) * change_idx + 1;
   end_idx = i * change_idx;
   labels(start_idx:end_idx) = {sprintf('%s', ith_label)};
   selected_points = points_matrices_test{i, 10};
   result_matrix = [result_matrix, selected_points]; % Concatenate the matrix 
end

%% MIN
total_points_min = size(points_matrices{1, 6}, 2);
num_labels = length(files)/3;
euclidean_sum_min = zeros(num_labels, total_points_min);
index_sum_min = euclidean_sum_min;
for i = 1:num_labels
    ith_point_matrices_min = points_matrices{i, 6};
    ith_euclidean_sum_min = zeros(1, total_points_min);    
    parfor j = 1:total_points_min        
        for k = 1:total_points_min
            if j ~= k
                ith_euclidean_sum_min(j) = euclidean_sum_min(j) + norm(ith_point_matrices_min(:, j) - ith_point_matrices_min(:, k));
            end
        end
    end
    ith_euclidean_sum_min = ith_euclidean_sum_min/max(ith_euclidean_sum_min);
    [ith_euclidean_sum_min, sortIdx] = sort(ith_euclidean_sum_min(1,:));
    euclidean_sum_min(i,:) = ith_euclidean_sum_min;
    index_sum_min(i,:) = sortIdx;
    clear ith_euclidean_sum 
    clear sortIdx
end
indexes_min = 1:total_points_min;
euclidean_sum_min = sum(euclidean_sum_min);
to_plot_min = [indexes_min;euclidean_sum_min];

x = to_plot_min(1, :); % Extract the first row for x-axis values
y = to_plot_min(2, :); % Extract the second row for y-axis values
dy = diff(y) ./ diff(x); % Calculate the first derivative
d2y = diff(dy) ./ diff(x(1:end-1)); % Calculate the second derivative
[~, change_idx_min] = min(abs(d2y)); % Find the index of the maximum absolute change in the second derivative

% Plot the data
figure;
plot(x, y, '-o'); % Original data
hold on;
plot(x(change_idx_min + 1), y(change_idx_min + 1), 'r*', 'MarkerSize', 10); % Inflection point
title('Data Plot with Most Significant Inflection Point min');
xlabel('X-axis');
ylabel('Y-axis');
grid on; % Adds a grid to the plot
legend('Data', 'Most Significant Inflection Point');

%% MEAN
total_points_mean = size(points_matrices{1, 7}, 2);
num_labels = length(files)/3;
euclidean_sum_mean = zeros(num_labels, total_points_mean);
index_sum_mean = euclidean_sum_mean;
for i = 1:num_labels
    ith_point_matrices_mean = points_matrices{i, 7};
    ith_euclidean_sum_mean = zeros(1, total_points_mean);    
    parfor j = 1:total_points_mean        
        for k = 1:total_points_mean
            if j ~= k
                ith_euclidean_sum_mean(j) = euclidean_sum_mean(j) + norm(ith_point_matrices_mean(:, j) - ith_point_matrices_mean(:, k));
            end
        end
    end
    ith_euclidean_sum_mean = ith_euclidean_sum_mean/max(ith_euclidean_sum_mean);
    [ith_euclidean_sum_mean, sortIdx] = sort(ith_euclidean_sum_mean(1,:));
    euclidean_sum_mean(i,:) = ith_euclidean_sum_mean;
    index_sum_mean(i,:) = sortIdx;
    clear ith_euclidean_sum 
    clear sortIdx
end
indexes_mean = 1:total_points_mean;
euclidean_sum_mean = sum(euclidean_sum_mean);
to_plot_mean = [indexes_mean;euclidean_sum_mean];

x = to_plot_mean(1, :); % Extract the first row for x-axis values
y = to_plot_mean(2, :); % Extract the second row for y-axis values
dy = diff(y) ./ diff(x); % Calculate the first derivative
d2y = diff(dy) ./ diff(x(1:end-1)); % Calculate the second derivative
[~, change_idx_mean] = min(abs(d2y)); % Find the index of the maximum absolute change in the second derivative

% Plot the data
figure;
plot(x, y, '-o'); % Original data
hold on;
plot(x(change_idx_mean + 1), y(change_idx_mean + 1), 'r*', 'MarkerSize', 10); % Inflection point
title('Data Plot with Most Significant Inflection Point mean');
xlabel('X-axis');
ylabel('Y-axis');
grid on; % Adds a grid to the plot
legend('Data', 'Most Significant Inflection Point');

%% MAX
total_points_max = size(points_matrices{1, 8}, 2);
num_labels = length(files)/3;
euclidean_sum_max = zeros(num_labels, total_points_max);
index_sum_max = euclidean_sum_max;
for i = 1:num_labels
    ith_point_matrices_max = points_matrices{i, 8};
    ith_euclidean_sum_max = zeros(1, total_points_max);    
    parfor j = 1:total_points_max        
        for k = 1:total_points_max
            if j ~= k
                ith_euclidean_sum_max(j) = euclidean_sum_max(j) + norm(ith_point_matrices_max(:, j) - ith_point_matrices_max(:, k));
            end
        end
    end
    ith_euclidean_sum_max = ith_euclidean_sum_max/max(ith_euclidean_sum_max);
    [ith_euclidean_sum_max, sortIdx] = sort(ith_euclidean_sum_max(1,:));
    euclidean_sum_max(i,:) = ith_euclidean_sum_max;
    index_sum_max(i,:) = sortIdx;
    clear ith_euclidean_sum 
    clear sortIdx
end
indexes_max = 1:total_points_max;
euclidean_sum_max = sum(euclidean_sum_max);
to_plot_max = [indexes_max;euclidean_sum_max];

x = to_plot_max(1, :); % Extract the first row for x-axis values
y = to_plot_max(2, :); % Extract the second row for y-axis values
dy = diff(y) ./ diff(x); % Calculate the first derivative
d2y = diff(dy) ./ diff(x(1:end-1)); % Calculate the second derivative
[~, change_idx_max] = min(abs(d2y)); % Find the index of the maximum absolute change in the second derivative

% Plot the data
figure;
plot(x, y, '-o'); % Original data
hold on;
plot(x(change_idx_max + 1), y(change_idx_max + 1), 'r*', 'MarkerSize', 10); % Inflection point
title('Data Plot with Most Significant Inflection Point max');
xlabel('X-axis');
ylabel('Y-axis');
grid on; % Adds a grid to the plot
legend('Data', 'Most Significant Inflection Point');

% Plot all points
points_matrices_test_2 = points_matrices_test;
for i = 1:num_labels
    points_matrices_test_2{i, 11} = [points_matrices{i, 6}(:,index_sum_min(i,1:change_idx_min)), points_matrices{i, 7}(:,index_sum_mean(i,1:change_idx_mean)), points_matrices{i, 8}(:,index_sum_max(i,1:change_idx_max))];   
end
total_points_all = size(points_matrices_test_2{1, 11}, 2);
num_labels = length(files)/3;
euclidean_sum_all = zeros(num_labels, total_points_all);
index_sum_all = euclidean_sum_all;
for i = 1:num_labels
    ith_point_matrices_all = points_matrices_test_2{i, 11};
    ith_euclidean_sum_all = zeros(1, total_points_all);    
    parfor j = 1:total_points_all        
        for k = 1:total_points_all
            if j ~= k
                ith_euclidean_sum_all(j) = euclidean_sum_all(j) + norm(ith_point_matrices_all(:, j) - ith_point_matrices_all(:, k));
            end
        end
    end
    ith_euclidean_sum_all = ith_euclidean_sum_all/max(ith_euclidean_sum_all);
    [ith_euclidean_sum_all, sortIdx] = sort(ith_euclidean_sum_all(1,:));
    euclidean_sum_all(i,:) = ith_euclidean_sum_all;
    index_sum_all(i,:) = sortIdx;
    clear ith_euclidean_sum_all 
    clear sortIdx
end
indexes_all = 1:total_points_all;
euclidean_sum_all = sum(euclidean_sum_all);
to_plot_all = [indexes_all;euclidean_sum_all];

x = to_plot_all(1, :); % Extract the first row for x-axis values
y = to_plot_all(2, :); % Extract the second row for y-axis values
dy = diff(y) ./ diff(x); % Calculate the first derivative
d2y = diff(dy) ./ diff(x(1:end-1)); % Calculate the second derivative
[~, change_idx_all] = max(abs(d2y)); % Find the index of the maximum absolute change in the second derivative

% Plot the data
figure;
plot(x, y, '-o'); % Original data
hold on;
plot(x(change_idx_all + 1), y(change_idx_all + 1), 'r*', 'MarkerSize', 10); % Inflection point
title('Data Plot with Most Significant Inflection Point all');
xlabel('X-axis');
ylabel('Y-axis');
grid on; % Adds a grid to the plot
legend('Data', 'Most Significant Inflection Point');

points_matrices_test_3 = points_matrices_test_2;
for i = 1:num_labels
    points_matrices_test_3{i, 12} = points_matrices_test_2{i, 9}(:,index_sum_min(i,1:change_idx_all));   
end

result_matrix_all = [];
points_per_sample = length(points_matrices_test_3{1, 12});
for i = 1:num_labels
   ith_label = strcat(string(points_matrices_test_3{i,1}),"-", string(points_matrices_test_3{i,4}));
   start_idx = (i-1) * points_per_sample + 1;
   end_idx = i * points_per_sample;
   labels(start_idx:end_idx) = {sprintf('%s', ith_label)};
   selected_points = points_matrices_test_3{i, 12};
   result_matrix_all = [result_matrix_all, selected_points]; % Concatenate the matrix 
end

% Plot all points mean of min, mean and max
points_matrices_test_4 = points_matrices_test_3;
for i = 1:num_labels
    mean_min = mean(points_matrices{i, 6}(:,index_sum_min(i,1:change_idx_min)),2);
    mean_mean = mean(points_matrices{i, 7}(:,index_sum_mean(i,1:change_idx_mean)),2);
    mean_max = mean(points_matrices{i, 8}(:,index_sum_max(i,1:change_idx_max)),2);
    points_matrices_test_4{i, 13} = [mean_min, mean_mean, mean_max];
end

result_matrix_all_mean = [];
points_per_sample_all_mean = size(points_matrices_test_4{1, 13},2); %length(points_matrices_test_4{1, 13});
for i = 1:num_labels
   ith_label = strcat(string(points_matrices_test_4{i,1}),"-", string(points_matrices_test_4{i,4}));
   start_idx = (i-1) * points_per_sample_all_mean + 1;
   end_idx = i * points_per_sample_all_mean;
   labels_all_mean(start_idx:end_idx) = {sprintf('%s', ith_label)};
   selected_points = points_matrices_test_4{i, 13};
   result_matrix_all_mean = [result_matrix_all_mean, selected_points]; % Concatenate the matrix 
end

points_matrices_test_5 = points_matrices_test_4;
for i = 1:num_labels
    mean_min = mean(points_matrices{i, 6}(:,index_sum_min(i,1:change_idx_min)),2);
    mean_mean = mean(points_matrices{i, 7}(:,index_sum_mean(i,1:change_idx_mean)),2);
    mean_max = mean(points_matrices{i, 8}(:,index_sum_max(i,1:change_idx_max)),2);
    points_matrices_test_5{i, 14} = [mean([mean_min, mean_mean, mean_max],2)];
end

result_matrix_all_mean_mean = [];
points_per_sample_all_mean_mean = size(points_matrices_test_5{1, 14},2); %length(points_matrices_test_4{1, 13});
for i = 1:num_labels
   ith_label = strcat(string(points_matrices_test_5{i,1}),"-", string(points_matrices_test_5{i,2}),"-", string(points_matrices_test_5{i,4}));
   start_idx = (i-1) * points_per_sample_all_mean_mean + 1;
   end_idx = i * points_per_sample_all_mean_mean;
   labels_all_mean_mean(start_idx:end_idx) = {sprintf('%s', ith_label)};
   selected_points = points_matrices_test_5{i, 14};
   result_matrix_all_mean_mean = [result_matrix_all_mean_mean, selected_points]; % Concatenate the matrix 
end

save('shap_subset_working_environment_2.mat', '-v7.3');

%% Save variables in files to be used in Python for example
csvwrite('result_matrix.csv', result_matrix); % Export the resulting matrix to a CSV file
csvwrite('result_matrix_all.csv', result_matrix_all); % Export the resulting matrix to a CSV file
flattenedArray = labels(:); % Flatten the cell array to a single column
labels_all_mean_2 = string(labels_all_mean(:));
labels_all_mean_3 = string(labels_all_mean_mean(:));
cell2csv('patietn_labels.csv', flattenedArray); % Write the flattened array to a CSV file
cell2csv('patietn_labels_all_mean.csv', labels_all_mean_2); % Write the flattened array to a CSV file
cell2csv('patietn_labels_all_mean_mean.csv', labels_all_mean_3); % Write the flattened array to a CSV file
writematrix(relevant_reactions', 'relevant_reactions_labels.csv');                     
%save(yourfile_male.mat, 'flattenedArray', 'labels', 'labels_all_mean', 'labels_all_mean_mean', 'points_matrices_test_4','relevant_reactions', 'result_matrix_all', 'result_matrix_all_mean', '-v7.3');

%% Plot
column_labels = cellstr(labels_all_mean_3);
cgo=clustergram(result_matrix_all_mean_mean,'Standardize','Row', 'Colormap', 'redbluecmap');
%set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_labels, 'ColumnLabelsRotate', 45)
set(cgo,'Linkage','complete','ColumnLabels', column_labels, 'ColumnLabelsRotate', 45)
%rm = struct('GroupNumber',{180,185,183,184,182},'Annotation',{'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5'},'Color',{'m','r','g','y','b'});
%set(cgo,'RowGroupMarker',rm)

%% Analyze endo-type abundance in top_44 and subset_7 reactions clusters
% top_44
clusters = 1:4; % x-axis representing Cluster 1, 2, 3, and 4
groupA = [0.2500, 0.0000, 0.0000, 0.7500];
groupB = [0.5000, 0.1667, 0.0833, 0.2500];
groupC = [0.1724, 0.2759, 0.3448, 0.2069];
groupD = [0.0476, 0.3333, 0.4762, 0.1429];

% Plot each group as a line
figure;
hold on; % Allow multiple lines on the same plot
plot(clusters, groupA, '-o', 'DisplayName', 'Endo- Group A', 'LineWidth', 1.5);
plot(clusters, groupB, '-s', 'DisplayName', 'Endo- Group B', 'LineWidth', 1.5);
plot(clusters, groupC, '-d', 'DisplayName', 'Endo- Group C', 'LineWidth', 1.5);
plot(clusters, groupD, '-^', 'DisplayName', 'Endo- Group D', 'LineWidth', 1.5);

% Add title and labels
title('Cluster analysis - Cluster');
xlabel('Cluster');
ylabel('Proportion');
legend('show'); % Display the legend
grid on; % Add grid for better readability
hold off;
 
% subset_7
clusters = 1:4; % x-axis representing Cluster 1, 2, 3, and 4
groupA = [0.0000, 0.2500, 0.5833, 0.1667];
groupB = [0.1667, 0.0833, 0.6667, 0.0833];
groupC = [0.2759, 0.3448, 0.2759, 0.1034];
groupD = [0.5476, 0.1905, 0.0238, 0.2381];

% Plot each group as a line
figure;
hold on; % Allow multiple lines on the same plot
plot(clusters, groupA, '-o', 'DisplayName', 'Endo- Group A', 'LineWidth', 1.5);
plot(clusters, groupB, '-s', 'DisplayName', 'Endo- Group B', 'LineWidth', 1.5);
plot(clusters, groupC, '-d', 'DisplayName', 'Endo- Group C', 'LineWidth', 1.5);
plot(clusters, groupD, '-^', 'DisplayName', 'Endo- Group D', 'LineWidth', 1.5);

% Add title and labels
title('Cluster analysis - Cluster');
xlabel('Cluster');
ylabel('Proportion');
legend('show'); % Display the legend
grid on; % Add grid for better readability
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Treatment simulation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine the reactions of interest
relevant_intracel_reactions = [7, 22, 23, 34, 2, 19, 29];

group_1 = points_matrices_test_5(find([points_matrices_test_5{:,2}] == 1),:);
[n_patient_1, ~] =  size(group_1);

max_val = [];
min_val = [];
for i = 1:n_patient_1
    ith_patient_flux = group_1{i,12};
    ith_patient_flux_relevant = ith_patient_flux(relevant_intracel_reactions,:);
    ith_max = max(ith_patient_flux_relevant, [], 2);
    ith_min = min(ith_patient_flux_relevant,[], 2);
    max_val = [max_val, ith_max];
    min_val = [min_val, ith_min];
end

%max_val = max(max_val, [], 2);
%min_val = min(min_val, [], 2);
max_val = mean(max_val, 2);
min_val = mean(min_val, 2);

group_4 = points_matrices_test_5(find([points_matrices_test_5{:,2}] == 4),:);
[n_patient_4, ~] =  size(group_4);

max_val_4 = [];
min_val_4 = [];
for i = 1:n_patient_4
    ith_patient_flux = group_4{i,12};
    ith_patient_flux_relevant = ith_patient_flux(relevant_intracel_reactions,:);
    ith_max = max(ith_patient_flux_relevant, [], 2);
    ith_min = min(ith_patient_flux_relevant,[], 2);
    max_val_4 = [max_val_4, ith_max];
    min_val_4 = [min_val_4, ith_min];
end

%max_val_4 = max(max_val_4, [], 2);
%min_val_4 = min(min_val_4, [], 2);
max_val_4 = mean(max_val_4, 2);
min_val_4 = mean(min_val_4, 2);