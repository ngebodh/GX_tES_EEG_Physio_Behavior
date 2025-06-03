%% GX_DataDownSample_CheckFiles
%% Downsampling Test Script
%
% This script is designed to test the output of downsampling
% and data extraction, specifically for EEG and ptracker data.
% It compares the downsampled and extracted data against a re-processed
% version of the original raw data to ensure data integrity and accuracy.
%
% The script performs the following checks for each dataset found:
% 1.  **EEG Time-Series Data Comparison:**
%     - Replicates the padding and downsampling logic from the main
%       downsampling script (`GX_DataDownSample.m`) on the original EEG data.
%     - Calculates Root Mean Squared Error (RMSE) and Mean Absolute Error (MAE)
%       between the re-processed original EEG channel and the downsampled EEG data.
%     - Logs pass/fail based on predefined RMSE and MAE thresholds.
% 2.  **Trigger Data Comparison:**
%     - Compares the total number of triggers.
%     - Compares trigger codes and types.
%     - Compares trigger offsets (timestamps) within a defined tolerance,
%       accounting for changes due to downsampling.
%     - Logs pass/fail for each trigger aspect.
% 3.  **Ptracker Data Comparison:**
%     - Replicates the import and resampling logic for the original ptracker CSV data.
%     - Calculates RMSE and MAE between the re-processed original ptracker
%       performance data and the downsampled ptracker performance data.
%     - Logs pass/fail based on predefined RMSE and MAE thresholds.
%
% INPUTS:
%   - `original_data_folder`: Path to the root directory containing original
%     EEG (.cnt) and ptracker (.csv) files.
%     Expected structure: `original_data_folder\dataset_id\dataset_id\ptracker-*.csv`
%     and `original_data_folder\dataset_id\*.cnt`
%   - `downsampled_data_folder`: Path to the directory containing the
%     downsampled MATLAB structures (`EEG_DS_Struct_XXXX.mat`).
%   - `EEG_RMSE_THRESHOLD`: Acceptable RMSE for EEG signal comparison (e.g., 0.01).
%   - `EEG_MAE_THRESHOLD`: Acceptable MAE for EEG signal comparison (e.g., 0.005).
%   - `PTRACKER_RMSE_THRESHOLD`: Acceptable RMSE for ptracker signal comparison (e.g., 0.05).
%   - `PTRACKER_MAE_THRESHOLD`: Acceptable MAE for ptracker signal comparison (e.g., 0.02).
%   - `TRIGGER_OFFSET_TOLERANCE`: Maximum allowed difference in trigger
%     sample offsets (e.g., 1 for +/- 1 sample due to rounding).
%
% OUTPUTS:
%   - A log file (`Downsampling_Test_Log_YYYYMMDD_HHMMSS.txt`) generated
%     in the current working directory, detailing the test results for
%     each dataset, including pass/fail status and quantitative metrics.
%
% DEPENDENCIES:
%   - `read_eep_cnt` function (from ANT EEG import functions).
%   - MATLAB's Signal Processing Toolbox (`resample` function).
%
% USAGE:
%   1. Set the `original_data_folder` and `downsampled_data_folder` variables
%      to your specific file paths.
%   2. Adjust the `_THRESHOLD` and `_TOLERANCE` variables as needed for your
%      data's characteristics and acceptable error margins.
%   3. Run the script in MATLAB.
%   4. Review the generated log file for test results.
%
% Author: Nigel Gebodh
% Date: May 28, 2025

%% Start
% Clear variables from previous runs, except for the paths if they were defined externally.
clearvars -except original_data_folder downsampled_data_folder

% --- User Inputs and Setup ---
% Define paths for original and downsampled data.
% IMPORTANT: Ensure these paths are correct for your system.
original_data_folder = 'F:\GX\Data\'; % Path to your original .cnt and .evt files
downsampled_data_folder = 'G:\GX_Dataset_DS_V3\Data_downsampled_05282025\'; % Path to your downsampled .mat files

% Define thresholds for passing tests. Adjust these values based on your data
% characteristics and acceptable error margins.
EEG_RMSE_THRESHOLD = 0.01;   % Root Mean Squared Error threshold for EEG signal comparison
EEG_MAE_THRESHOLD = 0.005;   % Mean Absolute Error threshold for EEG signal comparison
PTRACKER_RMSE_THRESHOLD = 0.05; % RMSE threshold for ptracker signal comparison
PTRACKER_MAE_THRESHOLD = 0.02; % MAE threshold for ptracker signal comparison
TRIGGER_OFFSET_TOLERANCE = 1; % Maximum allowed difference in trigger offset (samples).
                              % Set to 0 for exact match, 1 for allowing +/- 1 sample due to rounding.

% --- Setup Logger ---
% Create a unique log filename with a timestamp.
log_filename = sprintf('Downsampling_Test_Log_%s.txt', datestr(now, 'yyyymmdd_HHMMSS'));
% Define the full path for the log file (saves in the current working directory).
log_filepath = fullfile(pwd, log_filename);
% Open the log file for writing. 'w' mode overwrites existing files.
log_file_id = fopen(log_filepath, 'w');

% Check if the log file was opened successfully.
if log_file_id == -1
    error('Could not open log file for writing. Check permissions or path: %s', log_filepath);
end

% Write initial header information to the log file.
fprintf(log_file_id, '--- Downsampling Test Log ---\n');
fprintf(log_file_id, 'Date: %s\n', datestr(now));
fprintf(log_file_id, 'Original Data Folder: %s\n', original_data_folder);
fprintf(log_file_id, 'Downsampled Data Folder: %s\n', downsampled_data_folder);
fprintf(log_file_id, 'EEG RMSE Threshold: %.4f\n', EEG_RMSE_THRESHOLD);
fprintf(log_file_id, 'EEG MAE Threshold: %.4f\n', EEG_MAE_THRESHOLD);
fprintf(log_file_id, 'Ptracker RMSE Threshold: %.4f\n', PTRACKER_RMSE_THRESHOLD);
fprintf(log_file_id, 'Ptracker MAE Threshold: %.4f\n', PTRACKER_MAE_THRESHOLD);
fprintf(log_file_id, 'Trigger Offset Tolerance: %d samples\n', TRIGGER_OFFSET_TOLERANCE);
fprintf(log_file_id, '-----------------------------\n\n');

% --- Get list of all downsampled files to test ---
% Find all .mat files that follow the naming convention 'EEG_DS_Struct_XXXX.mat'.
downsampled_mat_files = dir(fullfile(downsampled_data_folder, 'EEG_DS_Struct_*.mat'));

% Check if any downsampled files were found.
if isempty(downsampled_mat_files)
    error('No downsampled .mat files found in the specified folder: %s. Please check the path and file naming convention.', downsampled_data_folder);
end

total_files = length(downsampled_mat_files);
fprintf('Found %d downsampled files to test.\n\n', total_files);
fprintf(log_file_id, 'Found %d downsampled files to test.\n\n', total_files);

% --- Loop through each downsampled file found ---
for file_idx = 1:total_files
    current_filename = downsampled_mat_files(file_idx).name;
    
    % Extract the dataset_id (e.g., '0101' from 'EEG_DS_Struct_0101.mat').
    % This assumes the ID is always between 'EEG_DS_Struct_' and '.mat'.
    name_parts = strsplit(current_filename, {'EEG_DS_Struct_', '.mat'});
    dataset_id = name_parts{2}; 

    fprintf('\n--- Testing Dataset: %s (%d/%d) ---\n', dataset_id, file_idx, total_files);
    fprintf(log_file_id, '\n--- Testing Dataset: %s (%d/%d) ---\n', dataset_id, file_idx, total_files);

    % Use a try-catch block to handle errors for individual datasets gracefully,
    % allowing the script to continue processing other files.
    try
        % --- Load Downsampled Data ---
        downsampled_file_path = fullfile(downsampled_data_folder, current_filename);
        fprintf('  Loading downsampled data from: %s\n', downsampled_file_path);
        fprintf(log_file_id, '  Loading downsampled data from: %s\n', downsampled_file_path);
        downsampled_data_struct = load(downsampled_file_path);
        DSamp = downsampled_data_struct.DSamp; % Access the DSamp structure

        % --- Load Original EEG Data (.cnt file) ---
        % Construct path for original .cnt file based on the assumed structure:
        % original_data_folder\dataset_id\GX_*.cnt
        original_eeg_path = fullfile(original_data_folder, dataset_id, filesep); % filesep handles / or \
        original_eeg_file_struct = dir(fullfile(original_eeg_path, '*.cnt'));

        if isempty(original_eeg_file_struct)
            error('Original .cnt file not found for dataset %s in %s. Skipping this dataset.', dataset_id, original_eeg_path);
        end
        original_eeg_filename = fullfile(original_eeg_path, original_eeg_file_struct(1).name);
        fprintf('  Loading original EEG data from: %s\n', original_eeg_filename);
        fprintf(log_file_id, '  Loading original EEG data from: %s\n', original_eeg_filename);
        
        % Read only 5 samples first to get nsample, then load the full file.
        Samp_temp = read_eep_cnt(original_eeg_filename, 1, 5);
        original_EEG = read_eep_cnt(original_eeg_filename, 1, Samp_temp.nsample); % Load entire original EEG
        
        % Handle the specific case for '1401' where a trigger is removed.
        if strcmp(dataset_id, '1401')
            if ~isempty(original_EEG.triggers) % Ensure triggers exist before attempting to remove
                original_EEG.triggers(1) = [];
                fprintf('  Note: Removed first trigger for dataset 1401 to match downsampling script behavior.\n');
                fprintf(log_file_id, '  Note: Removed first trigger for dataset 1401 to match downsampling script behavior.\n');
            else
                fprintf('  Warning: No triggers found for dataset 1401, cannot remove first trigger as per script logic.\n');
                fprintf(log_file_id, '  Warning: No triggers found for dataset 1401, cannot remove first trigger as per script logic.\n');
            end
        end

        % --- Load and Process Original Ptracker Data (.csv file) ---
        % Construct path for original ptracker-*.csv file based on the assumed structure:
        % original_data_folder\dataset_id\dataset_id\ptracker-*.csv
        original_ptracker_path = fullfile(original_data_folder, dataset_id, dataset_id, filesep);
        original_ptracker_file = fullfile(original_ptracker_path, ['ptracker-', dataset_id, '.csv']);

        if ~exist(original_ptracker_file, 'file')
            error('Original ptracker CSV file not found for dataset %s in %s. Skipping this dataset.', dataset_id, original_ptracker_path);
        end
        fprintf('  Loading original ptracker data from: %s\n', original_ptracker_file);
        fprintf(log_file_id, '  Loading original ptracker data from: %s\n', original_ptracker_file);

        % Replicate ptracker import logic from GX_DataDownSample.m script.
        delimiter = ',';
        startRow = 2;
        formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]'; % Matches original script's textscan format
        fileID = fopen(original_ptracker_file,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        
        original_ptracker_raw = [dataArray{[3,11]}]; % Extract time (column 3) and performance (column 11)
        % Normalize time: subtract first sample time and convert to seconds.
        original_ptracker_raw(:,1) = (original_ptracker_raw(:,1) - original_ptracker_raw(1,1)) ./ 1000;

        % Replicate ptracker resampling from GX_DataDownSample.m script.
        % The original script uses `resample(data, time, desiredFs, desiredFs, ScreenFs)`.
        % `desiredFs` is DSamp.ptrackerfs, and `ScreenFs` is hardcoded as 60.
        ScreenFs = 60; %Refresh rate of screen, hardcoded to account for dropped frames
        expected_ptracker_perf = resample(original_ptracker_raw(:,2), original_ptracker_raw(:,1), DSamp.ptrackerfs, DSamp.ptrackerfs, ScreenFs);
        expected_ptracker_time = [[0:length(expected_ptracker_perf)-1]./DSamp.ptrackerfs]';

        % --- Test EEG Time-Series Data ---
        fprintf(log_file_id, '\n  --- EEG Time-Series Test ---\n');
        
        channel_to_plot = 1; %Look at only the 1st channel

        fs_original_eeg = original_EEG.rate;
        fs_downsampled_eeg = DSamp.fs;
        DownSampleFactor_eeg = fs_original_eeg / fs_downsampled_eeg;

        % Replicate padding logic from your downsampling script for accurate comparison.
        
        padding_samples_eeg = fs_original_eeg; 
        DSpadremove_eeg = padding_samples_eeg / DownSampleFactor_eeg;

        original_channel_data = original_EEG.data(channel_to_plot, :);

        % Apply padding to the original channel data.
        padded_original_channel_data = [...
            ones(1, padding_samples_eeg) * original_channel_data(1), ... % Pad with first sample value
            original_channel_data, ...
            ones(1, padding_samples_eeg) * original_channel_data(end)    % Pad with last sample value
        ];

        % Resample the padded original EEG data to the downsampled rate.
        resampled_padded_data = resample(padded_original_channel_data, 1, DownSampleFactor_eeg);
        
        % Remove the padded segments from the resampled data, matching your script's logic.
        original_resampled_data_clean = resampled_padded_data(DSpadremove_eeg + 1 : end - DSpadremove_eeg);

        % Ensure both signals have the same length for direct comparison.
        % Trim to the minimum length if there are slight discrepancies.
        min_len_eeg = min(length(original_resampled_data_clean), length(DSamp.EEGdata(channel_to_plot,:)));
        signal_diff_original_eeg = original_resampled_data_clean(1:min_len_eeg);
        signal_diff_downsampled_eeg = DSamp.EEGdata(channel_to_plot, 1:min_len_eeg);

        % Calculate the difference between the two signals.
        eeg_signal_difference = signal_diff_original_eeg - signal_diff_downsampled_eeg;
        % Calculate Root Mean Squared Error (RMSE).
        eeg_rmse = sqrt(mean(eeg_signal_difference.^2));
        % Calculate Mean Absolute Error (MAE).
        eeg_mae = mean(abs(eeg_signal_difference));

        % Log the results for EEG time-series comparison.
        fprintf(log_file_id, '    EEG Channel %d RMSE: %.6f (Threshold: %.4f)\n', channel_to_plot, eeg_rmse, EEG_RMSE_THRESHOLD);
        fprintf(log_file_id, '    EEG Channel %d MAE: %.6f (Threshold: %.4f)\n', channel_to_plot, eeg_mae, EEG_MAE_THRESHOLD);

        if eeg_rmse < EEG_RMSE_THRESHOLD && eeg_mae < EEG_MAE_THRESHOLD
            fprintf(log_file_id, '    EEG Time-series data match: PASSED\n');
        else
            fprintf(log_file_id, '    EEG Time-series data match: FAILED\n');
        end

        % --- Test Trigger Data ---
        fprintf(log_file_id, '\n  --- Trigger Test ---\n');

        num_original_triggers = length(original_EEG.triggers);
        num_downsampled_triggers = length(DSamp.triggers);

        fprintf(log_file_id, '    Original Number of Triggers: %d\n', num_original_triggers);
        fprintf(log_file_id, '    Downsampled Number of Triggers: %d\n', num_downsampled_triggers);

        trigger_count_pass = (num_original_triggers == num_downsampled_triggers);
        if trigger_count_pass
            fprintf(log_file_id, '    Trigger count match: PASSED\n');
        else
            fprintf(log_file_id, '    Trigger count match: FAILED\n');
        end

        trigger_code_type_pass = true;
        % Only compare codes/types if the trigger counts match.
        if trigger_count_pass 
            for i = 1:num_original_triggers
                % Convert trigger codes to double for numerical comparison.
                original_code = str2double(original_EEG.triggers(i).code); 
                downsampled_code = str2double(DSamp.triggers(i).code); 
                original_type = original_EEG.triggers(i).type;
                downsampled_type = DSamp.triggers(i).type;

                if original_code ~= downsampled_code || original_type ~= downsampled_type
                    fprintf(log_file_id, '    Mismatch at trigger %d: Original Code=%d, Type=%d; Downsampled Code=%d, Type=%d\n', ...
                        i, original_code, original_type, downsampled_code, downsampled_type);
                    trigger_code_type_pass = false;
                    break; % Exit loop on first mismatch to avoid excessive logging for a single issue.
                end
            end
        else
            trigger_code_type_pass = false; % If counts don't match, full code/type match is impossible.
        end

        if trigger_code_type_pass
            fprintf(log_file_id, '    Trigger codes and types match: PASSED\n');
        else
            fprintf(log_file_id, '    Trigger codes and types match: FAILED\n');
        end

        trigger_offset_pass = true;
        % Only compare offsets if the trigger counts match.
        if trigger_count_pass 
            % Calculate expected downsampled offsets based on original offsets and downsampling factor.
            expected_downsampled_offsets = ceil([original_EEG.triggers.offset] * (fs_downsampled_eeg / fs_original_eeg));
            actual_downsampled_offsets = [DSamp.triggers.offset];

            for i = 1:num_original_triggers
                % Check if the absolute difference exceeds the defined tolerance.
                if abs(expected_downsampled_offsets(i) - actual_downsampled_offsets(i)) > TRIGGER_OFFSET_TOLERANCE
                    fprintf(log_file_id, '    Offset Mismatch at trigger %d: Expected=%d, Actual=%d (Tolerance: %d)\n', ...
                        i, expected_downsampled_offsets(i), actual_downsampled_offsets(i), TRIGGER_OFFSET_TOLERANCE);
                    trigger_offset_pass = false;
                    break; % Exit loop on first mismatch.
                end
            end
        else
            trigger_offset_pass = false; % If counts don't match, full offset match is impossible.
        end

        if trigger_offset_pass
            fprintf(log_file_id, '    Trigger offsets match (within tolerance): PASSED\n');
        else
            fprintf(log_file_id, '    Trigger offsets match: FAILED\n');
        end

        % --- Test Ptracker Data ---
        fprintf(log_file_id, '\n  --- Ptracker Data Test ---\n');

        % Ensure both ptracker signals have the same length for direct comparison.
        min_len_ptracker = min(length(expected_ptracker_perf), length(DSamp.ptrackerPerf));
        ptracker_diff_expected = expected_ptracker_perf(1:min_len_ptracker);
        ptracker_diff_actual = DSamp.ptrackerPerf(1:min_len_ptracker);

        % Calculate the difference, RMSE, and MAE for ptracker data.
        ptracker_signal_difference = ptracker_diff_expected - ptracker_diff_actual;
        ptracker_rmse = sqrt(mean(ptracker_signal_difference.^2));
        ptracker_mae = mean(abs(ptracker_signal_difference));

        % Log the results for ptracker data comparison.
        fprintf(log_file_id, '    Ptracker Performance RMSE: %.6f (Threshold: %.4f)\n', ptracker_rmse, PTRACKER_RMSE_THRESHOLD);
        fprintf(log_file_id, '    Ptracker Performance MAE: %.6f (Threshold: %.4f)\n', ptracker_mae, PTRACKER_MAE_THRESHOLD);

        if ptracker_rmse < PTRACKER_RMSE_THRESHOLD && ptracker_mae < PTRACKER_MAE_THRESHOLD
            fprintf(log_file_id, '    Ptracker data match: PASSED\n');
        else
            fprintf(log_file_id, '    Ptracker data match: FAILED\n');
        end

    catch ME % Catch any errors that occur during processing of a single dataset.
        fprintf(log_file_id, '  ERROR processing dataset %s: %s\n', dataset_id, ME.message);
        fprintf(log_file_id, '  Test for dataset %s: FAILED due to error.\n', dataset_id);
        warning('Error processing dataset %s: %s', dataset_id, ME.message); % Display warning in command window.
    end
    fprintf(log_file_id, '-------------------------------------\n');
end

% --- Finalize Logger ---
fprintf(log_file_id, '\n--- Test Summary ---\n');
fprintf(log_file_id, 'All tests completed. See log file for details.\n');
fclose(log_file_id); % Close the log file.
fprintf('\nAll tests completed. Results logged to: %s\n', log_filepath);


% --- END ---
