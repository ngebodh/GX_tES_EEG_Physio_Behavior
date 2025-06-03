%% GX_DataDownSample
%This script was written to downsample the EEG data to make the data more
%manageable. It converts the EEG data to MATLAB readable format and
%combines the CTT data with the EEG data. 
%
%INPUT:
%   -The location you want your data written to
%       Specified as a path i.e. D:\GX Project\Results\DataDownSampled\
%
%   -The location of your data. Note the matfiles should be in this folder too
%   inorder to write the stimulation type into the triggers. 
%       Specified as a path i.e. D:\GX Project\Data\
%
%   -The filename/numbers within the data folder
%
%   -Indicate if you want to load the default folders assigned in the script
%       or to use the ones you input above.
%
%   -The downsample factor. For example if my data are sampled at 2k Hz and I
%       want to downsample them to 250 Hz my downsample factor would be 8 (must be an integer) or 2000/250. 
%
%ASSUMED DATA FILE STRUCTURE:
%     |_Project
%      |----_Data
%           |----_0101
%                |----_0101
%                     |----ptracker-0101.csv
%                     |----ptracker-summary-0101..txt  <-(yes, double dot)
%                |----GX_01_2019-09-24_15-45-53.cnt
%                |----GX_01_2019-09-24_15-45-53.evt
%                |----MATLABfilestream0101924.mat
%                |----MATLABfilestream0101924.txt
%           |----0102
%           |----0103
%           :
%           :
%
%Requirements:
%-EEGLab
%-EEG import functions from ANT added to path:
%   https://www.ant-neuro.com/support/supporting-documentation-and-downloads
%
%Written by Nigel Gebodh in Oct. 2019.
%
%Updated:
% -05/28/2025:
% +Fixed bug on writing the ptracker timeseries out
% +Changed timers to display load times and elapsed
% +Changed mat file output to version 7.3
%
% -10/2/2023:
% +Fixed bug on writing montages to mat files.
% +Now default downsample mode takes all 4 digit numeric folder and
%  downsamples it when given path to data. 
%   
%
%

%% Clear Residuals 
clear all

%% Create Output Folder      


        %%%USER INPUT FOR SUBJECT NUMBER%%
    
        %Define the text file that PEBL and matlab will read

        dlg_title = 'Input';
        clear par
       
        while 1
            prompt = {sprintf('Step1. \nEnter path to save downsampled EEG data\n'),...
                      sprintf('\nStep 2.\nEnter path to get the EEG data to be downsampled\n'),...
                      sprintf('\nStep 3.(Option A or B)\nA)\nEnter dataset numbers  or use defaults (comma seperated)\n'),...
                      sprintf('B)\nUse default dataset numbers? \nDownsamples all files with 4 digits (like 0103) in "path to EEG data" folder \n\nOptions:\n1-Yes, downsample all files found in path (Step 1) \n0-No use files I provided above (Step 2)\n'),...
                      sprintf('\n\nStep 4.\nEnter downsample factor (integer, Original-fs/Desired-fs)')};
            dims = [1 60;1 60; 3 60; 1 60; 1 60];
            def = {'D:\GX Project\Results\DataDownSampled\',...
                   'D:\GX Project\Data\',...
                   '0101',...
                   '1',...
                   '2'};
               
            answer = inputdlg(prompt,dlg_title,dims,def);
            if ~isempty(answer)
                par.DatSaveLocIn = answer{1};
                par.DatLocIn = answer{2};
                
                par.DatasetsIncludedIn=strsplit(answer{3},',');
                par.UseDefaultDataIn=str2num(answer{4});
                par.DSfactorIn=str2num(answer{5});
                
                if ~exist(par.DatLocIn, 'dir')
%                     disp([{sprintf('ERROR in PATH TO DATA: \n\nThe path to the data enter does not exist.Please enter valid path.\nPath passed where data is was: %s', par.DatLocIn)}])
                    error_msg=sprintf(['ERROR in PATH TO DATA: \n\nSee dialog Step 2. \nThe path to the data enter does not exist.',...
                                      'Please enter valid path. \nAlso check that the path is formatted correctly and ends with a back slash.',...
                                      '\nPath you passed where data is was:\n %s'], par.DatLocIn);
                    error(error_msg)
%                     return
                    
                elseif ~exist(par.DatSaveLocIn,'dir')
%                     disp([{sprintf('ERROR in SAVE DATA PATH: \n\nThe path you entered to save data in does not exist.Please create is or enter valid path.\nPath passed to save data is: %s', par.DatSaveLocIn)}])
                    error_msg=sprintf(['ERROR in SAVE DATA PATH: \n\nSee dialog Step 1. \nThe path you entered to save data in does not exist.',...
                                       'Please create it or enter valid path. \nAlso check that the path is formatted correctly and ends with a back slash.',...
                                       '\nPath you passed to save data is:\n %s'], par.DatSaveLocIn);
                    error(error_msg)
%                     return
                else 
                    break
                end 
                

                break
            else 
                disp(['Input Canceled'])
                return
            end
        end




    
  
%% Define Data Locations and Files to Look At


% tic
% --- Start timing for overall script ---
t_script_start = tic;
% -------------------------------------


    %Where to save the data
    pathsave=par.DatSaveLocIn;%strcat('D:\GX Project\Results\DataDownSampled\');
    
    if ~exist(pathsave)
        mkdir(pathsave)
    end 

    %This is where all the data are stored
    DataLoc=par.DatLocIn;%'D:\GX Project\Data\'
    
    
    %These are the files that we want to look at
    if par.UseDefaultDataIn==1
       %These are the default files to load and downsample
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %                Pull Default Files
       %  Look through the folder and find all files that have numeric
       %   labels like 0101. 
       %       
               % Specify the folder path
               folderPath = DataLoc;

               % List all files in the folder
               fileList = dir(folderPath);

               % Initialize an empty cell array to store the filtered file names
               filteredFileNames = {};

               % Loop through the files
               for i = 1:length(fileList)
                   % Get the current file name
                   fileName = fileList(i).name;

                   % Check if the file name contains exactly four digits
                   if ~isempty(regexp(fileName, '\d{4}', 'once'))
                       % If it has four digits, add it to the filteredFileNames array
                       filteredFileNames{end+1} = fileName;
                   end
               end


                % Initialize an empty cell array to store the combined file names
                combinedFileNames = {};

                % Loop through the filtered file names and add them to the combinedFileNames array
                for i = 1:length(filteredFileNames)
                    combinedFileNames = [combinedFileNames, filteredFileNames{i}];
                end
               DatasetsIncluded=combinedFileNames; 

        %        DatasetsIncluded={'0101','0102','0103','0104','0201','0202'};

            else 
               DatasetsIncluded= par.DatasetsIncludedIn;
    end 
       %                Pull Default Files - End 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    
       
    %Strip off any spaces
    DatasetsIncluded=cellfun(@(x) regexprep(x, '\s+',''),DatasetsIncluded,'UniformOutput',false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Downsampling settings 
    DownSampleFactor=par.DSfactorIn; %Downsample by this factor
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for ii=1:length(DatasetsIncluded)    
    t_loop_start = tic;
    disp(['Running, downsampling dataset# ' num2str(ii),' of ',num2str(length(DatasetsIncluded))])  
    numcount=1;

    
%% _______________Import Performance Data__________________________________    

%Get the data path
disp(['Now running ' DatasetsIncluded{ii} '.....' ])
SelectedFle=strcat(DataLoc,DatasetsIncluded{ii},'\',DatasetsIncluded{ii},'\','ptracker-',DatasetsIncluded{ii},'.csv');

opts.SelectedVariableNames = [3,11]; 
opts.DataRange = '';

filename = SelectedFle;
delimiter = ',';
startRow = 2;
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

%We just want time and Performance 
ptrackerData{ii} = [dataArray{[3,11]}];%[dataArray{1:end-1}];
ptrackerData{ii}(:,1)=(ptrackerData{ii}(:,1)-ptrackerData{ii}(1,1))./1000; %Minus the 1st sample and convert to seconds

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

desiredFs = 100; 
ScreenFs = 60; 
ptrackerPerf{ii}=resample(ptrackerData{ii}(:,2),ptrackerData{ii}(:,1),desiredFs,desiredFs,ScreenFs);
ptrackerTime{ii}=[[0:length(ptrackerPerf{ii})-1]./desiredFs]'; 












%% _______________Import the EEG files ____________________________________

    %Define where each EEG file is
    GetFilesFrom=strcat(DataLoc,DatasetsIncluded{ii},'\');
    
    if ~exist( GetFilesFrom,'dir')
        numcount= numcount+1; %Added to keep the order of existing files
        return       
    end
    
    %Get the EEG file name to load
    Files=dir(fullfile(GetFilesFrom, '*.cnt'));
    
    %Define the filename of the EEG file to import
    filename= [char(Files(1).name)];
    
    
    EEG=[];%Create an empty EEG variable 
    %Create the full path to the EEG 
    PathData_EEprobe=[GetFilesFrom,filename];
    
    %Pull in just 5 samples to get information about the full EEG dataset. 
    Samp=read_eep_cnt(PathData_EEprobe,1,5);
    
    %Now pull in the whole EEG file. 
    EEG=read_eep_cnt(PathData_EEprobe,1,Samp.nsample);
    EEG.srate=EEG.rate;
    
    %There was an extra blank trigger of data for subject 1401. 
    %This removes that blank trigger to make concatination easy. 
    if strmatch(DatasetsIncluded{ii},'1401')
        EEG.triggers(1)=[];
    end 
    
    load_time = toc(t_loop_start);
    disp(['      EEG data ',DatasetsIncluded{ii},' loaded, ','time taken:', num2str(load_time), ' seconds'])


    
%% ______________Set Up Downsampling ______________________________________
    
    %Create empty variable to hold the downsampled data. 
    EEGDataDownSamp=[];
        
    fs=EEG.rate/DownSampleFactor; %Define the new sampling rate
    fsOld=EEG.rate;%Hold on to the old sampling rate
    DSpadding=ones(1,EEG.rate); %Add this much padding to signal to avoid edge artifacts when downsampling
    DSpadremove=length(DSpadding)/DownSampleFactor;%Remove these many samples after down sampling 
    Subj=DatasetsIncluded{ii};%Define the name of the file that we're downsampling 
    
    for nn=1:EEG.nchan%Loop through each channel in the EEG file
        DS=[]; %temp variable to hold downsampled output
        
        %Down sample operation
        %Padding added was the 1st and last sample of the EEG data.
        DS=resample([DSpadding*EEG.data(nn,1),EEG.data(nn,:),DSpadding*EEG.data(nn,end)],1,DownSampleFactor);
        EEGDataDownSamp(nn,:)=DS(DSpadremove+1:end-DSpadremove);%Remove the padded samples that were added to avoid edge artifacts 

      %Displays channel progression. 
%     disp(['Channel#', num2str(nn),' ', num2str(toc)])
    end 
    
%% __________ Convert triggers & Make Structure ___________________________
    
    %We will store the output as a structure
    clear DSamp mm Mont Montages

    %Make a timevector with the new sampling rate
    t2=[0:length(EEGDataDownSamp(1,:))-1]./fs;
    
    MatFiles=dir(fullfile(GetFilesFrom, '*.mat'));
    if ~isempty(MatFiles)
        load(strcat(GetFilesFrom,MatFiles.name),'Montages');
         if str2num(DatasetsIncluded{ii}(1:2))>10 %This checks if subject number is more than 10 (>10 is new protocol)
             %Check EEG for stim start triggers '16' and create an array
             %with all the stim type of that size. 
            Mont=repmat(Montages(1),sum(str2num(vertcat(EEG.triggers.code))==16),1)'; 
         else
%             Mont=repmat(Montages,sum(str2num(vertcat(EEG.triggers.code))==16),1);
            Mont=repmat(Montages,4,1); %NG-9/29/2023
            Mont=Mont(:)';
         end
        mm=1; %Counter
        AddStimType=1;%Flag to add stim type to Trigger struct
    end
    


    %Loop through each line and add the downsampled trigger codes/times to the
    %structure.
    for nn=1:length(vertcat(EEG.triggers.time))
        %Downsampled trigger time
     DSamp.triggers(nn).time   = t2(ceil([EEG.triggers(nn).offset]*(fs/fsOld)));
        %Downsampled trigger sample
     DSamp.triggers(nn).offset =ceil([EEG.triggers(nn).offset]*(fs/fsOld));  
        %Inherit the trigger code from EEG struct
     DSamp.triggers(nn).code   =EEG.triggers(nn).code;
        %Inherit the trigger type from EEG struct
     DSamp.triggers(nn).type   =EEG.triggers(nn).type;
     
             TrigType=EEG.triggers(nn).type;
             if  TrigType==2
             DSamp.triggers(nn).Label='Block Start';

             elseif TrigType==16
             DSamp.triggers(nn).Label='Stim Start';
             
             if AddStimType==1
                 if mm<=length(Mont)
                     DSamp.triggers(nn).StimType=upper(Mont{mm});
                     mm=mm+1;
                 else
                     DSamp.triggers(nn).StimType='Mistrigger';
                 end
             end 
             
             elseif TrigType==32
             DSamp.triggers(nn).Label='Stim Stop';

             end 
    end     
        
    %Write other usefull variables to the struct.
    DSamp.EEGdata=EEGDataDownSamp; %Downsamples EEG data
    DSamp.fs=fs;                %Downsampled fs
    DSamp.fsOld=fsOld;          %Data original fs
    DSamp.time=t2;              %Downsampled time vector
    DSamp.label=EEG.label;      %EEG channel labels
    DSamp.nchan=EEG.nchan;      %EEG channel numbers
    DSamp.rate=fs;              %Duplicate downsampled fs
    DSamp.npt=length(EEGDataDownSamp(1,:));%Length of downsampled data
    DSamp.Subj=DatasetsIncluded{ii};%Subject file number
    DSamp.ptrackerPerf=ptrackerPerf{ii}; %NG-05/26/2025 
    DSamp.ptrackerTime=ptrackerTime{ii}; %NG-05/26/2025 
    DSamp.ptrackerfs = desiredFs;
    
    
    %Save the output 
    save([pathsave,'EEG_DS_Struct_',DatasetsIncluded{ii}],'DSamp','-v7.3')
    load_time =  toc(t_loop_start) - load_time;
    disp(['      Downsampled EEG data for ',DatasetsIncluded{ii},' saved, ','time taken:', num2str(load_time) ' seconds'])


    
    %% Check Output
    % t2=[0:length(EEGDataDownSamp(1,:))-1]./fs;
    % figure; plot(EEG.time/1000,EEG.data(1,:)); hold on; plot(t2,EEGDataDownSamp(1,:))
    
load_time = toc(t_script_start);
fprintf(['      Elapsed time: %s seconds\n' ], num2str(round(load_time,3)))
end 





%% Display when done

elapsed_time = toc(t_script_start);
fprintf(['\n\n _______________________________________________________________________ \n' ])
fprintf(['\n\n      Time taken to extract all data: %s (H:M:S:Ms)      \n' ], datestr(elapsed_time/86400, 'HH:MM:SS:FFF'))
fprintf(['\n\n      DONE!      \n\n\n' ])
fprintf(['\n\n _______________________________________________________________________ \n' ])

%%








    
