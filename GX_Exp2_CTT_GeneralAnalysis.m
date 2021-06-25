%% GX_Exp2_CTT_GeneralAnalysis
% This script was written to examine the CTT data for Experiment 2 of the
% GX project. The primary goal is to extract all the trials and sort them
% by stimulation condition. Then get an aggregate measure of performance in
% relation to stimulation.
% 
% 
% Written by: Nigel Gebodh
% Date: January 2020
%
%
% Requirements:
% -Raincloudplots toolbox: 
% * https://wellcomeopenresearch.org/articles/4-63
% * https://peerj.com/preprints/27137v1.pdf
% -ANT Neuro file importer functions:
% * https://www.ant-neuro.com/support/supporting-documentation-and-downloads
%- Needs this file in same directory to pull montages GX_SubjectMontages_TaskDesign
%Internal:
%Some aspects of this code were taken from: GX_ZscoredPerfData.m
%


%% Clear Residuals 
clear all
close all

tic


%% Set Flags
SveAllpics=0; % Save the figure output? 0=No, 1=Yes
closefigs=1;  % Close all the figures periodocally? 0=No, 1=Yes
matlab_version='2019b';


%Double check versions. 
[ver]=version;
ver(end-14:end-10)
    if strmatch(ver(end-14:end-10),matlab_version)
        disp('Versions match moving on')
    else 
        error("SCRIPT ERROR: The MATLAB version you assigned in 'matlab_version' does not match the MATLAB verions detected")
    end 

%% Create Results folder        

    %This is where all the flagged figures will be saved. 
    %NOTE: All items in the folder will be deleted before saving new items!
%     pathsave=strcat('D:\GX Project\Results\DataOutput_Exp2_CTT\');
    pathsave=strcat('D:\GX\Results\DataOutput_Exp2_CTT\05292021\');
    
    prefix = strcat(pathsave);

    if SveAllpics==1 %1-Save output pics, 0-Don'd save output pics
        
        existance=exist(strcat(pathsave,'FigOutput'));
        if existance==0
            [s,m,mm]=mkdir(pathsave,'FigOutput');
            prefix = strcat(pathsave,'FigOutput','\');
        else
            delete([pathsave 'FigOutput\*.fig'])
            delete([pathsave 'FigOutput\*.png'])
            delete([pathsave 'FigOutput\*.pdf'])
            delete([pathsave 'FigOutput\*.eps']) 
            %rmdir([pathsave,'FigOutput'],'s'); %To erase the folder
            prefix = strcat(pathsave,'FigOutput','\');
        end
    end 

              
        
%% Define Data Locations and Files to Look At

%This is where all the data are stored
% DataLoc='D:\GX Project\Data\'
DataLoc='D:\GX\Data\'
%These are the files that we want to look at
DatasetsIncluded={'1101','1102',...
                  '1201','1202',...
                  '1301','1302',...
                  '1401','1402',...
                  '1501','1502',...
                  '1601','1602',...
                  '1801','1802',...
                  '1901','1902',...
                  '2001','2002',...
                  '2101','2102',...
                  '2201','2202',...
                  '2301','2302',...
                  '2401','2402',...
                  '2501','2502',...
                  '2601','2602'};


for ii=1:length(DatasetsIncluded)
SelectedFle=strcat(DataLoc,DatasetsIncluded{ii},'\',DatasetsIncluded{ii},'\','ptracker-',DatasetsIncluded{ii},'.csv');
filename = SelectedFle;



if strmatch(ver(end-14:end-10),'2018a')
%% Matlab version 2018a Import
    opts.SelectedVariableNames = [3,11]; 
    opts.DataRange = '';
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
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType',...
                'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    % Close the text file.
    fclose(fileID);
    clear opts

   
    %We just want time and Performance 
    ptrackerData{ii} = [dataArray{[3,11]}];%[dataArray{1:end-1}];
    ptrackerData{ii}(:,1)=(ptrackerData{ii}(:,1)-ptrackerData{ii}(1,1))./1000; %Minus the 1st sample and convert to seconds
    
elseif strmatch(ver(end-14:end-10),'2019b')    
%% Matlab version 2019b Import
    % Setup the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 13);

    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["subnum", "trial", "time", "posX", "posY", "userdeltaX", "userdeltaY", "timeDelta", "targetDeltaX", "targetDeltaY", "deviation", "mouseD1", "mouseD2"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    dataArray = readtable(filename, opts);

    % Convert to output type
    dataArray= table2array(dataArray);

    % Clear temporary variables
    clear opts

    
    
    %We just want time and Performance 
    ptrackerData{ii} = dataArray(:,[3,11,12,13, 4, 5]);%[dataArray{1:end-1}];
    ptrackerData{ii}(:,1)=(ptrackerData{ii}(:,1)-ptrackerData{ii}(1,1))./1000; %Minus the 1st sample and convert to seconds
else
    error("Data import not supported for this version of MATLAB please import data to dataArray manually and use 'Generate Script' option.")
end 
% 
% opts.SelectedVariableNames = [3,11]; 
% opts.DataRange = '';
% 
% 
% 
% delimiter = ',';
% startRow = 2;
% 
% % For more information, see the TEXTSCAN documentation.
% formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
% 
% % Open the text file.
% fileID = fopen(filename,'r');
% 
% % Read columns of data according to the format.
% % This call is based on the structure of the file used to generate this
% % code. If an error occurs for a different file, try regenerating the code
% % from the Import Tool.
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% 
% % Close the text file.
% fclose(fileID);
% 
% %We just want time and Performance 
% ptrackerData{ii} = [dataArray{[3,11,12,13, 4, 5]}];%[dataArray{1:end-1}]; 3-Time steps, Mouse Velocity-
% ptrackerData{ii}(:,1)=(ptrackerData{ii}(:,1)-ptrackerData{ii}(1,1))./1000; %Minus the 1st sample and convert to seconds




%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

desiredFs = 100; 
ScreenFs = 60; 
ptrackerPerf{ii}=resample(ptrackerData{ii}(:,2),ptrackerData{ii}(:,1),desiredFs,desiredFs,ScreenFs); % Radial deviation 
ptrackerVelocity{ii}=resample(ptrackerData{ii}(:,3),ptrackerData{ii}(:,1),desiredFs,desiredFs,ScreenFs); %Radial Velocity 
ptrackerAcceleration{ii}=resample(ptrackerData{ii}(:,4),ptrackerData{ii}(:,1),desiredFs,desiredFs,ScreenFs); %Radial  Acceleration
ptrackerXpos{ii}=resample(ptrackerData{ii}(:,5),ptrackerData{ii}(:,1),desiredFs,desiredFs,ScreenFs); %X position
ptrackerYpos{ii}=resample(ptrackerData{ii}(:,6),ptrackerData{ii}(:,1),desiredFs,desiredFs,ScreenFs); %Y position
ptrackerTime{ii}=[[0:length(ptrackerPerf{ii})-1]./desiredFs]'; 


clear ptrackerData



toc
 
 
 %% Getting the EEG 
Chans=[1:32];
numcount=ii;
        %Define where each EEG file is 
%         GetFilesFrom=strcat('D:\GX Project\Data\' ,DatasetsIncluded{ii},'\');
        GetFilesFrom=strcat(DataLoc ,DatasetsIncluded{ii},'\');
        if ~exist( GetFilesFrom)
            numcount= numcount+1; %Added to keep the order of existing files
            disp(['....Subject file not detected in folder: ' GetFilesFrom])
            disp(['..Skipping subject file: ' DatasetsIncluded{ii}])
            return
            
        end 
        %Get the EEG file name to load
         Files=dir(fullfile(GetFilesFrom, '*.cnt')); 



Files=dir(fullfile(GetFilesFrom, '*.cnt')); 
 
filename= [char(Files(1).name)];

EEG=[];
PathData_EEprobe=[GetFilesFrom,filename];    
   
Samp=read_eep_cnt(PathData_EEprobe,1,5); 

EEG=read_eep_cnt(PathData_EEprobe,1,Samp.nsample);
EEG.srate=2000;
EEG.nbchan=length(Chans);
EEG.etc=[];
EEG.trials=[];

%Remove extra trigger
if DatasetsIncluded{ii}=='1401', EEG.triggers(1)= []; end 

DataEEG{numcount}=EEG.data([Chans],:);
% DataEEG{numcount}=EEG.data([33,35],:);

AllEvents{numcount}=[EEG.triggers.offset];
AllEventsCode{numcount}={EEG.triggers.code};
AllEventsTime{numcount}=[EEG.triggers.time];

fs{numcount}=2000;   %EEG.rate;                 %Get the sampling rate
nSmp=[0:size(DataEEG{numcount},2)-1];%Created a vector the same size as the samples
t{numcount}=(nSmp)/fs{numcount};                %Created a time vector in sec
clear nSmp
N=size(DataEEG{numcount},2);

ref = [1:32];              %Electrodes to reference to
nchan=32;


clear DataEEG Adj_topoly_Each Samp baselineT 
BLcorDC{1,numcount}=[];
EEG.data=[];
EEG.time=[];
 
 
%% Set up events

    %Create a matrix of montages. 
    MatFiles=dir(fullfile(GetFilesFrom, '*.mat'));
    if ~isempty(MatFiles)
        load(strcat(GetFilesFrom,MatFiles.name),'Montages');
        MontHold=repmat(Montages,4,1);
        Mont=upper(MontHold(:)');

    end
     
 
 % _____________Looking At Each Stimulation Trial__________________________
 
 %Find all the Stim on Triggers
if  strcmp(DatasetsIncluded{ii},'1401')==1, AllEventsCode{ii}(1)=[]; end 

 Evnt_Stimstrt=AllEvents{ii}(find(str2num(vertcat(AllEventsCode{ii}{:}))==16));
 
 if DatasetsIncluded{ii}=='0102'
     Evnt_Stimstrt=Evnt_Stimstrt(1:end-1)
     Mont=Mont(1:length( Evnt_Stimstrt))
     MontAll(ii,1:length(Mont))=Mont;
 elseif DatasetsIncluded{ii}=='0101'
     Evnt_Stimstrt=Evnt_Stimstrt(1:end-1)
     Mont=Mont(1:length( Evnt_Stimstrt))
     MontAll(ii,1:length(Mont))=Mont;
 else 
      MontAll(ii,1:length(Mont))=Mont;
 end 
 
 clear Emp
 if sum(cellfun(@isempty,{MontAll{ii,:}}))>0
        Emp=find(cellfun(@isempty,{MontAll{ii,:}})); 
        for tt=1:length(Emp)
            MontAll{ii,Emp(tt)}='';
        end 
 end 
 
 
 
 
 %Pull the EEG events out so that we can sort and use them. 
 % We will use the EEG events to look a the behavioral data. 
 Evnt_Stimstrt2=(Evnt_Stimstrt-AllEvents{ii}(1)).*(desiredFs/fs{1});%
 Evnt_BlockStart=AllEvents{ii}(find(str2num(vertcat(AllEventsCode{ii}{:}))==2));
 Evnt_BlockStart2=(Evnt_Stimstrt-AllEvents{ii}(1)).*(desiredFs/fs{1});%
 startT=60*1.75;
 endT=60*2.5;
  
%Here we sort all the behavioral data into trials. 
 for mm=1:20;length(Evnt_Stimstrt);
 clear enUp enLw
 pta1=Evnt_Stimstrt(mm)-(startT*fs{1});
 pta2=((Evnt_Stimstrt2(mm))-(startT*desiredFs));
 
 ptb1=Evnt_Stimstrt(mm)+(endT*fs{1});
 ptb2=((Evnt_Stimstrt2(mm))+(endT*desiredFs));
 
 Tseg=-startT:1/fs{1}:endT;
 Tseg2=-startT:1/desiredFs:endT;
 
 

        PerfSorted{ii}(mm,:)=ptrackerPerf{ii}(pta2:ptb2,1);
        VelocitySorted{ii}(mm,:)=ptrackerVelocity{ii}(pta2:ptb2,1);
        AccelerationSorted{ii}(mm,:)=ptrackerAcceleration{ii}(pta2:ptb2,1);
        XPosSorted{ii}(mm,:)=ptrackerXpos{ii}(pta2:ptb2,1);
        YPosSorted{ii}(mm,:)=ptrackerYpos{ii}(pta2:ptb2,1);
 end     
        
   disp(['Done with file ' DatasetsIncluded{ii} ])     
 end 

toc

for ii=1:length(DatasetsIncluded)
 clear Emp
 if sum(cellfun(@isempty,{MontAll{ii,:}}))>0
        Emp=find(cellfun(@isempty,{MontAll{ii,:}})); 
        for tt=1:length(Emp)
            MontAll{ii,Emp(tt)}='';
        end 
 end 
end

%Change this later. 
%This just overwrites the montage map that we created earlier. 
clear MontAll
GX_SubjectMontages_TaskDesign2

% return 

%% Pull out performance
clr=[0 1 0; 1 0 0; 0 0 1];
AA=vertcat(DatasetsIncluded{:});
NumUniqueSubjs=length(unique(str2num(AA(:,1:2))));

clear  Varib CoeffVariation PercenMoreThanThres PerfPulledInMeanDevi PerfPulledInMeanDeviPerChange PerfPulledInMean
    DatInTlim={find( Tseg2==0)-(100*30):find( Tseg2==0); find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35); find( Tseg2==0)+(100*(40)):find( Tseg2==0)+(100*(40+30))};
    Ttrial=-30:1/100:60;
    DatInTlimTtrial={Ttrial(find(Ttrial==-30):find(Ttrial==0));Ttrial(find(Ttrial==0):find(Ttrial==30)); Ttrial(find(Ttrial==30):find(Ttrial==60))}
shifttleft=-25;
 for ii=1:length(DatasetsIncluded)
     for mm=1:size(PerfSorted{ii},1)
 figure; 
 subplot(1,3,1)
 
 for rr=1:3
     clear DatIn1
     meanDatIn=mean(PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0)));
     stdDatIn=std(PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0)));
     
     if rr==1, DatIn1=((PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0))) -meanDatIn)/stdDatIn;
     elseif rr==2, DatIn1=((PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)))-meanDatIn)/stdDatIn;
     elseif rr==3, DatIn1=((PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))))-meanDatIn)/stdDatIn;
     end
     

     
     
     DatMnMx(rr,:)=[min(DatIn1) ,max(DatIn1)];
     DatMnMx2(rr,:)=[min(stdDatIn) ,max(stdDatIn)];
      plot((DatInTlimTtrial{rr}),DatIn1,'Color',clr(rr,:),'Linewidth',2)
      hold on

        
 end 
             line([1 1]*0,[-1 1]*(max(DatMnMx(:,2))+max(DatMnMx2(:,2)))*2, 'Color','k', 'LineStyle','--','Linewidth',2)

             line([1 1]*30,[-1 1]*(max(DatMnMx(:,2))+max(DatMnMx2(:,2)))*2, 'Color','k', 'LineStyle','--','Linewidth',2)
             axis tight
             ylim([min(DatMnMx(:,1))-max(DatMnMx2(:,1)) max(DatMnMx(:,2))+max(DatMnMx2(:,2))])
             
             ylabel(['Z-Scored Deviation'])
             xlabel(['Time(sec)'])
             title(['Z-Scored Deviation' ])
             
             %Compute the Non-Zscored Integral
             
             PerfPulledInMeanDevi{ii,mm,1}= mean(PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0)));
             PerfPulledInMeanDevi{ii,mm,2}= mean(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)));
             PerfPulledInMeanDevi{ii,mm,3}= mean(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))));
             PerfPulledInMeanDeviPerChange{ii,mm}=100*((PerfPulledInMeanDevi{ii,mm,2}-PerfPulledInMeanDevi{ii,mm,1})/PerfPulledInMeanDevi{ii,mm,1});
             
             hold on
             
              txt = ['Percent Change \mu:' num2str(round(PerfPulledInMeanDeviPerChange{ii,mm})) '%'];
              text(shifttleft,max(DatMnMx(:,2))+max(DatMnMx2(:,2))-2,txt,'FontSize',12)


              txt = ['Before \mu:' num2str(round(PerfPulledInMeanDevi{ii,mm,1})), ''];
              text(shifttleft,max(DatMnMx(:,2))+max(DatMnMx2(:,2))-5,txt,'FontSize',12, 'Color',[0 1 0])

              txt = ['During \mu:' num2str(round(PerfPulledInMeanDevi{ii,mm,2})), ''];
              text(shifttleft,max(DatMnMx(:,2))+max(DatMnMx2(:,2))-7,txt,'FontSize',12, 'Color',[1 0 0])

              txt = ['After \mu:' num2str(round(PerfPulledInMeanDevi{ii,mm,3})), ''];
              text(shifttleft,max(DatMnMx(:,2))+max(DatMnMx2(:,2))-9,txt,'FontSize',12, 'Color',[0 0 1])
                
 subplot(1,3,2)
    h1 = raincloud_plot(((PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0))) -meanDatIn)/stdDatIn, 'box_on', 1, 'color', [0 1 0 ], 'alpha', 0.5,...
          'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,... %0.35
          'box_col_match', 1);
    line([1 1]*mean(((PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0))) -meanDatIn)/stdDatIn),[-1 1]*5, 'Color',[0 1 0 ], 'LineStyle','--')
      
    %During Stim 
      h2 = raincloud_plot(((PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)))-meanDatIn)/stdDatIn, 'box_on', 1, 'color', [1 0 0 ], 'alpha', 0.5,...
          'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
          'box_col_match', 1);
       line([1 1]*mean(((PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)))-meanDatIn)/stdDatIn),[-1 1]*5, 'Color',[1 0 0 ], 'LineStyle','--')
    %Post Stim 
      h3 = raincloud_plot(((PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))))-meanDatIn)/stdDatIn, 'box_on', 1, 'color', [0 0 1 ], 'alpha', 0.5,...
          'box_dodge', 1, 'box_dodge_amount', .95, 'dot_dodge_amount', 1.15,...
          'box_col_match', 1);
       line([1 1]*mean(((PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))))-meanDatIn)/stdDatIn),[-1 1]*5, 'Color',[0 0 1 ], 'LineStyle','--')

      LimtYMax= max([h1{1,1}.YData,h2{1,1}.YData,h3{1,1}.YData]);
      LimtYMax= LimtYMax+( LimtYMax*0.05);
      LimtYMin= min([h1{1,2}.YData,h2{1,2}.YData,h3{1,2}.YData]);
      LimtYMin=LimtYMin+(LimtYMin*0.05);

      set(gca,'YLim', [LimtYMin LimtYMax]);
      xlabel(['Z-Scored Deviation'])
      set(gca,'ytick',[])
      title(['Subj-' DatasetsIncluded{ii}, ' Trial-' num2str(mm),' Mont-' MontAll{ii,mm} ])
     
      
    subplot(1,3,3)  
    
clear DatInPlot1 DatInPlot2 DatInPlot3
     DatInPlot1=((PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0))) -meanDatIn)/stdDatIn;
     DatInPlot2=((PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)))-meanDatIn)/stdDatIn;
     DatInPlot3=((PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))))-meanDatIn)/stdDatIn;
    
     Zthresh=1.5;
     PercenMoreThanThres{ii,mm,1}=(sum(sort((DatInPlot1),'descend')>Zthresh)/length(DatInPlot1))*100;
     PercenMoreThanThres{ii,mm,2}=(sum(sort((DatInPlot2),'descend')>Zthresh)/length(DatInPlot2))*100;
     PercenMoreThanThres{ii,mm,3}=(sum(sort((DatInPlot3),'descend')>Zthresh)/length(DatInPlot3))*100;
     
     
     PerfPulledInZscoredMedian{ii,mm,1}=median(DatInPlot1);
     PerfPulledInZscoredMedian{ii,mm,2}=median(DatInPlot2);
     PerfPulledInZscoredMedian{ii,mm,3}=median(DatInPlot3);
     
     PerfPulledInZscoredMean{ii,mm,1}=mean(DatInPlot1);
     PerfPulledInZscoredMean{ii,mm,2}=mean(DatInPlot2);
     PerfPulledInZscoredMean{ii,mm,3}=mean(DatInPlot3);
     
     PerfPulledInMedian{ii,mm,1}= median(PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0)));
     PerfPulledInMedian{ii,mm,2}= median(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)));
     PerfPulledInMedian{ii,mm,3}= median(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))));
     
     %Arithmetic Mean
     PerfPulledInMean{ii,mm,1}= mean(PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0)));
     PerfPulledInMean{ii,mm,2}= mean(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)));
     PerfPulledInMean{ii,mm,3}= mean(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))));
  
     
      plot((sort((DatInPlot1),'descend')),(1/length(DatInPlot1):1/length(DatInPlot1):1)*100,'Color',[0 1 0 ],'Linewidth',2); hold on;
      plot((sort((DatInPlot2),'descend')),(1/length(DatInPlot1):1/length(DatInPlot1):1)*100,'Color',[1 0 0 ],'Linewidth',2); hold on
      plot((sort((DatInPlot3),'descend')),(1/length(DatInPlot1):1/length(DatInPlot1):1)*100,'Color',[0 0 1 ],'Linewidth',2);
      line([1 1]*Zthresh,[0 1]*100, 'Color',[0 0 1 ], 'Color',[0 0 0],'LineStyle','--')
      xlabel('z-score')
      ylabel('Percent (%)')
      ZscoreToDev=(Zthresh*stdDatIn)+meanDatIn;
      txt = ['Z-score ' num2str(Zthresh) ' = Devi ' num2str(round(ZscoreToDev))];
      text(0.5,80,txt,'FontSize',12)
      txt = (sum(sort((DatInPlot1),'descend')>Zthresh)/length(DatInPlot1))*100;
      
      txt = ['Before:' num2str(round(PercenMoreThanThres{ii,mm,1},2)), '%'];
      text(1.65,65,txt,'FontSize',12, 'Color',[0 1 0])
      
      txt = ['During:' num2str(round(PercenMoreThanThres{ii,mm,2},2)), '%'];
      text(1.65,60,txt,'FontSize',12, 'Color',[1 0 0])
      
      txt = ['After:' num2str(round(PercenMoreThanThres{ii,mm,3},2)), '%'];
      text(1.65,55,txt,'FontSize',12, 'Color',[0 0 1])
      axis tight
      title(['Percent Samples Above Threshold'])

      
    fname=[ 'Subj ZScored Perf-' DatasetsIncluded{ii} '-Trial'  num2str(mm),' Mont-' MontAll{ii,mm} ];
    set(gcf,'Name',fname,'Position', [376 400 1356 359],'PaperPositionMode','auto')
    
     if SveAllpics==1
           h = gcf;
           saveas(h,strcat(prefix,fname,'.fig'),'fig');
           saveas(h,strcat(prefix,fname,'.png'),'png');
           print(h,'-dpng', [prefix,fname], '-r600');
%            print(h,'-depsc', [prefix,fname], '-r600');
%            print(h,'-dpdf', [prefix,fname], '-r600');
        
     end
       close all
      if closefigs==1, close all,  end 

     end 
    
 end 

 
 
%Here we sort all the data that we gathered
 clear VarianceRatio CoeffVariationRatio PercenMoreThanThresPerDiff PrePostMeanDiffPerfPulledInMedian PrePostMeanDiffPerfPulledInMean PrePostMeanDiffPerfPulledInMeanPercent
 MontageMat2={'F30','M30',};
 
 for ii=1:length(DatasetsIncluded)
  for mm=1:size(PerfSorted{ii},1)
     if ~isempty(MontAll{ii,mm})
     IndMont=find(contains( MontageMat2,MontAll{ii,mm}));
     IndSubj=str2num(DatasetsIncluded{ii}(1:2));
     
     
     AA=contains(MontAll,MontAll{ii,mm});
     sum(AA(ii,:)==1);
     IndxCol=find(AA(ii,:));
     for IndTrial=1:sum(AA(ii,:)==1)
        
     %This is the X% increase or decrease (the percent diffence)
     PercenMoreThanThresPerDiff{IndSubj,IndMont,IndTrial}=((PercenMoreThanThres{ii,IndxCol(IndTrial),2}-PercenMoreThanThres{ii,IndxCol(IndTrial),1})/PercenMoreThanThres{ii,IndxCol(IndTrial),1})*100;
     
     PercenMoreThanThresPerDiffPrePostMean{IndSubj,IndMont,IndTrial}=((PercenMoreThanThres{ii,IndxCol(IndTrial),2}-mean([PercenMoreThanThres{ii,IndxCol(IndTrial),1},PercenMoreThanThres{ii,IndxCol(IndTrial),3}]))/...
                                                                     mean([PercenMoreThanThres{ii,IndxCol(IndTrial),1},PercenMoreThanThres{ii,IndxCol(IndTrial),3}]))*100;
     
     
      %This is the X% increase or decrease (the percent diffence) for
      %just the non-thresholded data.
     PercenDiffPerfPulledInMean{IndSubj,IndMont,IndTrial}=((PerfPulledInMean{ii,IndxCol(IndTrial),2}-PerfPulledInMean{ii,IndxCol(IndTrial),1})/PerfPulledInMean{ii,IndxCol(IndTrial),1})*100;
     GenDiffPerfPulledInMean{IndSubj,IndMont,IndTrial}=(PerfPulledInMean{ii,IndxCol(IndTrial),2}-PerfPulledInMean{ii,IndxCol(IndTrial),1});
     %Take the mean of the pre post and subtract it from during.
     PrePostMeanDiffPerfPulledInMean{IndSubj,IndMont,IndTrial}=(PerfPulledInMean{ii,IndxCol(IndTrial),2}-mean([PerfPulledInMean{ii,IndxCol(IndTrial),1},PerfPulledInMean{ii,IndxCol(IndTrial),3}]));
     
     PrePostMeanDiffPerfPulledInMeanPercent{IndSubj,IndMont,IndTrial}=(PerfPulledInMean{ii,IndxCol(IndTrial),2}-mean([PerfPulledInMean{ii,IndxCol(IndTrial),1},PerfPulledInMean{ii,IndxCol(IndTrial),3}]))./...
                                                               mean([PerfPulledInMean{ii,IndxCol(IndTrial),1},PerfPulledInMean{ii,IndxCol(IndTrial),3}]);
                                                           
     PercenDiffPerfPulledInMedian{IndSubj,IndMont,IndTrial}=((PerfPulledInMedian{ii,IndxCol(IndTrial),2}-PerfPulledInMedian{ii,IndxCol(IndTrial),1})/PerfPulledInMedian{ii,IndxCol(IndTrial),1})*100;
     GenDiffPerfPulledInMedian{IndSubj,IndMont,IndTrial}=(PerfPulledInMedian{ii,IndxCol(IndTrial),2}-PerfPulledInMedian{ii,IndxCol(IndTrial),1});
     %Take the mean of the pre post and subtract it from during.
     PrePostMeanDiffPerfPulledInMedian{IndSubj,IndMont,IndTrial}=(PerfPulledInMedian{ii,IndxCol(IndTrial),2}-mean([PerfPulledInMedian{ii,IndxCol(IndTrial),1},PerfPulledInMedian{ii,IndxCol(IndTrial),3}]));
     

     end 
     end 
  end 
 end
    
 
 %%% 
 clear PercenMoreThanThresPerDiffMean  CoeffVariationRatioMean  VarianceRatioMean PercenMoreThanThresPerDiffMedian
 clear GenDiffPerfPulledInMeanPooled GenDiffPerfPulledInMedianPooled PerDiffPerfPulledInMeanPooled PerDiffPerfPulledInMedianPooled
 clear PerDiffPerfPulledInMeanZPooled PerDiffPerfPulledInMedianZPooled
 

clear GenDiffAccPulledInMeanPooled GenDiffAccPulledInMedianPooled PerDiffAccPulledInMeanPooled PerDiffAccPulledInMedianPooled
clear GenDiffVelPulledInMeanPooled GenDiffVelPulledInMedianPooled PerDiffVelPulledInMeanPooled PerDiffVelPulledInMedianPooled
 
AA=vertcat(DatasetsIncluded{:});
NumUniqueSubjsNums=(unique(str2num(AA(:,1:2))));

  for ii=1:size(PercenMoreThanThresPerDiff,1)-10;%NumUniqueSubjs

      iii=ii+10;
  if ~isempty(vertcat(PercenMoreThanThresPerDiff{iii,1,:}))
  for mm=1:length(MontageMat2)
     if ~isempty(MontAll{ii,mm})
     IndMont=find(contains( MontageMat2,MontAll{ii,mm}));
     IndSubj=str2num(DatasetsIncluded{ii}(1:2));
     
     
     AA=contains(MontAll,MontAll{ii,mm});
     sum(AA(ii,:)==1);
     IndxCol=find(AA(ii,:));

        
     %This is the X% increase or decrease (the percent diffence)
     PercenMoreThanThresPerDiffMean{ii,mm}=nanmean([PercenMoreThanThresPerDiff{iii,mm,:}]);
     PercenMoreThanThresPerDiffMedian{ii,mm}=nanmedian([PercenMoreThanThresPerDiff{iii,mm,:}]);
     
     
     %General differnece in performance 
     PerDiffPerfPulledInMeanPooled{ii,mm}=nanmean([PercenDiffPerfPulledInMean{iii,mm,:}]);
     PerDiffPerfPulledInMedianPooled{ii,mm}=nanmedian([PercenDiffPerfPulledInMedian{iii,mm,:}]);
     
     GenDiffPerfPulledInMeanPooled{ii,mm} =nanmean([GenDiffPerfPulledInMean{iii,mm,:}]);
     GenDiffPerfPulledInMedianPooled{ii,mm} =nanmedian([GenDiffPerfPulledInMedian{iii,mm,:}]);
     
     PrePostMeanDiffPerfPulledInMeanPooled{ii,mm}=nanmean([PrePostMeanDiffPerfPulledInMean{iii,mm,:}]);
     PrePostMeanDiffPerfPulledInMedianPooled{ii,mm}=nanmedian([PrePostMeanDiffPerfPulledInMedian{iii,mm,:}]);
     

     end 
  end 
  end 
  end
 
%  return 
 
 
 
 %%
AA=vertcat(DatasetsIncluded{:});
NumUniqueSubjsNums=(unique(str2num(AA(:,1:2))));
clear DatOut 
PlotThese={'PercenMoreThanThresPerDiff{ii,nn,:}',1,[0,100],'Trials With Improvments More than Thres','Trials w/','Mean Deviation Improvement During Stim PrePost MeanCorrection-All'};
% DatOut=ones(2,2,2);

for kk=1:2%6
    figure; 
    clear cc;cc=1;
for ii=min(NumUniqueSubjsNums):max(NumUniqueSubjsNums);
    
     for nn=1:2;
         
         clear FF; 
         if kk==1
         FF=(vertcat(PercenMoreThanThresPerDiffPrePostMean{ii,nn,:}));%This is sorted!
         fname=['All trials and subjects-PercentMorethanThreshold'];
         MultOne=-1;
         elseif kk==2
             
         %PrePostMeanDiffPerfPulledInMeanPercent
         FF=100*(vertcat(PrePostMeanDiffPerfPulledInMeanPercent{ii,nn,:}));  %This is sorted!

         fname=['All trials and subjects-PercenDiffinMeanDeviation-Perf'];
         MultOne=-1;
         elseif kk==3
         FF=(vertcat(PrePostMeanDiffVelPulledInMean{ii,nn,:}));  
         fname=['All trials and subjects-PercenDiffinMeanDeviation-Velocity PrePost Corr'];
         MultOne=1;
         elseif kk==4
         FF=(vertcat(PrePostMeanDiffAccPulledInMean{ii,nn,:}));  
         fname=['All trials and subjects-PercenDiffinMeanDeviation-Acceleration PrePost Corr']; 
         MultOne=1;
         elseif kk==5
         FF=(vertcat(PercenDiffVelPulledInMean{ii,nn,:}));  
         fname=['All trials and subjects-PercenDiffinMeanDeviation-Velocity'];
         MultOne=1;
         elseif kk==6
         FF=(vertcat(PercenDiffAccPulledInMean{ii,nn,:}));  
         fname=['All trials and subjects-PercenDiffinMeanDeviation-Acceleration']; 
         MultOne=1;
         end
         
    if ~isempty(FF)
        subplot(size(DatasetsIncluded,2)/2,2,cc),

         imagesc(MultOne*FF'), 
         maxFF=max(abs(FF));
         if kk<2 && maxFF>100, maxFF=100; end 
         caxis([-1*maxFF maxFF]); 
         set(gca,'ytick',[])
         colormap(flipud(redblue))
         hb=colorbar('location','eastoutside')

         ylabel(hb, '%');
         if nn==2, ThisMont='M30'; else;  ThisMont='F30'; end 
         ylabel([{['Subj-' num2str(DatasetsIncluded{cc}(1:2))]},{['-' ThisMont]}])
         if cc>size(DatasetsIncluded,2)-2
             xlabel('Trials')
         end 
         if contains(MontAll{cc},'F30'),jj=1; else, jj=2; end 
         DatOut(ii-10,jj,kk)=(sum((FF)<0)/20)*100;
           cc=cc+1;
     end 
       
         
     end, 
end

set(gcf,'Name',fname, 'Position',[ 898         215         936        1043])

         if SveAllpics==1
               h = gcf;
               saveas(h,strcat(prefix,fname,'.fig'),'fig');
               saveas(h,strcat(prefix,fname,'.png'),'png');
               print(h,'-dpng', [prefix,fname], '-r600');

           end
           if closefigs==1, close all,  end 

PlotThese={'_',1,[0,100],'Trials With Improvments (%)','% Trials w/ Reduced Extreme Events','Percent of trials with reduced extreme events';...
           '_',1,[0,100],'Trials With Improvments (%)','% Trials w/ Reduced Mean Deviation','Percent of trials with reduced mean deviation';...
           '_',1,[0,100],'Trials With Improvments (%)','% Trials w/ Reduced Extreme Events','Percent of trials with reduced extreme events';...
           '_',1,[0,100],'Trials With Improvments (%)','% Trials w/ Reduced Mean Deviation','Percent of trials with reduced mean deviation';
           '_',1,[0,100],'Trials With Improvments (%)','% Trials w/ Reduced Extreme Events','Percent of trials with reduced extreme events';...
           '_',1,[0,100],'Trials With Improvments (%)','% Trials w/ Reduced Mean Deviation','Percent of trials with reduced mean deviation'};
rr=kk;
clear DatOut2
if kk==1,
    DatOut(find(DatOut(:,1)==0),:)=[];
else
    ind=find(DatOut(:,1,kk)==0);

    DatOut(ind,:,:)=[];
end

DatOut2=DatOut;
DatIn=(DatOut2(:,:,kk)).*PlotThese{1,2};

 clear x y tttext
    figure;
    y = repmat(1:size(DatIn,1),size(DatIn,2),1)'; % generate y-coordinates
    x =repmat([1:2],length(NumUniqueSubjsNums),1); % generate x-coordinates
    tttext=num2cell((( round(DatIn,2))));
    tttext = cellfun(@num2str, tttext, 'UniformOutput', false); % convert to string
    imagesc(1*((DatIn))); 
    text(x(:), y(:), tttext, 'HorizontalAlignment', 'Center','fontsize',14,'Color',[249 166 2]./255)
    caxis([PlotThese{rr,3}]); 
    ylabel('Subjects')
    xlabel('Stimulation Conditions')
    set(gca,'XTick',[1:length(MontageMat2)],'XTickLabels', MontageMat2)
    set(gca,'YTick',[1:length(NumUniqueSubjsNums)],'YTickLabels', num2str(NumUniqueSubjsNums))
    set(gca,'Fontsize',16)
    hold on
    plot(repmat(3.5,1,12), 0:11,'Color','k','Linewidth',3)
    plot(repmat(6.5,1,12), 0:11,'Color','k','Linewidth',3)
   
    hb=colorbar;
    ylabel(hb, PlotThese{rr,4},'Fontsize',16);
    colormap(flipud(bone))
    title([PlotThese{rr,5}])

    fname=[ PlotThese{rr,6}];
    set(gcf,'Name',fname) 
    
    
         if SveAllpics==1
               h = gcf;
               saveas(h,strcat(prefix,fname,'.fig'),'fig');
%                saveas(h,strcat(prefix,fname,'.png'),'png');
               print(h,'-dpng', [prefix,fname], '-r600');
               print(h,'-depsc', [prefix,fname], '-r600');

           end
           if closefigs==1, close all,  end 
    
end 
 disp('done')
 

 %% Behavior Before During After Stimulation
AAnames=vertcat(DatasetsIncluded{:});
NumUniqueSubjsNums=(unique(str2num(AAnames(:,1:2))));

 clrs = cbrewer('seq', 'Blues', 20, 'pchip');
 bxpltclr =[241, 161, 4]./255;%[255,20,147]./255; [1, 146, 255]./255;%
 bxpltclr2 =[255, 233, 121]./255; %[0, 255, 249]./255;%
 bxpltmedclr =[165, 129, 5]./255;
 Meanlinesclr=[255, 96, 0]./255;%[0, 116, 63]./255; 


for kk=1
for ii=1:length(AAnames);clear AA BB CC; 
    if kk==1
    AA=([PerfSorted{1,ii}(:,find( Tseg2==0)-(100*30):find( Tseg2==0))]);
    BB=([PerfSorted{1,ii}(:,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35))]);
    CC=(PerfSorted{1,ii}(:,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))));
    dattype='Deviation';
    
    elseif kk==2
    AA=([VelocitySorted{1,ii}(:,find( Tseg2==0)-(100*30):find( Tseg2==0))]);
    BB=([VelocitySorted{1,ii}(:,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35))]);
    CC=(VelocitySorted{1,ii}(:,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))));
    dattype='Vel';
    elseif kk==3
    AA=([AccelerationSorted{1,ii}(:,find( Tseg2==0)-(100*30):find( Tseg2==0))]);
    BB=([AccelerationSorted{1,ii}(:,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35))]);
    CC=(AccelerationSorted{1,ii}(:,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))));
    dattype='Acc';
    end
    
    
    if mod(ii,2)>0,  figure;    if strcmp(MontAll{ii},'F30')==1; gg=1; else gg=2; end; subplot(1,2,gg), else, if strcmp(MontAll{ii},'F30')==1; gg=1; else gg=2; end; subplot(1,2,gg), end 

    for ss=1:20
      plot(1:3,([mean(AA(ss,:),2),mean(BB(ss,:),2),mean(CC(ss,:),2)])','-','linewidth',2, 'color', clrs(ss,:)); hold on
      colormap(clrs);
%       caxis([0.5, 20]); 
%       colorbar
    end 

      line([1 1]*1,mean([mean(AA,2)])'+([-1 1].*(std([mean(AA,2)]')./sqrt(20))),'Color',Meanlinesclr,'linewidth',2); hold on
      line([1 1]*2,mean([mean(BB,2)])'+([-1 1].*(std([mean(BB,2)]')./sqrt(20))),'Color',Meanlinesclr,'linewidth',2); hold on
      line([1 1]*3,mean([mean(CC,2)])'+([-1 1].*(std([mean(CC,2)]')./sqrt(20))),'Color',Meanlinesclr,'linewidth',2); hold on
    hold on
     
    plot(1:3,mean([mean(AA,2),mean(BB,2),mean(CC,2)])','*-','linewidth',3, 'color',Meanlinesclr)
    
    hold on
    
    hh=boxplot([mean(AA,2),mean(BB,2),mean(CC,2)],'notch','off','BoxStyle','outline','Widths',0.2,'symbol',''); hold on;

    set(findobj(hh,'type','line','Tag','Upper Whisker','Tag','Lower Whisker'),'LineStyle','-');
    clear jj; jj=findobj(hh,'type','line');
    set(jj([1 2 8 9 15 16]),'LineStyle','-');
    hh_out = findobj(hh,'Tag','Box'); 
        for j=1:length(hh_out) 
        patch(get(hh_out(j),'XData'),get(hh_out(j),'YData'),bxpltclr2,'EdgeColor',bxpltclr2,'FaceAlpha',.5 ,'LineStyle','-');
        end 
    set(hh,{'linew'},{3},{'color'},{bxpltclr})

    title(['Subj:' DatasetsIncluded{ii}(1:2)  '-' MontAll{ii} '-' dattype ' Across Trials'])
    ylim([10 max(max(([mean(AA,2),mean(BB,2),mean(CC,2)])'))+5])
    xlim([0.75 3.25])
    ylabel([ dattype])
    xlabel('Time Period')
    set(gca,'XTick',[1:3],'XTickLabels', {'Pre', 'During','Post'})
    
    if mod(ii,2)==0
        
    fname=['Subj-' DatasetsIncluded{ii}(1:2) '-Both Montages-'  dattype ' Across Trials-Boxplot and Line'];
    set(gcf,'Name',fname) 

         if SveAllpics==1
               h = gcf;
               saveas(h,strcat(prefix,fname,'.fig'),'fig');
%                saveas(h,strcat(prefix,fname,'.png'),'png');
               print(h,'-dpng', [prefix,fname], '-r600');
               print(h,'-deps', [prefix,fname], '-r600');
               print(h,'-dpdf', [prefix,fname], '-r600');

           end
           if closefigs==1, close all,  end 
           
     end 
    
end 
end 



%%
  
 
AA=vertcat(DatasetsIncluded{:});
NumUniqueSubjsNums=(unique(str2num(AA(:,1:2))));


 
 
 %Performance of the mean and medians of Z-scored and  non-Z-scored
  
  %Non-Zscored- Mean
 %%% Looking at All mean deviations from the mean of pre post
 %Essentially mean(Pre_i ,Post_i)-During_i ; i=1:4
 clear  PlotThese
 PlotThese={'PercenMoreThanThresPerDiffMean'  ,-1,[-1,1].*max(max(abs(cell2mat(eval('PercenMoreThanThresPerDiffMean'))))),'Improvement In Extreme Events(%)','Mean Improvment in Extreme Events','Mean Improvments in Extreme Events-All';...
            'PercenMoreThanThresPerDiffMedian',-1,[-1,1].*max(max(abs(cell2mat(eval('PercenMoreThanThresPerDiffMedian'))))),'Improvement In Extreme Events(%)','Median Improvment in Extreme Events','Median Improvments in Extreme Events-All';...
            'PerDiffPerfPulledInMeanPooled'   ,-1,[-1,1].*max(max(abs(cell2mat(eval('PerDiffPerfPulledInMeanPooled' ))))),'Deviation Change(%)','Mean Deviation During Stim','Mean Deviation Improvement During Stim PrePost MeanCorrection-All';...
            'PerDiffPerfPulledInMedianPooled' ,-1,[-1,1].*max(max(abs(cell2mat(eval('PerDiffPerfPulledInMedianPooled'))))),'Deviation Change(%)','Median Deviation During Stim','Median Deviation Improvement During Stim PrePost MeanCorrection-All'};
% PercenMoreThanThresPerDiffMean
 
 for rr=1:size(PlotThese,1)
 DatIn=cell2mat(eval(PlotThese{rr,1}))*PlotThese{1,2};
%  ind=find(DatOut(:,1,kk)==0);
%  DatOut(ind,:,:)=[];
 
 clear x y tttext
    figure;
    y = repmat(1:size(DatIn,1),size(DatIn,2),1)'; % generate y-coordinates
    x =repmat([1:2],length(NumUniqueSubjsNums),1); % generate x-coordinates
    tttext=num2cell((round( DatIn,2)));
    tttext = cellfun(@num2str, tttext, 'UniformOutput', false); % convert to string
    imagesc(1*((DatIn))); 
    text(x(:), y(:), tttext, 'HorizontalAlignment', 'Center','fontsize',14,'Color',[249 166 2]./255)
    caxis([PlotThese{rr,3}]); 
    ylabel('Subjects')
    xlabel('Stimulation Conditions')
    set(gca,'XTick',[1:length(MontageMat2)],'XTickLabels', MontageMat2)
    set(gca,'YTick',[1:length(NumUniqueSubjsNums)],'YTickLabels', num2str(NumUniqueSubjsNums))
    set(gca,'Fontsize',16)
    hold on
    plot(repmat(3.5,1,12), 0:11,'Color','k','Linewidth',3)
    plot(repmat(6.5,1,12), 0:11,'Color','k','Linewidth',3)

    % scatter(MedianPerfTrial(:,[1,4,7,2,5,8,3,6,9],1)<15)
    
    hb=colorbar;
    ylabel(hb, PlotThese{rr,4},'Fontsize',16);
    colormap(flipud(redblue))
    title([PlotThese{rr,5}])
    fname=[ PlotThese{rr,6}];
    set(gcf,'Name',fname) 

         if SveAllpics==1
               h = gcf;
               saveas(h,strcat(prefix,fname,'.fig'),'fig');
               saveas(h,strcat(prefix,fname,'.png'),'png');
               print(h,'-dpng', [prefix,fname], '-r600');
               print(h,'-depsc', [prefix,fname], '-r600');
    %            print(h,'-dpdf', [prefix,fname], '-r600');
           end
           if closefigs==1, close all,  end 
 end 
 
 
%Plotting the means of F30 and M30


figure
 for rr=3
 DatIn=cell2mat(eval(PlotThese{rr,1}))*PlotThese{1,2};
 DatIn_SE=  nanstd(DatIn)./sqrt(sum(~isnan(DatIn),1));
 errorbar([1:2],mean(DatIn),DatIn_SE*1.96)
  xlim([0.5 2.5])
 set(gca,'XTick',[1:2],'XTickLabels', MontageMat2)
ylabel('Change(%)')
 
 end 
 
 
 
 disp('done')
 
 
 figure;
 errorbar([1:2],mean(cell2mat(Testt)),   std(cell2mat(Testt))./    sqrt(sum(~isnan(cell2mat(Testt)),1)) *1.96 )
   xlim([0.5 2.5])
 set(gca,'XTick',[1:2],'XTickLabels', MontageMat2)
ylabel('Change(%)')
 
Testt2= cell2mat(Testt);
%Make a blue colorbar to input to figure;  
 figure;
%  plot(1,1,'color', clrs(1,:));
hb=colorbar;
ylabel(hb, PlotThese{rr,4},'Fontsize',16);
colormap(flipud( cbrewer('seq', 'Blues', 20, 'pchip')))
fname='Dummy_fig_bluecolorbar';
set(gcf,'Name',fname) 
if SveAllpics==1
   h = gcf;
   saveas(h,strcat(prefix,fname,'.fig'),'fig');
%    saveas(h,strcat(prefix,fname,'.png'),'png');
   print(h,'-dpng', [prefix,fname], '-r600');
   print(h,'-depsc', [prefix,fname], '-r600');
   print(h,'-dpdf', [prefix,fname], '-r600');
end
if closefigs==1, close all,  end 
 
 
 
 
 
