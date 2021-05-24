%% GX_Exp1_CTT_GeneralAnalysis
% This script was written to examine the CTT data for Experiment 1 of the
% GX project. The primary goal is to extract all the trials and sort them
% by stimulation condition. Then get an aggregate measure of performance in
% relation to stimulation.
% 
% 
% Written by: Nigel Gebodh
% Date: October 2019
%
%
% Requirements:
% -Raincloudplots toolbox: 
% * https://wellcomeopenresearch.org/articles/4-63
% * https://peerj.com/preprints/27137v1.pdf
% -ANT Neuro file importer functions:
% * https://www.ant-neuro.com/support/supporting-documentation-and-downloads
%
%Internal:
%Some aspects of this code were taken from: GX_ZscoredPerfData.m
%
%% Clear Residuals 
clear all
close all

tic


%% Set Flags
%Set if to save the figure (1-Yes save, 0-No don't save) or to close all
%plotted figures (1-Yes plot, 0-No don't plot). 
SveAllpics=0; 
closefigs=0;
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


%     pathsave=strcat('D:\GX Project\Results\DataOutput_Exp1_CTT\');
    pathsave=strcat('D:\GX\Results\DataOutput_Exp1_CTT\');
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
            %             rmdir([pathsave,'FigOutput'],'s'); %To erase the folder
            prefix = strcat(pathsave,'FigOutput','\');
        end
    end 

              
        
%% Define Data Locations and Files to Look At

%This is where all the data are stored
% DataLoc='D:\GX Project\Data\'
DataLoc='D:\GX\Data\';
%These are the files that we want to look at
DatasetsIncluded={'0102','0103','0104',...
                  '0201','0202',...
                  '0301','0302','0303',...
                  '0401','0402','0403',...
                  '0501','0504','0505',...
                  '0601','0602','0603',...
                  '0701','0702','0703',...
                  '0801','0802','0803',...
                  '0901','0902','0903',...
                  '1001','1002','1003'};


              
Dat05mA=DatasetsIncluded;
Dat20mA={};
enUpAll=[]; enLwAll=[];


for ii=1:length(DatasetsIncluded)
SelectedFle=strcat(DataLoc,DatasetsIncluded{ii},'\',DatasetsIncluded{ii},'\','ptracker-',DatasetsIncluded{ii},'.csv');
filename = SelectedFle;%'D:\GX Project\Data\0101\0101\ptracker-0101.csv';

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
    ptrackerData{ii} = dataArray(:,[3,11]);%[dataArray{1:end-1}];
    ptrackerData{ii}(:,1)=(ptrackerData{ii}(:,1)-ptrackerData{ii}(1,1))./1000; %Minus the 1st sample and convert to seconds
else
    error("Data import not supported for this version of MATLAB please import data to dataArray manually and use 'Generate Script' option.")
end 









%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

desiredFs = 100; %Desired upsampled sampling frequency for CTT
ScreenFs = 60; 
ptrackerPerf{ii}=resample(ptrackerData{ii}(:,2),ptrackerData{ii}(:,1),desiredFs,desiredFs,ScreenFs);
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


% Some info below is redundant because it's based off some generic import code 
Files=dir(fullfile(GetFilesFrom, '*.cnt')); 
filename= [char(Files(1).name)];

EEG=[];
PathData_EEprobe=[GetFilesFrom,filename];    
   
%We're just going to pull in 5 samples because we just want the triggers
%because we're not looking at the EEG data yet and all the triggers are
%stored in the EEG files. 
Samp=read_eep_cnt(PathData_EEprobe,1,5); %Pull out 5 samples
EEG=read_eep_cnt(PathData_EEprobe,1,Samp.nsample);
EEG.srate=2000;
EEG.nbchan=length(Chans);
EEG.etc=[];
EEG.trials=[];


DataEEG{numcount}=EEG.data([Chans],:);

%Get all the trigger events
AllEvents{numcount}=[EEG.triggers.offset];
AllEventsCode{numcount}={EEG.triggers.code};
AllEventsTime{numcount}=[EEG.triggers.time];

fs{numcount}=2000;                   %Get the sampling rate
nSmp=[0:size(DataEEG{numcount},2)-1];%Created a vector the same size as the samples
t{numcount}=(nSmp)/fs{numcount};     %Created a time vector in sec
clear nSmp

clear DataEEG Adj_topoly_Each Samp baselineT 
BLcorDC{1,numcount}=[];
EEG.data=[];
EEG.time=[];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %Save each subjects EEG and performance data
% % BLcorDC_out=BLcorDC{ii};
% % ptrackerPerf_out=ptrackerPerf{ii};
% % ptrackerTime_out=ptrackerTime{ii};
% % AllTrigs=EEG.triggers;
% % save(strcat(prefix,DatasetsIncluded{ii},'_EEG.mat'),'BLcorDC_out', 'ptrackerPerf_out','ptrackerTime_out','AllTrigs')
% % clear ptrackerPerf_out ptrackerTime_out BLcorDC_out AllTrigs
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
 
 
 
%% Sort out the montages

    MatFiles=dir(fullfile(GetFilesFrom, '*.mat'));
    if ~isempty(MatFiles)
        load(strcat(GetFilesFrom,MatFiles.name),'Montages');
        MontHold=repmat(Montages,4,1);
        Mont=upper(MontHold(:)');
%         MontAll{ii}= Mont; 
    end
     
 
 % _____________Looking At Each Stimulation Trial__________________________
 
 %Find all the Stim on Triggers (code:16)
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
 
 
 
 
 
 Evnt_Stimstrt2=(Evnt_Stimstrt-AllEvents{ii}(1)).*(desiredFs/fs{1});%
 Evnt_BlockStart=AllEvents{ii}(find(str2num(vertcat(AllEventsCode{ii}{:}))==2));
 Evnt_BlockStart2=(Evnt_Stimstrt-AllEvents{ii}(1)).*(desiredFs/fs{1});%
 startT=60*1.75;
 endT=60*2.5;
 
 

%  NumStimTrials=length(Evnt_Stimstrt);
 for mm=1:length(Evnt_Stimstrt)
 clear enUp enLw
 pta1=Evnt_Stimstrt(mm)-(startT*fs{1});
 pta2=(Evnt_Stimstrt2(mm))-(startT*desiredFs);
 
 ptb1=Evnt_Stimstrt(mm)+(endT*fs{1});
 ptb2=(Evnt_Stimstrt2(mm))+(endT*desiredFs);
 
 Tseg=-startT:1/fs{1}:endT;
 Tseg2=-startT:1/desiredFs:endT;
 
        PerfSorted{ii}(mm,:)=ptrackerPerf{ii}(pta2:ptb2,1);

 end     
        
   disp(['Done with file ' DatasetsIncluded{ii} ])     
 end 


for ii=1:length(DatasetsIncluded)
 clear Emp
 if sum(cellfun(@isempty,{MontAll{ii,:}}))>0
        Emp=find(cellfun(@isempty,{MontAll{ii,:}})); 
        for tt=1:length(Emp)
            MontAll{ii,Emp(tt)}='';
        end 
 end 
end 









%% Look through each trial

clr=[0 1 0; 1 0 0; 0 0 1];
% close all

clear  Varib CoeffVariation PercenMoreThanThres
    DatInTlim={find( Tseg2==0)-(100*30):find( Tseg2==0); find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35); find( Tseg2==0)+(100*(40)):find( Tseg2==0)+(100*(40+30))};
    Ttrial=-30:1/100:60;
    DatInTlimTtrial={Ttrial(find(Ttrial==-30):find(Ttrial==0));Ttrial(find(Ttrial==0):find(Ttrial==30)); Ttrial(find(Ttrial==30):find(Ttrial==60))}

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
%              hold on
             line([1 1]*30,[-1 1]*(max(DatMnMx(:,2))+max(DatMnMx2(:,2)))*2, 'Color','k', 'LineStyle','--','Linewidth',2)
             axis tight
             ylim([min(DatMnMx(:,1))-max(DatMnMx2(:,1)) max(DatMnMx(:,2))+max(DatMnMx2(:,2))])
             
             ylabel(['Z-Scored Deviation'])
             xlabel(['Time(sec)'])
             title(['Z-Scored Deviation' ])
    
    
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
    Varib{ii,mm,1}=var(PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0)));  
    Varib{ii,mm,2}=var(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)));    
    Varib{ii,mm,3}=var(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30)))); 
    
    CoeffVariation{ii,mm,1}=std(PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0)))/mean(PerfSorted{1,ii}(mm,find( Tseg2==0)-(100*30):find( Tseg2==0)));  
    CoeffVariation{ii,mm,2}=std(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)))/mean(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*5):find( Tseg2==0)+(100*35)));    
    CoeffVariation{ii,mm,3}=std(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30))))/mean(PerfSorted{1,ii}(mm,find( Tseg2==0)+(100*40):find( Tseg2==0)+(100*(40+30)))); 
     

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
 clear VarianceRatio CoeffVariationRatio PercenMoreThanThresPerDiff PrePostMeanDiffPerfPulledInMedian PrePostMeanDiffPerfPulledInMean
 MontageMat2={'F0','F5','F30','M0','M5','M30','P0','P5','P30'};
 
 for ii=1:length(DatasetsIncluded)
  for mm=1:size(PerfSorted{ii},1)
     if ~isempty(MontAll{ii,mm})
     IndMont=find(contains( MontageMat2,MontAll{ii,mm}));
     IndSubj=str2num(DatasetsIncluded{ii}(1:2));
     
     
     AA=contains(MontAll,MontAll{ii,mm});
     sum(AA(ii,:)==1);
     IndxCol=find(AA(ii,:));
     for IndTrial=1:sum(AA(ii,:)==1)
        
     
      %This is the X% increase or decrease (the percent diffence) for
      %just the non-thresholded data.
     PercenDiffPerfPulledInMean{IndSubj,IndMont,IndTrial}=((PerfPulledInMean{ii,IndxCol(IndTrial),2}-PerfPulledInMean{ii,IndxCol(IndTrial),1})/PerfPulledInMean{ii,IndxCol(IndTrial),1})*100;
% %      GenDiffPerfPulledInMean{IndSubj,IndMont,IndTrial}=(PerfPulledInMean{ii,IndxCol(IndTrial),2}-PerfPulledInMean{ii,IndxCol(IndTrial),1});
% % 
% %      %Take the mean of the pre post and subtract it from during.
% %      PrePostMeanDiffPerfPulledInMean{IndSubj,IndMont,IndTrial}=(PerfPulledInMean{ii,IndxCol(IndTrial),2}-mean([PerfPulledInMean{ii,IndxCol(IndTrial),1},PerfPulledInMean{ii,IndxCol(IndTrial),3}]));
     

     end 
     end 
  end 
 end
    
 
 %%% 
 clear PercenMoreThanThresPerDiffMean  CoeffVariationRatioMean  VarianceRatioMean 
 clear GenDiffPerfPulledInMeanPooled GenDiffPerfPulledInMedianPooled PerDiffPerfPulledInMeanPooled PerDiffPerfPulledInMedianPooled
 clear PerDiffPerfPulledInMeanZPooled PerDiffPerfPulledInMedianZPooled PerDiffPerfPulledInStdPooled
 
  for ii=1:10
  for mm=1:length(MontageMat2)
     if ~isempty(MontAll{ii,mm})
     IndMont=find(contains( MontageMat2,MontAll{ii,mm}));
     IndSubj=str2num(DatasetsIncluded{ii}(1:2));
     
     
     AA=contains(MontAll,MontAll{ii,mm});
     sum(AA(ii,:)==1);
     IndxCol=find(AA(ii,:));

     
     PerDiffPerfPulledInMeanPooled{ii,mm}=nanmean([PercenDiffPerfPulledInMean{ii,mm,:}]);
     PerDiffPerfPulledInStdPooled{ii,mm}=nanstd([PercenDiffPerfPulledInMean{ii,mm,:}]);
     
     end 
  end 
  end
 

 
%% % ____________Plotting Performance as Imagesc with digits____________________________           
           
AA=vertcat(DatasetsIncluded{:});
NumUniqueSubjsNums=(unique(str2num(AA(:,1:2))));
NumConds=9;
clear  PlotThese
 PlotThese={cell2mat(PerDiffPerfPulledInMeanPooled)  ,-1,[-1,1].*max(max(abs(cell2mat(PerDiffPerfPulledInMeanPooled)))),...
            'Change(%)','Percent Change in Performance','Mean Percent Change in Performance During Stim PrePost-All-Imagesc'};
 
% PercenMoreThanThresPerDiffMean
 
 for rr=1:size(PlotThese,1)
 DatIn=(PlotThese{rr,1})*PlotThese{1,2};

 clear x y tttext
    figure;
    y = repmat(1:size(DatIn,1),size(DatIn,2),1)'; % generate y-coordinates
    x = repmat([1:NumConds],length(NumUniqueSubjsNums),1); % generate x-coordinates
    tttext=num2cell(round(DatIn,1));
    tttext = cellfun(@num2str, tttext, 'UniformOutput', false); % convert to string
    imagesc(1*((DatIn))); 
    text(x(:), y(:), tttext, 'HorizontalAlignment', 'Center','fontsize',14,'Color',[0 0 0])  %[249 166 2]./255 %[139,0,139]./255
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
 
 
 
 
 
%TODO  : Clean up boxplots 
 
 % Plot a boxplot of the agregate performance 
 
 
figure; 
boxplot(DatIn)
ylabel('Change(%)')
set(gca,'XTick',[1:9],'XTickLabels', MontageMat2)

DatIn_SE = nanstd(DatIn)./sqrt(sum(~isnan(DatIn),1));
figure; 
bar(nanmean(DatIn,1))
hold on
% plot(1:9, DatIn,'k*')
% hold on
errorbar(1:9, nanmean(DatIn),DatIn_SE*1.95); xlim([0.5 9.5])
set(gca,'XTick',[1:9],'XTickLabels', MontageMat2)
ylabel('Change(%)')


 clrs = cbrewer('seq', 'Blues', 20, 'pchip');
 bxpltclr =[241, 161, 4]./255;%[255,20,147]./255; [1, 146, 255]./255;%
 bxpltclr2 =[255, 233, 121]./255; %[0, 255, 249]./255;%
 bxpltmedclr =[165, 129, 5]./255;
 Meanlinesclr=[255, 96, 0]./255;%[0, 116, 63]./255; 
 
 figure
  plot(1:9, DatIn, '*'); hold on
  hh=boxplot(DatIn,'notch','off','BoxStyle','outline','Widths',0.2,'symbol',''); hold on;

    set(findobj(hh,'type','line','Tag','Upper Whisker','Tag','Lower Whisker'),'LineStyle','-');
    clear jj; jj=findobj(hh,'type','line');
%     set(jj([1 2 8 9 15 16]),'LineStyle','-');
    hh_out = findobj(hh,'Tag','Box'); 
        for j=1:length(hh_out) 
        patch(get(hh_out(j),'XData'),get(hh_out(j),'YData'),bxpltclr2,'EdgeColor',bxpltclr2,'FaceAlpha',.5 ,'LineStyle','-');
        end 
    set(hh,{'linew'},{3},{'color'},{bxpltclr})
 
 
 
 
 disp('done')          
           
           
