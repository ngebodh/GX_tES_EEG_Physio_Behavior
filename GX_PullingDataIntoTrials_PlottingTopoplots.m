%% GX_Pilot Data Pull out Perf to Segments
%
%Code written to chunck out the EEG and performance data into trials.
%The code chunks all the traisl regardless of stim type. 
%
% 
% Written by: Nigel Gebodh Dec. 2019
%
%Design
%{
Chunk the EEG data into trials
Chuck the performance data into trials
The stim types are indicated with the 'MontAll' variable


%}

%% Clear Residuals 
clear all
close all


%% Set Flags

SaveTrialedData=1; %Save the trialed data 0-No don't save, 1- Yes, save.
PlotTopoplots=1; %Defines if you want topoplots for each stim trial plotted
SveAllpics=1; 
closefigs=1;
ClearMatfiles=1;
PhaseDesignAll=1;%Run both phase 1 and 2? 1-Yes, 0-No. 
PhaseDesign=2;   %If run specific phase, select phase, either 1 or 2 
Daterec='07022020';


%% Create Results folder        

if PhaseDesignAll==0
    if PhaseDesign==2
        pathsave=strcat(['D:\GX Project\Results\DataChunkedtoTrials_Phase2_' Daterec '\']);
    else 
       pathsave=strcat(['D:\GX Project\Results\DataChunkedtoTrials_Phase1_' Daterec '\']);
    end 
    DesignLoop=PhaseDesign;
else
    pathsave=strcat(['D:\GX Project\Results\DataChunkedtoTrials_AllPhases_' Daterec '\']);
    DesignLoop=1:2;
end 


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
            
            if ClearMatfiles==1
            delete([pathsave 'FigOutput\*.mat'])
            end
            
            %             rmdir([pathsave,'FigOutput'],'s'); %To erase the folder
            prefix = strcat(pathsave,'FigOutput','\');
        end
    end 

for PhaseDesign= DesignLoop             
        
%% Define Data Locations and Files to Look At

%This is where all the data are stored
DataLoc='D:\GX Project\Data\';
%These are the files that we want to look at
if PhaseDesign==2
    DatasetsIncluded={...
                        '1101','1102',...
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
    
elseif PhaseDesign==1
    DatasetsIncluded={...
        '0101','0102','0103','0104',...
        '0201','0202',...
        '0301','0302','0303',...
        '0401','0402','0403',...
        '0501','0504','0505',...
        '0601','0602','0603',...
        '0701','0702','0703',...
        '0801','0802','0803',...
        '0901','0902','0903',...
        '1001','1002','1003'};
end

if PhaseDesign==2
    clear MontAll
    GX_SubjectMontages_TaskDesign2
else 
    %This needs to be fixed. I need to pass the proper montage matrix
    MontageMat2={'F0','F5','F30','M0','M5','M30','P0','P5','P30'};
   load(['D:\GX Project\Results\DataChunkedtoTrials_12162019\FigOutput\','MontAllCompiled.mat'])
end


for ii=1:length(DatasetsIncluded)
    disp(['Now running ' DatasetsIncluded{ii} '.....' ])
SelectedFle=strcat(DataLoc,DatasetsIncluded{ii},'\',DatasetsIncluded{ii},'\','ptracker-',DatasetsIncluded{ii},'.csv');
% opts.Sheet = '2007';
opts.SelectedVariableNames = [3,11]; 
opts.DataRange = '';%'2:11';

% preview(SelectedFle,opts)
% preview('D:\GX Project\Data\0101\0101\ptracker-0101.csv',opts)
% D:\GX Project\Data\0101\0101\ptracker-0101.csv
% 
% M = readmatrix(SelectedFle,opts)
filename = SelectedFle;%'D:\GX Project\Data\0101\0101\ptracker-0101.csv';
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

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

desiredFs = 100; 
ScreenFs = 60; 
ptrackerPerf{ii}=resample(ptrackerData{ii}(:,2),ptrackerData{ii}(:,1),desiredFs,desiredFs,ScreenFs);
ptrackerTime{ii}=[[0:length(ptrackerPerf{ii})-1]./desiredFs]'; 


clear ptrackerData
tic


toc
 
 
 %% Getting the EEG 
Chans=[1:32];
numcount=ii;
        %Define where each EEG file is 
        GetFilesFrom=strcat('D:\GX Project\Data\' ,DatasetsIncluded{ii},'\');
        if ~exist( GetFilesFrom)
            numcount= numcount+1; %Added to keep the order of existing files 
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

DataEEG{numcount}=EEG.data([Chans],:);

    %There was an extra blank trigger of data for subject 1401. 
    %This removes that blank trigger to make concatination easy. 
    if strmatch(DatasetsIncluded{ii},'1401')
        EEG.triggers(1)=[];
    end 

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


Loc4Chans=['Standard-10-10-Cap33_V6.loc'];
EEG.chanlocs = readlocs(Loc4Chans);
 
    
        %Baseline and drift correct the data
        Adj_topoly_Each{numcount}=DataEEG{numcount};
        baselineT=   find(t{numcount}>1 & t{numcount}<5);  %was 60 find(t>-10 & t<50);  % for baseline correction 
        BLamp = mean(Adj_topoly_Each{numcount}(:,baselineT),2); % record baseline amplitude (t<0) for each channel
        BLcorDC{numcount} = Adj_topoly_Each{numcount}(:,:) - repmat(BLamp,[1,length(t{numcount})]); % baseline correction
        BLcorDC{1,numcount}(33:34,:)=EEG.data(33:34,:);

% Quickly Look for potentially bad electrodes
for nn=1:32 
    PP(nn,:)=sum(( BLcorDC{1,numcount}(nn,1:2000*60*10)).^2); 
%     PP(nn,:)=sum((DSamp.data(nn,:)).^2);
end


% % % % 
% % % % %% Create Montage File
% % % % if PhaseDesign==1
% % % % 
% % % %     %Create a matrix of montages. 
% % % %     MatFiles=dir(fullfile(GetFilesFrom, '*.mat'));
% % % %     if ~isempty(MatFiles)
% % % %         load(strcat(GetFilesFrom,MatFiles.name),'Montages');
% % % %         MontHold=repmat(Montages,4,1);
% % % %         Mont=upper(MontHold(:)');
% % % % 
% % % %     end
% % % %      
% % % %  
% % % %  % _____________Looking At Each Stimulation Trial__________________________
% % % %  
% % % %  %Find all the Stim on Triggers
% % % % % if  strcmp(DatasetsIncluded{ii},'1401')==1, AllEventsCode{ii}(1)=[]; end 
% % % % 
% % % %  Evnt_Stimstrt=AllEvents{ii}(find(str2num(vertcat(AllEventsCode{ii}{:}))==16));
% % % %  
% % % %  if DatasetsIncluded{ii}=='0102'
% % % %      Evnt_Stimstrt=Evnt_Stimstrt(1:end-1)
% % % %      Mont=Mont(1:length( Evnt_Stimstrt))
% % % %      MontAll(ii,1:length(Mont))=Mont;
% % % %  elseif DatasetsIncluded{ii}=='0101'
% % % %      Evnt_Stimstrt=Evnt_Stimstrt(1:end-1)
% % % %      Mont=Mont(1:length( Evnt_Stimstrt))
% % % %      MontAll(ii,1:length(Mont))=Mont;
% % % %  else 
% % % %       MontAll(ii,1:length(Mont))=Mont;
% % % %  end 
% % % %  
% % % %  clear Emp
% % % %  if sum(cellfun(@isempty,{MontAll{ii,:}}))>0
% % % %         Emp=find(cellfun(@isempty,{MontAll{ii,:}})); 
% % % %         for tt=1:length(Emp)
% % % %             MontAll{ii,Emp(tt)}='';
% % % %         end 
% % % %  end 
% % % % 
% % % % end
%%

    

%%%_______________________Fix Bad Electrodes_______________________________

        if strcmp(DatasetsIncluded{ii},'1501') || strcmp(DatasetsIncluded{ii},'1502')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [3]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end

        if strcmp(DatasetsIncluded{ii},'0202')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [26]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end
       
        if strcmp(DatasetsIncluded{ii},'0402')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [21]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end 
        
        if strcmp(DatasetsIncluded{ii},'0403')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [21,23]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end        
        
        if strcmp(DatasetsIncluded{ii},'0501')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [5]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end
        
        
        if strcmp(DatasetsIncluded{ii},'0601')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [5]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end

        if strcmp(DatasetsIncluded{ii},'0701')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [26]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end        
        
        if strcmp(DatasetsIncluded{ii},'0702')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [26,30]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end  
        
        if strcmp(DatasetsIncluded{ii},'0703')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [27]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end                  
        
        if strcmp(DatasetsIncluded{ii},'0801')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [5,14]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end
        
        if strcmp(DatasetsIncluded{ii},'0802')|| strcmp(DatasetsIncluded{ii},'0803')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [5, 14]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end
        
        
        if strcmp(DatasetsIncluded{ii},'0901') || strcmp(DatasetsIncluded{ii},'0902') || strcmp(DatasetsIncluded{ii},'0903')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [3]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end
        
        if strcmp(DatasetsIncluded{ii},'1002')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [14,25]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end
        
        if strcmp(DatasetsIncluded{ii},'1003')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [7,25]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end

        if strcmp(DatasetsIncluded{ii},'2002') || strcmp(DatasetsIncluded{ii},'2001')
            dataOut = fillBadChannelsFast( BLcorDC{1,numcount}(1:32,:)', [25]);
            BLcorDC{1,numcount}(1:32,:)=dataOut';
        end      
        
        

clear DataEEG Adj_topoly_Each Samp baselineT
% EEG.data=[];
EEG.time=[];


%{

 %% ____Look at Whole Recording __________________
 
%  [b,a]=butter(3,[[0.5, 40]/1000],'bandpass')
%  figure; 
%  Addthis=AllEventsTime{ii}(1);
%  plot(t{ii},filtfilt(b,a,BLcorDC{ii}(15,:))); 
%  ylabel(['EEG Voltage(\muV)'])
% 
%  
%  hold on; 
%  yyaxis right
%  plot(ptrackerTime{ii}(:,1)+ Addthis, ptrackerPerf{ii}(:,1));
%  xlabel(['Time (Sec)'])
%  ylabel(['Performance'])
%  title(['Whole Recording-Channel C3'])
% 
%     axis tight  
%     fname=[ 'Subj-' DatasetsIncluded{ii} '-Whole Session Rec'];
%     set(gcf,'Name',fname,'Position',[521         436        1062         559],'PaperPositionMode','auto')
%      if SveAllpics==1
%             h = gcf;
%             saveas(h,strcat(prefix,fname,'.fig'),'fig');
%             saveas(h,strcat(prefix,fname,'.png'),'png');
%             print(h,'-dpng', [prefix,fname], '-r600');
%             print(h,'-depsc', [prefix,fname], '-r600');
%             print(h,'-dpdf', [prefix,fname], '-r600');
%      end 
%         
%      if closefigs==1, close all,  end
 
 %}
 
 
%% 


 
 % _____________Looking At Each Stimulation Trial__________________________
 
 %Find all the Stim on Triggers

 Evnt_Stimstrt=AllEvents{ii}(find(str2num(vertcat(AllEventsCode{ii}{:}))==16));
 
 
 EEGout.fs=fs{ii};
Perfout.fs=desiredFs;
if PhaseDesign==1
    TrialTotals=length(Evnt_Stimstrt);
else
    TrialTotals=20;
end

 
 
 if PhaseDesign==1; trls=4; else, trls=20;end 
     MatFiles=dir(fullfile(GetFilesFrom, '*.mat'));
     if ~isempty(MatFiles)
         load(strcat(GetFilesFrom,MatFiles.name),'Montages');
         MontHold=repmat(Montages,trls,1);
         Mont=upper(MontHold(:)');
         %         MontAll{ii}= Mont;
     end
     
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
 
 
 
 
 
 Evnt_Stimstrt2=(Evnt_Stimstrt-AllEvents{ii}(1)).*(desiredFs/fs{ii});%Adjust for delay between EEG and task start
 Evnt_BlockStart=AllEvents{ii}(find(str2num(vertcat(AllEventsCode{ii}{:}))==2));
 Evnt_BlockStart2=(Evnt_Stimstrt-AllEvents{ii}(1)).*(desiredFs/fs{ii});%

 


 for mm=1:TrialTotals
   
     if mm<=length(Evnt_Stimstrt)
%EEG data including physio  
EEGout.PreStim{mm,1}   =BLcorDC{ii}(1:34   ,    Evnt_Stimstrt(mm)-(30*fs{ii}):Evnt_Stimstrt(mm));% 
EEGout.PreStim{mm,2}   =MontAll{ii,mm};

EEGout.DuringStim{mm,1}=BLcorDC{ii}(1:34   ,    Evnt_Stimstrt(mm)+(5.25*fs{ii}):Evnt_Stimstrt(mm)+((5.25+30)*fs{ii}));
EEGout.DuringStim{mm,2}=MontAll{ii,mm};

EEGout.PostStim{mm,1}  =BLcorDC{ii}(1:34   ,    Evnt_Stimstrt(mm)+(41.5*fs{ii}):Evnt_Stimstrt(mm)+((41.5+30)*fs{ii}));%Here we did 41 because there seems to be some spectral lekage into the post stim
EEGout.PostStim{mm,2}  =MontAll{ii,mm};




%Performance 
Perfout.PreStim{mm,1}        =ptrackerPerf{ii}(Evnt_Stimstrt2(mm)-(30*desiredFs):Evnt_Stimstrt2(mm),1); 
Perfout.PreStim{mm,2}        =MontAll{ii,mm};

Perfout.DuringStim{mm,1}     =ptrackerPerf{ii}(Evnt_Stimstrt2(mm)+(5.25*desiredFs):Evnt_Stimstrt2(mm)+((5.25+30)*desiredFs),1); 
Perfout.DuringStim{mm,2}     =MontAll{ii,mm};

Perfout.PostStim{mm,1}       =ptrackerPerf{ii}(Evnt_Stimstrt2(mm)+(41.5*desiredFs):Evnt_Stimstrt2(mm)+((41.5+30)*desiredFs),1); 
Perfout.PostStim{mm,2}       =MontAll{ii,mm};

    
     
     if PlotTopoplots==1
         
        
     figure;
     clear h1 pxx
     for jj=1:3
         clear datin datpull datmont locs ElecSaturated
         ElecSaturated=[];ElecNoisy=[];
         subplot(2,4,jj)
         magamp=1e2;
         
         datmont=EEGout.PreStim{mm,2}; %Data montage
             if lower(datmont(1))=='f'
                 channum=9;  %FC5
             elseif lower(datmont(1))=='m'
                 channum=15; %C3
             elseif lower(datmont(1))=='p'
                 channum=20; %CP5
             end
         
         
         if jj==1 %Pre Stim
             datpull=EEGout.PreStim{mm,1}(1:32,:);
             datin=mean(datpull,2);
             labels='Pre Stim';

             
         elseif jj==2 %During Stim 
             %Normalize to  the pre stimulation data
             datpull=EEGout.DuringStim{mm,1}(1:32,:)-mean(EEGout.PreStim{mm,1}(1:32,:),2);
          
             
            if str2num(datmont(2))==0
                locs=1:length(datpull);
            else
                [~,locs]=findpeaks(1*datpull(channum,:), 'MinPeakDistance',50,'MinPeakProminence',500,'Annotate','extents');
            end
             
             magamp=2e4;
             datin=mean(datpull(:,locs),2);
             labels='During Stim';
             
         elseif jj==3 %Post Stim 
             datpull=EEGout.PostStim{mm,1}(1:32,:);
             datin=mean(datpull,2);
             labels='Post Stim';
             
         end
         
         
        
         topoplot((datin)./1000,...
             EEG.chanlocs,'headrad',0.5,'plotrad',0.59,'style','map','electrodes','off','shading','interp');
             title(['Sub ',DatasetsIncluded{ii}, ' ',datmont,' Trial ', num2str(mm)]);
             colorbar
             
         %Plot all EEG electrodes             
         subplot(2,4,jj+4)
         seglen=length(datpull(1,:));
         satthreh=seglen*0.10;
         xlen=1:seglen;
         sat=1; eenoisy=1;
         for ee=1:32; 
             plot(xlen, repmat(-(magamp*ee),1, length(datpull(1,:))),':k' );
             hold on;
             plot(xlen, (datpull(ee,:)-mean(datpull(ee,:)))-(magamp*ee)  ); hold on;axis tight; 
             yticks(ee)=-(magamp*ee);
           
                 %Check for saturation    
                 if sum(diff(datpull(ee,:))==0)>satthreh  
                     ElecSaturated{sat}=EEG.chanlocs(ee).labels;
                     sat=sat+1;
                 end 

                 %Check for noise
                 clear pxx2 pxx16 ff; 
                 [pxx16, ~]=pwelch(datpull(16,:), 1000,500,1000,2000); 
                 [pxx2, ~]=pwelch(datpull(ee,:), 1000,500,1000,2000);
                 bndpw_Cz=bandpower(pxx16,2000,[55 65]);
                 bndpw_EE=bandpower(pxx2,2000,[55 65]);

                 if (bndpw_EE/bndpw_Cz)>3.5
                     ElecNoisy{eenoisy}=EEG.chanlocs(ee).labels;
                     eenoisy=eenoisy+1;
                 end 

         end
         xlabel('Samples');
         set(gca,'YTick',[fliplr(yticks)],'YTickLabel',fliplr({EEG.chanlocs.labels}));
         
         
         if ~isempty(ElecSaturated)
             %{
         %Code snip from Frank Zalkow
         %Source: https://www.mathworks.com/matlabcentral/answers/28537-set-the-yticklabel-to-different-colors
             %}
             % Get the current tick labels
             ticklabels = get(gca,'YTickLabel');
             % Create an empty list to be filled in
             ticklabels_new = cell(size(ticklabels));
             for i = 1:length(ticklabels)
                 %If theres saturation set the color to red
                 if sum(contains(ElecSaturated, ticklabels{i}))>=1
                     ticklabels_new{i} = ['\color{red} ' ticklabels{i}];
                 else
                     %If no Saturation leave as black
                     ticklabels_new{i} = ['\color{black} ' ticklabels{i}];
                 end
             end
         elseif ~isempty(ElecNoisy)
             % Get the current tick labels
             ticklabels = get(gca,'YTickLabel');
             % Create an empty list to be filled in
             ticklabels_new = cell(size(ticklabels));
             for i = 1:length(ticklabels)
                 %If theres noise set the color to red
                 if sum(contains(ElecNoisy, ticklabels{i}))>=1
                     ticklabels_new{i} = ['\color{magenta} ' ticklabels{i}];
                elseif ~isempty(ElecSaturated)
                    if sum(contains(ElecSaturated, ticklabels{i}))>=1
                     ticklabels_new{i} = ['\color{red} ' ticklabels{i}];
                    end 
                else
                     %If no nosie leave as black
                     ticklabels_new{i} = ['\color{black} ' ticklabels{i}];
                 end
             end
           
         else 
             ticklabels_new = get(gca,'YTickLabel');
         end
         
         % set the tick labels
         set(gca, 'YTickLabel', ticklabels_new);
         
         title([{['Sub ',DatasetsIncluded{ii}, ' ',labels]}]);
         
         %Plot the welch spectrum
         subplot(2,4,[4,8])
%          [pxx,ff]=pwelch(datpull(channum,:),2000,1000, 2000,2000);
         [pxx{jj},ff]=pwelch(datpull(16,:),2000,1000, 2000,2000);
         h1{jj,:}=plot(ff, db(pxx{jj}),'linewidth',2);
         hold on
         ylabel(['PSD (dB/Hz)']);
         xlabel(['Frequency(Hz)']);
         title([{['Sub ',DatasetsIncluded{ii}, ' ',datmont,' Trial ', num2str(mm)]},{['Chan: ' EEG.chanlocs(16).labels]}]);
         set(gca, 'XScale', 'log');
         xlim([0,100]);
         
         
     end
     hh=legend([h1{1:3}],{['\color[rgb]{' num2str(h1{1}.Color) '}Pre'],...
                          ['\color[rgb]{' num2str(h1{2}.Color) '}During'],...
                          ['\color[rgb]{' num2str(h1{3}.Color) '}Post']},'Location','Northwest');
     
     fname=['Topo DC- Sub ',DatasetsIncluded{ii}, '-',datmont,' Trial ', num2str(mm)];
     set(gcf,'Name',fname,'Position',[391  153 1587 580]);
     
     if SveAllpics==1
         h = gcf;
         saveas(h,strcat(prefix,fname,'.fig'),'fig');
         saveas(h,strcat(prefix,fname,'.png'),'png');
         saveas(h,strcat(prefix,fname,'.pdf'),'pdf');
         print(h,'-dpng', [prefix,fname], '-r600');
         
     end
     if closefigs==1, close all,  end
     
     
     
     
     figure
     for jj=1:3
       %Plot the welch spectrum

         h1{jj,:}=plot(ff, db(pxx{jj}),'linewidth',2);
         hold on
         ylabel(['PSD (dB/Hz)']);
         xlabel(['Frequency(Hz)']);
         title([{['Sub ',DatasetsIncluded{ii}, ' ',datmont,' Trial ', num2str(mm)]},{['Chan: ' EEG.chanlocs(16).labels]}]);
         set(gca, 'XScale', 'log');
         xlim([0,100]);
         
         
     end
     ylim([-50 150])
     hh=legend([h1{1:3}],{['\color[rgb]{' num2str(h1{1}.Color) '}Pre'],...
                          ['\color[rgb]{' num2str(h1{2}.Color) '}During'],...
                          ['\color[rgb]{' num2str(h1{3}.Color) '}Post']},'Location','Northwest');
     
    fname=['Topo DC- Sub ',DatasetsIncluded{ii}, '-',datmont,' Trial ', num2str(mm),'-JustSpectrum'];
    set(gcf,'Name',fname,'Position',[1000        1107         257         231]);

    if SveAllpics==1
        h = gcf;
        saveas(h,strcat(prefix,fname,'.fig'),'fig');
        saveas(h,strcat(prefix,fname,'.png'),'png');
        saveas(h,strcat(prefix,fname,'.pdf'),'pdf');
        print(h,'-dpng', [prefix,fname], '-r600');

    end
    if closefigs==1, close all,  end
    
    
     
    
    
    
    
     figure;
     clear h1 
     
         clear datin datpull datmont locs ElecSaturated
                  
         magamp=1e2;
         
         datmont=EEGout.PreStim{mm,2}; %Data montage
             if lower(datmont(1))=='f'
                 channum=5;  %F3;
                 map_lims=[25];
             elseif lower(datmont(1))=='m'
                 channum=15; %C3
                 map_lims=[18];
             elseif lower(datmont(1))=='p'
                 channum=20; %CP5
                 map_lims=[18];
             end
         
         
             %During Stim 
             %Normalize to  the pre stimulation data
             datpull=EEGout.DuringStim{mm,1}(1:32,:)-mean(EEGout.PreStim{mm,1}(1:32,:),2);
             
             if str2num(datmont(2))==0
                  locs=1:length(datpull);
             else
                [~,locs]=findpeaks(1*datpull(channum,:), 'MinPeakDistance',50,'MinPeakProminence',500,'Annotate','extents');
             end 

             magamp=2e4;
             datin=mean(datpull(:,locs),2);
             labels='During Stim';

         
         
        
         topoplot((datin)./1000,...
             EEG.chanlocs,'headrad',0.5,'plotrad',0.59,'style','map','electrodes','off','shading','interp','maplimits',[-1 1].*map_lims);
             title(['Sub ',DatasetsIncluded{ii}, ' ',datmont,' Trial ', num2str(mm)]);
             colorbar
             

     fname=['Topo DC- Sub ',DatasetsIncluded{ii}, '-',datmont,' Trial ', num2str(mm) ,'-Justtopo'];
     set(gcf,'Name',fname,'Position',[1061         426         757         695]);
     
     if SveAllpics==1
         h = gcf;
         saveas(h,strcat(prefix,fname,'.fig'),'fig');
%          saveas(h,strcat(prefix,fname,'.png'),'png');
%          saveas(h,strcat(prefix,fname,'.pdf'),'pdf');
         print(h,'-dpng', [prefix,fname], '-r300');
         print(h,'-dpdf', [prefix,fname], '-r600');
         
     end
     if closefigs==1, close all,  end
     
     
     
     
     
     
     
     
     
     
     end 
     
     end 
     
 end 

 
BLcorDC{ii}=[];
Evnt_StimstrtAll{ii,:}=Evnt_Stimstrt;

SelectedFileNum=DatasetsIncluded{ii};

if SaveTrialedData==1
save(strcat(prefix,DatasetsIncluded{ii},'_EEG.mat'), 'EEGout','Perfout','EEG','MontAll','SelectedFileNum','-v7');
end 
 
end  % End of each subject loading 



end %End of Phase design














%%%%%%%%%%%%%    END    %%%%%%%%%%%%%