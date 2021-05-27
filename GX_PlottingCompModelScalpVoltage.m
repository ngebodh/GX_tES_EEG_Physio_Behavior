% %% GX_Plotting Computational Model Output Voltage
% Here we plot the scalp voltage that was predicted by an MRI derived computational
% head model. Data were derived by Comsol in a bespoke manner and saved to a .matfile 
%
%
% Written by: Nigel Gebodh
% Date: June 2020

%% Clear Residuals
clearvars
close all


%% Flags 

SveAllpics=1; 
closefigs=1;
Daterec='Model_05262020';


%% Setting Save Path 
  
  pathsave=strcat(['D:\GX\Results\ModelTopoplots\' Daterec '\']);
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
            
%             if ClearMatfiles==1
%             delete([pathsave 'FigOutput\*.mat'])
%             end
            
            %             rmdir([pathsave,'FigOutput'],'s'); %To erase the folder
            prefix = strcat(pathsave,'FigOutput','\');
        end
    end 


%% Load Data

load('CompModelEEG.mat')


%Get Topoplot Locations
Loc4Chans=['Standard-10-10-Cap33_V6.loc']; 
EEG.chanlocs = readlocs(Loc4Chans);



%% Plot Topoplots

topo_labels={'Frontal','Motor','Parietal'};
stim_amp={'1.0 mA','0.5 mA'};
map_lims=[25,18,18]

cc=1
for jj=1:2
 for ii=1:3,
figure;
datin= [[CompModelEEG{:,ii+2}]-CompModelEEG{33,ii+2}]/jj;
topoplot(datin,...
         EEG.chanlocs,'headrad',0.5,'plotrad',0.59,'style','map','electrodes','off','shading','interp','maplimits',[-1 1].*map_lims(ii));
         colorbar
         
title([topo_labels{ii} ' Stim Amp:' stim_amp{jj}])   
cc=cc+1;


fname=['ComputationalModelScalpVoltage' topo_labels{ii} ' Stim Amp-' num2str(jj)];
     set(gcf,'Name',fname,'Position',[1061         426         757         695]);
     

     if SveAllpics==1
         h = gcf;
         saveas(h,strcat(prefix,fname,'.fig'),'fig');
         print(h,'-dpng', [prefix,fname], '-r600');
         print(h,'-dpdf', [prefix,fname], '-r600');
         
     end
     if closefigs==1, close all,  end

 end
end 

     
     
