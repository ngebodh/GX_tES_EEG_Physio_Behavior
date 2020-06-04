%% GX_PlottingDemographicInfo
%This script was written to plot the demographic and questionnaire data
%that were collected for the GX project
%
% Written by: Nigel Gebodh
% Date: May 2020
%
% Requirements:
% -GX_Subject Info & Behavioral Data.xlsx
%
%Updates:
%-June 2020: 
%   Added additional figure for questionnaires. 
%


%% Clear residuals 

clear all
close all


SveAllpics=0; 
closefigs =0; 


 pathsave=strcat('D:\GX\Results\DemographicInfo\');

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











%% Pull in the data
%{
Here we pull in the data from the spreadsheets that they were put into


%}
%Where is the spreadsheet? 
DatLoc='D:\GX\Data\GX_Subject Info & Behavioral Data.xlsx';

%Get the demographic data
clear AA
AA=importfile(DatLoc, 'Demographic',[2 26]);
AA(:,1)=[];%Clear the blank column

%Get the questionnaire data
clear BB
BB=importfile(DatLoc, 'All Behavioral',[2 90]);

clear BB2
BB2=BB; 



%Exclude subject 'KP'    
Indx=find([AA{:,1}] ~= 17 );
    Gen=vertcat(AA{Indx,5}); 
    Age=str2double(vertcat(AA{Indx,4})); 
    Height=str2double(vertcat(AA{Indx,6})); 
    Weight=str2double(vertcat(AA{Indx,7})); 

    clear ind
    ind=cellfun(@isnumeric,{AA{:,4}});
    
%Pull out the data that we will be plotting. Refer to spreadsheet for column names.      
SubjInit=cellfun(@char,{AA{ind,2}},'UniformOutput',false)';
Gens=cellfun(@char,{AA{ind,5}},'UniformOutput',false)';
AllSess=vertcat(BB{find(cellfun(@isnumeric,{BB{:,2}})==1),6});
IndAllSess=find(cellfun(@isnumeric,{BB{:,2}})==1);

for ii=1:length(SubjInit)
gg=[vertcat(AllSess{:})==SubjInit{ii,:}];
IndInsert=IndAllSess(find((gg(:,1).*gg(:,2))==1));

    for jj = 1:length(IndInsert)
     BB2{IndInsert(jj),80}=Gens(ii,:);
    end 
end 

 % Gather Questionnaire Data
validInd=find(cellfun(@isnumeric,{BB{:,1}})==1); 
ExInd=validInd(vertcat(BB{validInd,1})==17);



for ii=1:size(BB,2)
    %Make these rows empty
    BB2{ExInd,ii}=[];
end 

AllSessMF=vertcat(BB2{find(cellfun(@isnumeric,{BB{:,2}})==1),80});    
    
    
%% Put everything together and plot it
    
    
%%% M/F by Age
x1 = Age(Gen=='M');
x2 = Age(Gen=='F');
x3 = [x1', x2']';
x = [x1; x2; x3];
g = [zeros(length(x1), 1); ones(length(x2), 1); 2*ones(length(x3), 1)]; 


%Setting up Colors
        Fclr =[123, 189, 69]./255;
        Fclr2 =[202, 229, 181]./255;
        Mclr =[0, 148, 210]./255;
        Mclr2 =[153, 212, 237]./255;
        
        Allclrs=[164, 85, 255]./255;
        bxpltclr= [4, 170, 208]./255;
        bxpltclr2 = [4, 170, 208]./255;
        hisclr = [0.5 0.5 0.5];
figure; 

    %Setting up location of axes in figure
%      hp1 = uipanel('position',[0.1331*1.05 0.1501-0.1 0.1535*1.45 0.7749+0.1]);   
%      hp2 = uipanel('position',[0.3361*1.05 0.1503-0.1 0.1566*1.45 0.7747+0.1]);
%      hp3 = uipanel('position',[0.5422*1.05 0.1503-0.1 0.1566*1.45 0.7747+0.1]);
%      hp4 = uipanel('position',[0.7484*1.05 0.1503-0.1 0.1566*1.45 0.7747+0.1]);

     hp1 = uipanel('position',[0.005706742627348,0.0501,0.222575000000004,0.8749]);   
     hp2 = uipanel('position',[0.241913042895449,0.0503,0.227070000000003,0.8747]);
     hp3 = uipanel('position',[0.489953431635396,0.0503,0.227070000000001,0.8747]);
     hp4 = uipanel('position',[0.737026434316362,0.0503,0.22707,0.8747]);

%     subplot(1,4,1)
%     axes('Parent',hp1)
        scatterhist(nan(1,length([x3])),[x3], 'Kernel','on','Marker','none','Direction','out','Location','NorthEast','Color',hisclr,'Parent',hp1);
        hold on
        plot(1, x1, '*', 'linewidth', 3,'Color', Mclr2);
            
        plot(2, x2, '*', 'linewidth', 3,'Color', Fclr2);
            
        
        plot(3, x1, '*', 'linewidth', 3,'Color', Mclr2);
        plot(3, x2, '*', 'linewidth', 3,'Color', Fclr2);
          
       
        hh=boxplot(x, g,'BoxStyle','outline','Widths',0.2,'symbol','');
        

        
            set(findobj(hh,'type','line','Tag','Upper Whisker','Tag','Lower Whisker'),'LineStyle','-');
            clear jj; jj=findobj(hh,'type','line');
            set(jj([1 2 8 9 15 16]),'LineStyle','-');
            set(jj([1 2 3 4 5]),'Color',Mclr);
            set(jj([8 9 10 11 12]),'Color',Fclr);
            set(jj([15 16 17 18 19]),'Color',[0.5 0.5 0.5]);
            hh_out = findobj(hh,'Tag','Box');
%             for j=1:length(hh_out)
%                 patch(get(hh_out(j),'XData'),get(hh_out(j),'YData'),bxpltclr2,'EdgeColor',bxpltclr2,'FaceAlpha',.25 ,'LineStyle','-');
%             end
            set(hh,{'linew'},{2}) %,{'color'},{bxpltclr}
            
             plot(1, mean(x1), '^', 'linewidth', 1,'Color',  Mclr,'MarkerFaceColor',Mclr);
             plot(2, mean(x2), '^', 'linewidth', 1,'Color',  Fclr,'MarkerFaceColor',Fclr);
             plot(3, mean(x3), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[0 0 0]);
             
       
        
        xlabel(['Gender']);
        ylabel(['Age(years)']) ;
        set(gca, 'XTickLabels',{'Males', 'Females', 'All Subjects'});
        set(gca, 'box','off','XAxisLocation','bottom','YAxisLocation','left');      
        
        
        
        
%%% M/F by Weight and Height        
clear x1 x2
x1 = Height(Gen=='M');
y1 = Weight(Gen=='M'); 
x2 = Height(Gen=='F'); 
y2 = Weight(Gen=='F'); 


%     subplot(1,4,2)
    axes('Parent',hp2);
        scatterhist([x1;x2],[y1;y2], 'Kernel','on','Marker','none','Direction','out','Location','NorthEast','Color',hisclr,'Parent',hp2);
        hold on
        plot(x1, y1, '*', 'linewidth', 3,'Color', Mclr2);
        hold on
        plot(x2, y2, '*', 'linewidth', 3,'Color', Fclr2);
        hold on
        
             plot(mean(x1),mean(y1), '^', 'linewidth', 1,'Color',  Mclr,'MarkerFaceColor',Mclr);
             hold on
             plot(mean(x2),mean(y2), '^', 'linewidth', 1,'Color',  Fclr,'MarkerFaceColor',Fclr);
             hold on
%              plot(mean([x1;x2]),mean([y1;y2]), '^', 'linewidth', 1,'Color', [0 0 0]) 
             
                     hold on
                line(mean([x1;x2]) + [-std([x1;x2]) std([x1;x2])], [1, 1].*mean([y1;y2]),'Color', [0.5 0.5 0.5]);
                hold on
                line([1, 1].*mean([x1;x2]),mean([y1;y2]) + [-std([y1;y2]) std([y1;y2])],'Color', [0.5 0.5 0.5]);
                hold on 
                plot(median([x1;x2]),median([y1;y2]), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
                plot(mean([x1;x2]),mean([y1;y2]), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[0 0 0]);
                
        
        xlabel(['Height (cm)']);
        ylabel(['Weight (kg)']);
        set(gca, 'box','off','XAxisLocation','bottom','YAxisLocation','left');         
        
        
        
    
%Pull out the data that we'll be using. Refer to spreadsheet for column names.   
clear indM; indM= IndAllSess(vertcat(AllSessMF{:})=='M');
clear indF; indF= IndAllSess(vertcat(AllSessMF{:})=='F');
SleepHrs  = double(vertcat(BB2{:,12}));
    SleepHrsM= double(vertcat(BB2{indM,12}));
    SleepHrsF= double(vertcat(BB2{indF,12}));

SleepQual = double(vertcat(BB2{:,13}));
    SleepQualM= double(vertcat(BB2{indM,13}));
    SleepQualF= double(vertcat(BB2{indF,13}));
    
KSSpre    = double(vertcat(BB2{:,47}));
    KSSpreM    = double(vertcat(BB2{indM,47}));
    KSSpreF    = double(vertcat(BB2{indF,47}));
    
KSSpost   = double(vertcat(BB2{:,53}));
    KSSpostM   = double(vertcat(BB2{indM,53}));
    KSSpostF   = double(vertcat(BB2{indF,53}));


    
%     subplot(1,4,3)
%     axes('Parent',hp3)
%         plot(SleepHrs, SleepQual, '*', 'linewidth', 3,'Color', Allclrs)
        scatterhist(SleepHrs,SleepQual, 'Kernel','on','Marker','none','Direction','out','Location','NorthEast','Color',hisclr,'Parent',hp3);
        hold on
        plot(SleepHrsM, SleepQualM, '*', 'linewidth', 3,'Color', Mclr2);
        hold on
        plot(SleepHrsF, SleepQualF, '*', 'linewidth', 3,'Color', Fclr2);
            hold on
            plot(nanmean(SleepHrsM), nanmean(SleepQualM), '^', 'linewidth', 1,'Color', Mclr,'MarkerFaceColor',Mclr);
            hold on
            plot(nanmean(SleepHrsF), nanmean(SleepQualF), '^', 'linewidth', 1,'Color', Fclr,'MarkerFaceColor',Fclr);
            
        hold on
        line(nanmean(SleepHrs) + [-nanstd(SleepHrs) nanstd(SleepHrs)], [1, 1].*nanmean(SleepQual),'Color', [0.5 0.5 0.5]);
        hold on
        line([1, 1].*nanmean(SleepHrs),nanmean(SleepQual) + [-nanstd(SleepQual) nanstd(SleepQual)],'Color', [0.5 0.5 0.5]);
        hold on 
        plot(nanmedian(SleepHrs), nanmedian(SleepQual), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
        plot(nanmean(SleepHrs), nanmean(SleepQual), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[0 0 0]);

        
        xlabel(['Sleep Hours']);
        ylabel(['Sleep Quality'])  ;  
        xlim([2.5 10]);
        ylim([2.5 10]);
%         set(gca,'YTick',[0 , 5, 10], 'YTickLabels',{'Worst', 'Normal', 'Best'})
        set(gca, 'box','off','XAxisLocation','bottom','YAxisLocation','left'); 
        
        
        
        
        
%     subplot(1,4,4)
       scatterhist(KSSpre,KSSpost, 'Kernel','on','Marker','none','Direction','out','Location','NorthEast','Color',hisclr,'Parent',hp4);
        
        hold on 
        plot(KSSpreM, KSSpostM, '*', 'linewidth', 3,'Color', Mclr2);
        hold on
        plot(KSSpreF, KSSpostF, '*', 'linewidth', 3,'Color', Fclr2);
        hold on
                plot(nanmean(KSSpreM), nanmean(KSSpostM), '^', 'linewidth', 1,'Color', Mclr,'MarkerFaceColor',Mclr);
                hold on
                plot(nanmean(KSSpreF), nanmean(KSSpostF), '^', 'linewidth', 1,'Color', Fclr,'MarkerFaceColor',Fclr);
                hold on
        
        hold on
        line(nanmean(KSSpre) + [-nanstd(KSSpre) nanstd(KSSpre)], [1, 1].*nanmean(KSSpost),'Color', [0.5 0.5 0.5]);
        hold on
        line([1, 1].*nanmean(KSSpre),nanmean(KSSpost) + [-nanstd(KSSpost) nanstd(KSSpost)],'Color', [0.5 0.5 0.5]);
        plot(nanmedian(KSSpre), nanmedian(KSSpost), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
        plot(nanmean(KSSpre), nanmean(KSSpost), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[0 0 0]);
       

        xlabel(['KSS Pre']);
        ylabel(['KSS Post']) ;   
        xlim([0.75 9.25]);
        ylim([0.75 9.25]);
        
        fname='All Demo and PrePost Questionnaire';
        hh=gcf;
        set(hh,'name',fname,'Position', [ -2064         518        1865         466]);
        set(gca, 'box','off','XAxisLocation','bottom','YAxisLocation','left'); 
        
        
         if SveAllpics==1
               h = gcf;
               saveas(h,strcat(prefix,fname,'.fig'),'fig');
               saveas(h,strcat(prefix,fname,'.png'),'png');
               print(h,'-dpng', [prefix,fname], '-r600');
               print(h,'-depsc', [prefix,fname], '-r600');
    %            print(h,'-dpdf', [prefix,fname], '-r600');
           end
         if closefigs==1, close all,  end 
        


%% Part 2


%Pull out the data that we'll be using. Refer to spreadsheet for column names.   
clear indM; indM= IndAllSess(vertcat(AllSessMF{:})=='M');
clear indF; indF= IndAllSess(vertcat(AllSessMF{:})=='F');
% SleepHrs  = double(vertcat(BB2{:,12}));
%     SleepHrsM= double(vertcat(BB2{indM,12}));
%     SleepHrsF= double(vertcat(BB2{indF,12}));
% 
% SleepQual = double(vertcat(BB2{:,13}));
%     SleepQualM= double(vertcat(BB2{indM,13}));
%     SleepQualF= double(vertcat(BB2{indF,13}));

%Mood
clear Col; Col=50;    
KSS_Mood_pre    = double(vertcat(BB2{:,Col}));
    KSS_Mood_preM    = double(vertcat(BB2{indM,Col}));
    KSS_Mood_preF    = double(vertcat(BB2{indF,Col}));
clear Col; Col=56;      
KSS_Mood_post   = double(vertcat(BB2{:,Col}));
    KSS_Mood_postM   = double(vertcat(BB2{indM,Col}));
    KSS_Mood_postF   = double(vertcat(BB2{indF,Col}));

%Anxiety
clear Col; Col=51;    
KSS_Anx_pre    = double(vertcat(BB2{:,Col}));
    KSS_Anx_preM    = double(vertcat(BB2{indM,Col}));
    KSS_Anx_preF    = double(vertcat(BB2{indF,Col}));
clear Col; Col=57;     
KSS_Anx_post   = double(vertcat(BB2{:,Col}));
    KSS_Anx_postM   = double(vertcat(BB2{indM,Col}));
    KSS_Anx_postF   = double(vertcat(BB2{indF,Col}));

%Energy
clear Col; Col= 52;   
KSS_Enrgy_pre    = double(vertcat(BB2{:,Col}));
    KSS_Enrgy_preM    = double(vertcat(BB2{indM,Col}));
    KSS_Enrgy_preF    = double(vertcat(BB2{indF,Col}));
clear Col; Col=58;     
KSS_Enrgy_post   = double(vertcat(BB2{:,Col}));
    KSS_Enrgy_postM   = double(vertcat(BB2{indM,Col}));
    KSS_Enrgy_postF   = double(vertcat(BB2{indF,Col}));
    
%Pain
clear Col; Col= 49;   
KSS_Pain_pre    = double(vertcat(BB2{:,Col}));
    KSS_Pain_preM    = double(vertcat(BB2{indM,Col}));
    KSS_Pain_preF    = double(vertcat(BB2{indF,Col}));
clear Col; Col= 55;    
KSS_Pain_post   = double(vertcat(BB2{:,Col}));
    KSS_Pain_postM   = double(vertcat(BB2{indM,Col}));
    KSS_Pain_postF   = double(vertcat(BB2{indF,Col}));    


figure; 

    %Setting up location of axes in figure
%      hp1 = uipanel('position',[0.1331*1.05 0.1501-0.1 0.1535*1.45 0.7749+0.1]);   
%      hp2 = uipanel('position',[0.3361*1.05 0.1503-0.1 0.1566*1.45 0.7747+0.1]);
%      hp3 = uipanel('position',[0.5422*1.05 0.1503-0.1 0.1566*1.45 0.7747+0.1]);
%      hp4 = uipanel('position',[0.7484*1.05 0.1503-0.1 0.1566*1.45 0.7747+0.1]);

     
     hp1 = uipanel('position',[0.005706742627348,0.0501,0.222575000000004,0.8749]);   
     hp2 = uipanel('position',[0.241913042895449,0.0503,0.227070000000003,0.8747]);
     hp3 = uipanel('position',[0.489953431635396,0.0503,0.227070000000001,0.8747]);
     hp4 = uipanel('position',[0.737026434316362,0.0503,0.22707,0.8747]);

     Vars={'Mood','Anx', 'Enrgy','Pain'};
     VarsLabel={'Mood','Anxiety', 'Energy','Pain'};
     PlotThese={'KSS_Mood_pre ','KSS_Mood_pre','KSS_Mood_preM','KSS_Mood_preF','KSS_Mood_postM','KSS_Mood_postF',}
     
     
     for ii = 1:length(Vars)
         clear KSSpre KSSpost KSSpreM KSSpostM KSSpreF KSSpostF
        KSSpre=eval(['KSS_' Vars{ii} '_pre']);
            KSSpreM=eval(['KSS_' Vars{ii} '_preM']);
            KSSpreF=eval(['KSS_' Vars{ii} '_preF']);
        
        
        KSSpost=eval(['KSS_' Vars{ii} '_post']);
                KSSpostM=eval(['KSS_' Vars{ii} '_postM']);
                KSSpostF=eval(['KSS_' Vars{ii} '_postF']);

       scatterhist(KSSpre,KSSpost, 'Kernel','on','Marker','none','Direction','out','Location','NorthEast','Color',hisclr,'Parent',eval(['hp' num2str(ii)]));
        
        hold on 
        plot(KSSpreM, KSSpostM, '*', 'linewidth', 3,'Color', Mclr2);
        hold on
        plot(KSSpreF, KSSpostF, '*', 'linewidth', 3,'Color', Fclr2);
        hold on
                plot(nanmean(KSSpreM), nanmean(KSSpostM), '^', 'linewidth', 1,'Color', Mclr,'MarkerFaceColor',Mclr);
                hold on
                plot(nanmean(KSSpreF), nanmean(KSSpostF), '^', 'linewidth', 1,'Color', Fclr,'MarkerFaceColor',Fclr);
                hold on
        
        hold on
        line(nanmean(KSSpre) + [-nanstd(KSSpre) nanstd(KSSpre)], [1, 1].*nanmean(KSSpost),'Color', [0.5 0.5 0.5]);
        hold on
        line([1, 1].*nanmean(KSSpre),nanmean(KSSpost) + [-nanstd(KSSpost) nanstd(KSSpost)],'Color', [0.5 0.5 0.5]);
        plot(nanmedian(KSSpre), nanmedian(KSSpost), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
        plot(nanmean(KSSpre), nanmean(KSSpost), '^', 'linewidth', 1,'Color', [0 0 0],'MarkerFaceColor',[0 0 0]);
       

        xlabel([VarsLabel{ii} ' Pre']);
        ylabel([VarsLabel{ii} ' Post']) ;
        if contains(Vars{ii},'Mood')
        % %xlim([0.75 9.25]);
%         ylim([3.5 8.5]);
        end 
        set(gca, 'box','off','XAxisLocation','bottom','YAxisLocation','left');      

     end 
        fname='All Demo and PrePost Questionnaire-Part2';
        hh=gcf;
        set(hh,'name',fname,'Position', [ -2064         518        1865         466]);
        set(gca, 'box','off','XAxisLocation','bottom');
        
        if SveAllpics==1
               h = gcf;
               saveas(h,strcat(prefix,fname,'.fig'),'fig');
               saveas(h,strcat(prefix,fname,'.png'),'png');
               print(h,'-dpng', [prefix,fname], '-r600');
               print(h,'-depsc', [prefix,fname], '-r600');
    %            print(h,'-dpdf', [prefix,fname], '-r600');
           end
           if closefigs==1, close all,  end 
