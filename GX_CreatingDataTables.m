

clear PercenDiffPerfPulled_shaped  dummy  
 for trl=1:4 %Sort by trials
for ss=1:10 %Subjects   
    if ss==1
  dummy(:,1)=PercenDiffPerfPulledInMean(ss,:,trl)';
    else 
  dummy=[dummy;PercenDiffPerfPulledInMean(ss,:,trl)'];
    end
end 

PercenDiffPerfPulled_shaped(:,trl)=dummy;
clear dummy
 end      
 
 
 %Put together Non UI table
 for ii=1:4
 clear dummy
 dummy = PercenDiffPerfPulled_shaped(:,ii);
 dummy =dummy(:)
 mask =cellfun(@(C) all(isnan(C)),dummy)
 dummy(mask)= {' '};
 
     switch(ii)
         case 1
         Trial_1=dummy;
         case 2
         Trial_2=dummy;
         case 3
         Trial_3=dummy;
         case 4
         Trial_4=dummy;
     end 
 end 

 
 
 %This adds the average deviation across trials
  clear dummy
 dummy =PerDiffPerfPulledInMeanPooled';
 dummy = num2cell(round(-1.*vertcat(dummy{:}),1));
  mask =cellfun(@(C) all(isnan(C)),dummy)
 dummy(mask)= {' '};
 Mean_Deviation =dummy;
PercenDiffPerfPulled_shaped = [PercenDiffPerfPulled_shaped,dummy];

 %This adds the standard deviation deviation across trials
 clear dummy
 dummy =PerDiffPerfPulledInStdPooled';
 dummy =dummy(:)
 mask =cellfun(@(C) all(isnan(C)),dummy)
 dummy(mask)= {' '};
Stand_dev =round(dummy,1);
 PercenDiffPerfPulled_shaped = [PercenDiffPerfPulled_shaped,dummy];


%This adds the stimulation type
Stim_type =repmat(MontageMat2',10,1)
PercenDiffPerfPulled_shaped = [Stim_type,PercenDiffPerfPulled_shaped];

%This adds the subject numbers
clear dummy
dummy =[NumUniqueSubjsNums.*ones(1,9)]';
dummy = num2cell(dummy(:))
Subj_N =dummy;
PercenDiffPerfPulled_shaped = [dummy,PercenDiffPerfPulled_shaped];



   data =PercenDiffPerfPulled_shaped;%squeeze(reshape(PercenDiffPerfPulledInMean,90,4));%PercenDiffPerfPulledInMean(:,:,1);
   columnName =   {'Subject Number','Stim. Type', 'Trial 1','Trial 2','Trial 3','Trial 4','Mean Deviation','Std.'};
%    fh = figure;
uifigure
   t = uitable('Units','normalized','Position',...
      [0 0 1 1], 'Data', data,... 
      'ColumnName', columnName,...
      'RowName',[]);
     t = uitable(uifigure, 'Data', data,... 
      'ColumnName', columnName,...
      'RowName',[]);
%    figPos = get(fh,'Position');
%    tableExtent = get(t,'Extent');
%    set(fh,'Position',[figPos(1:2), figPos(3:4).*tableExtent(3:4)]);
 
s1 = uistyle;
s1.BackgroundColor = [jet(3)];
addStyle(t,s1,'cell',[2 4;2 5;3 4])








 tt =table(Subj_N,Trial_1,Stim_type, Trial_2,Trial_3,Trial_4,Mean_Deviation,Stand_dev)
 uitable(uifigure, 'Data',tt,'ColumnName', columnName)

 %TO DO : Add the bandpower measures to the table once calculated. 
 %Can use xlswrite to export table
