




%% Load the data
Data_phase = {'Phase_1_EEG_SpecInfo.mat';'Phase_2_EEG_SpecInfo.mat'};
Monts_Phase1 ={'F0','F5','F30','M0','M5','M30','P0','P5','P30'};
Monts_Phase2 ={'F30','M30'};

PeakPSD =[];
PeakFreq =[]; 
PeakPSD_DurinPre_Ratio =[];
for load_phase=2%:2





load(['D:\GX\Results\DataChunkedtoTrials_AllPhases_06142021\FigOutput\', Data_phase{load_phase}])


if load_phase==1
    monts_phase_match =Monts_Phase1;
else
    monts_phase_match =Monts_Phase2;
end


for ind=1:size(EEG_SpecInfo,1)

    subj =str2num(EEG_SpecInfo{ind,2}(1:2));
    
    
    
    
dummy_mont =vertcat(EEG_SpecInfo{ind,1}(:,2))    
dummy_mont_uniq =unique(dummy_mont)

for mm=1:length(dummy_mont_uniq)
    dummy_mont_uniq_ind = strmatch(dummy_mont_uniq{mm},dummy_mont)
    mont_ind=strmatch(dummy_mont_uniq{mm},monts_phase_match);
    
    
    for kk=1:length(dummy_mont_uniq_ind)
       PeakPSD(mont_ind,kk,subj)= EEG_SpecInfo{ind,1}{dummy_mont_uniq_ind(kk),4}
       PeakFreq(mont_ind,kk,subj)= EEG_SpecInfo{ind,1}{dummy_mont_uniq_ind(kk),5}
       PeakPSD_DurinPre_Ratio(mont_ind,kk,subj)= EEG_SpecInfo{ind,1}{dummy_mont_uniq_ind(kk),1}
    end 
    
end

end 


end 
