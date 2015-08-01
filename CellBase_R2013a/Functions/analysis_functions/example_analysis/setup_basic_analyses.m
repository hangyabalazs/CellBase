addanalysis(@LRatio2,'property_names',{'ID_PC','Lr_PC'},'arglist',{'feature_names' {'WavePC1' 'Energy'}})
addanalysis(@LRatio2,'property_names',{'ID_amp','Lr_amp'},'arglist',{'feature_names' {'Amplitude' 'Energy'}})

for k=1:length(CELLIDLIST)
    cellid=CELLIDLIST{k};
    try
        ST = loadcb(cellid,'STIMSPIKES');
        if isequal(findcellstr(ST.events(:,1),'PulseOn'),0)
            prealignSpikes(CELLIDLIST(k),'FUNdefineEventsEpochs',...
                @defineEventsEpochs_pulseon,'filetype','stim',...
                'ifsave',1,'ifappend',1)
        end
    end
end

addanalysis(@nbisstim,'property_names',{'Hindex','D_KL'})
addanalysis(@spikeshapecorr,'property_names',{'R'})