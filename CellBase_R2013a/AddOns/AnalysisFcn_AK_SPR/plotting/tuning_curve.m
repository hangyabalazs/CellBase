function [stim, F, Fse] =  tuning_curve(Freq,OdorValveID);
%
%
%

stim = unique(OdorValveID);

for i=1:length(stim)
   ftemp = Freq(OdorValveID == stim(i)); 
   F(i)   = mean(ftemp);
   Fsd(i) = std(ftemp);
   Fse(i) = std(ftemp)/sqrt(max(1,length(ftemp)-1));
end
