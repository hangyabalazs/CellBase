cellid=pv_same_tetrode(4)  ;
tseg=findSegs2(cellid);
[stimes,seltsind,selisi]=extractSegSpikes(cellid,tseg);



stimes = TS*1e3;

isi = stimes(2:end)-stimes(1:end-1);
bins200=logspace(0,4,200);
% bins40=logspace(0,4,40);
bins=linspace(0,100,100);
count=histc(isi,bins200);
%  count=histc(isi,bins40);
figure(1)
clf;
stairs(bins200,count);
ba0=gca;
set(ba0,'XScale','log');

