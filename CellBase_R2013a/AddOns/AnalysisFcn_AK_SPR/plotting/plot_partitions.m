function plot_partitions(part1,colors1,x,width)
%
% plot_partitions(part1,colors1,x,width)
%
% colors1 = {'c','b','g','r','m','k'};
%
numcolors = length(colors1);

hold on;
for i = 1:size(part1,1)
    colnum = mod(i-1,numcolors)+1;
    plot([x x],[part1(i,1) part1(i,2)],'LineWidth',width,'Color',colors1{colnum});
end