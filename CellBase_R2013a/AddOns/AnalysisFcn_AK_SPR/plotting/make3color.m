function map=make3color(c1,c2,c3,varargin)


cd1=[(c2(1)-c1(1))/31 (c2(2)-c1(2))/31 (c2(3)-c1(3))/31];
cd2=[(c3(1)-c2(1))/31 (c3(2)-c2(2))/31 (c3(3)-c2(3))/31];

if nargin > 3
    mapsize = varargin{1};
else
    mapsize = 64;
end

for (i=0:30)
   map(i+1,:)= c1 + cd1*i;
end
for (i=1:31)
   map(i+32,:)=c2+cd2*i;
end


map(1,:)  = c1;
map(32,:) = c2;
map(64,:) = c3;
