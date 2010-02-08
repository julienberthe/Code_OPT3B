i=1;
for x=-10:0.1:10
[pds,dpds]=poids(x,60,0.1,'spline quadratique');
xxx(i)=pds;
xxxx(i)=dpds;
i=i+1;
end

plot(-10:0.1:10,xxx)
figure;
plot(-10:0.1:10,xxxx)
