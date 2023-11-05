function b=behav_split(behav,frames,interval)

startt=1;

b={};
for i=1:floor((size(behav.position,1)-frames)/interval)+1
    b{i}.position=behav.position(startt:startt+frames-1,:);
    b{i}.time=behav.time(startt:startt+frames-1);

    startt=startt+interval;
end
