% cross time trajectory

% circle: 
time_range={[1:10800],[5400:16200],[10800:21600]};
figure;
for tk=1:3
    subplot(3,1,tk)
    plot(behavIndividuals{2}.position(time_range{tk}(1):time_range{tk}(end),1),behavIndividuals{2}.position(time_range{tk}(1):time_range{tk}(end),2))
end

% square: 
time_range={[1:10800],[5400:16200],[10800:21600]};
figure;
for tk=1:3
    subplot(3,1,tk)
    plot(behavIndividuals{5}.position(time_range{tk}(1):time_range{tk}(end),1),behavIndividuals{5}.position(time_range{tk}(1):time_range{tk}(end),2))
end