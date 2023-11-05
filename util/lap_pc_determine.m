function [pc,all_infoScore,all_infoScore_norm,all_coherence]=lap_pc_determine(neuron_laps,behav_laps,binsize)

%young
pc={};
all_infoScore={};
all_infoScore_norm={};
all_coherence={};

for i=1:size(neuron_laps,1)
    for j=1:size(neuron_laps,2)
        % merge laps
        n=[];
        b=[];


        n.C=[];
        n.S=[];
        n.time=[];
        b.pos=[];
        b.time=[];

        for k=1:length(neuron_laps{i,j})
            n.C=[n.C,neuron_laps{i,j}{k}.C];
            n.S=[n.S,neuron_laps{i,j}{k}.S];
            n.time=[n.time;neuron_laps{i,j}{k}.time];
            b.pos=[b.pos;behav_laps{i,j}{k}.position];
            b.time=[b.time;behav_laps{i,j}{k}.time];
        end

        [pc{i,j},all_infoScore{i,j},all_infoScore_norm{i,j},all_coherence{i,j}] = permutingSpike_adapt_032623(n,b.pos,b.time,'S',0,binsize,10,'all',0);  

%         [pc{i,j}] = permutingPeak_032623(n,b.pos,b.time,'S',25,10);
    end
end