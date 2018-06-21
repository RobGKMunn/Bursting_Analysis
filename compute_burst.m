%%% takes an index of spikes (the bin number where they occur), and a train
%%% of binned spikecounts and returns the number of spikes in bursts
%%% (burst), the number of spikes in each 300ms increment after each spike
%%% (num_spikes_300ms) and a spiketrain containing only the "burst" spikes
%%% (spiketrainburst). Assumes a 20ms bin interval, and therefore defines a
%%% "burst" as a spike that is preceeded by a 300ms (20x15) quiescent
%%% period, with at least one other spike in the following 300ms. 

%%% Written By Robert Munn, Stanford University, 2018


function [burst,num_spikes_300ms,spiketrainburst] = compute_burst(spike_index,spiketrain)

burst = 0;
num = 0;
for g = 1:length(spike_index)
for u = spike_index(g)  
    if u-15 <= 0
        if sum(spiketrain(u+1:u+15)) >= 1   
           burst = burst + sum(spiketrain(u+1:u+15))+1;
           spiketrain_burst(u:u+15) = spiketrain(u:u+15);
           num = sum(spiketrain(u:u+15));
        else
           num = sum(spiketrain(u:u+15));
           spiketrainburst(u:u+15) = 0;
        end
        continue
    end
   if u+15 >= length(spiketrain)
                if sum(spiketrain(u+1:length(spiketrain))) >= 1
                    burst = burst + sum(spiketrain(u+1:length(spiketrain)))+1;
                    spiketrainburst(u:length(spiketrain)) = spiketrain(u:length(spiketrain));
                    num = sum(spiketrain(u:length(spiketrain)));
                else
                     num = sum(spiketrain(u:length(spiketrain)));
                    spiketrainburst(u:length(spiketrain)) = 0;
                end
    continue
   end
   
   if sum(spiketrain(u-15:u-1)) == 0
            if sum(spiketrain(u+1:u+15)) >= 1   
            burst = burst + sum(spiketrain(u+1:u+15))+1;
            spiketrainburst(u:u+15) = spiketrain(u:u+15);
            num = sum(spiketrain(u:u+15));
            else
                num = sum(spiketrain(u:u+15));
                spiketrainburst(u:u+15) = 0;
            end
   end    
        num_spikes_300ms(g,1) = num;
end
end
end