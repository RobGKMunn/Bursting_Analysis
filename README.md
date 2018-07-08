# Bursting_Analysis
Code used to find bursts in spiketrains and then produce behavioral correlates for the full train vs. just the bursts. Bursts are defined 

The master function do_bursting_analysis.m takes input in the form of a spreadsheet with columns "Sessions" "Tetrode" and "Unit". Sessions should contain a link to the databasemaker .mat file output. The script loads the appropriate .pos and .spk files, and rotates the position so that all the environments are aligned the same way. Head direction and speed are calculated. This data is used to break the spiketrain up into segments when the animal is facing each of the quadrants (NSEW).

The interspike interval is calculated for each spiketrain, and bursts are identified. Bursts are defined as a spike that has been preceeded by at least 300ms of quiescence, and is followed by at least one more spike in the subsequent 300ms. The script then calculates the number and percentage of spikes fired as bursts (vs single spikes). Spiketrains are generated from these indices that include only the "burst" spikes.

Behavior decoding is done on both the complete and burst spiketrains (i.e Mean directional vector length and direction, correlation between running speed and firing rate etc).

Relies on inpaint_nans.m from the MATLAB file exchange to interpolate values in speed vectors
