# Bursting_Analysis
Code used to find bursts in spiketrains and then produce behavioral correlates for the full train vs. just the bursts. Bursts are defined as any spike preceded by at least 300ms quiescence and followed by at least one other spike within 300ms

Relies on inpaint_nans.m from the MATLAB file exchange to interpolate values in speed vectors
