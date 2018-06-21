
%%% This function takes the usual axona posfile and spikefile inputs (posx, cellTS etc), 
%%% computes the direction vector, and breaks the spiketrains into directional portions.
%%% It then iterates through the spiketrains and finds "bursts". These are
%%% defined as spikes preceeded by a 300ms quiescent period, accompanied by
%%% at least one other spike within 300ms. It returns: The interspike
%%% interval overall and for each direction, the percent of spikes fired in
%%% bursts, the Fano Factor and accompanying 95% gamma distribution
%%% confidence intervals and therefore whether the cell is more, less or
%%% the same variance as a Poisson process.

%%% Written by Robert Munn, Stanford University 2018

function [oldscore,oldscorex,oldscorey,oldslope,oldslopey,oldslopex,oldburstscore,oldburstscorex,oldburstscorey,oldslopeburst,oldslopeburstx,oldslopebursty,...
    pref_dir_circ,pref_dir_burst_circ,spkDir,rotation,speedScore_overall,slope_overall,speedScore_burst_overall,slope_burst_overall,speedScore_north,slope_north,speedScore_south,slope_south,speedScore_east,slope_east,speedScore_west,slope_west,...
    speedScore_burst_north,slope_burst_north,speedScore_burst_south,slope_burst_south,speedScore_burst_east,slope_burst_east,speedScore_burst_west,slope_burst_west,...
    speedScore,hdAxis,hd_fr,hdAxisburst,hd_frburst,mvl,mvlburst,mv_arg,mv_argburst,pref_angle,pref_angleburst,...
    num_overall_burst_spikes,num_south_burst_spikes,num_north_burst_spikes,num_west_burst_spikes,num_east_burst_spikes,...
    speed,hdDir,spiketrain_burst_overall,isi_overall,FF_overall,gamma_ci_overall,Poisson_overall,percent_overall_burst,...
    spiketrain_burst_north,isi_north,FF_north,gamma_ci_north,Poisson_north,percent_north_burst,...
    spiketrain_burst_south,isi_south,FF_south,gamma_ci_south,Poisson_south,percent_south_burst,...
    spiketrain_burst_east,isi_east,FF_east,gamma_ci_east,Poisson_east,percent_east_burst,...
    spiketrain_burst_west,isi_west,FF_west,gamma_ci_west,Poisson_west,percent_west_burst] = bursting_analysis(posx,posx2,posy,posy2,post,cellTS)

percent_west_burst = 0;
percent_east_burst = 0;
percent_north_burst = 0;
percent_south_burst = 0;
n = 1;

if length(posy2) > length(posy)
    posy2 = posy2(1:length(posy));
end
if length(posy) > length(posy2)
    posy = posy(1:length(posy2));
end
if length(posx2) > length(posx)
    posx2 = posx2(1:length(posx));
end
if length(posx) > length(posx2)
    posx = posx(1:length(posx2));
end

posy = -posy;
posy2 = -posy2; %the output from database maker is the mirror image of tint. This flips the output to be the same

minX = nanmin(posx); maxX = nanmax(posx);
minY = nanmin(posy); maxY = nanmax(posy);
xLength = maxX - minX;
yLength = maxY - minY;

posx = (posx - min(posx)); % Scale to boxsize (100cm)
posy = (posy - min(posy));

if xLength > yLength
    scalefac = 100/xLength;
posx = posx*scalefac;
posx2 = posx2*scalefac;
posy = posy*scalefac;
posy2 = posy2*scalefac;
elseif yLength > xLength
        scalefac = 100/yLength;
posx = posx*scalefac;
posx2 = posx2*scalefac;
posy = posy*scalefac;
posy2 = posy2*scalefac;
end
    
minX = nanmin(posx); maxX = nanmax(posx);
minY = nanmin(posy); maxY = nanmax(posy);
xLength = maxX - minX;
yLength = maxY - minY;


if yLength > xLength+10 % if the stable wall (cue wall) is west % rotate 90 degrees ccw
    posx_new = -posy; posx2_new = -posy2; posy = posx; posy2 = posx2;
    posx = posx_new; posx2 = posx2_new;
    rotation = 1;
else
    rotation = 0;
end
% Low speed threshold (cm/s)
lowSpeedThreshold = 1.5; 
% High speed threshold (cm/2)
highSpeedThreshold = 120;

if lowSpeedThreshold > 0 || highSpeedThreshold > 0
    % Calculate the speed at each timepoint
    speed = speed2D(posx,posy,post);
    speedx = speed1D(posx,post);
    speedy = speed1D(posy,post);

    bad_index = find(speed < lowSpeedThreshold | speed > highSpeedThreshold);
    % Remove the segments that have too high or too low speed
    posx(bad_index) = NaN;
    posy(bad_index) = NaN;
    posx2(bad_index) = NaN;
    posy2(bad_index) = NaN;
    speed(bad_index) = NaN;
    speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
    speed = gauss_smoothing(speed,10); % gaussian kernel smoothing
end
% align spikes %
roundCellTS = round(cellTS*100)/100;

hdDir = compute_hd(posx,posy);
spiketrain = computeFR(roundCellTS,post);
spiketrain(bad_index) = NaN; % take out spikes at bad movement indices

spiketrain_south = post; % initialise
spiketrain_north = post;
spiketrain_east = post;
spiketrain_west = post;
overall_id = find(spiketrain > 0);

if max(overall_id) > length(hdDir)
    diff = max(overall_id) - length(hdDir);
    overall_id = overall_id(1:end-diff);
end

speedScore = corr(speed,spiketrain,'rows','complete');
spkDir = hdDir(overall_id);
spkDir(isnan(spkDir)) = [];
spkDir_rads = deg2rad(spkDir);
%do circular analysis
[pref_dir_circ,~] = circ_vmpar(spkDir_rads);
pref_dir_circ = rad2deg(pref_dir_circ + pi);
hd = hdstat(spkDir,hdDir,0.02, 1);
hdAxis = hd.axis; hd_fr = hd.ratemap;
% compute mean vector length and argument
exp_dir = exp(-1i*hdAxis);
rayleigh_vector = pi/(numel(hdAxis)*sin(pi/numel(hdAxis)))*sum(hd_fr.*exp_dir)/sum(hd_fr);
mvl = sqrt(real(rayleigh_vector)^2+imag(rayleigh_vector)^2);
mv_arg = -atan2(imag(rayleigh_vector),real(rayleigh_vector))*180/pi;
if mv_arg < 0
    mv_arg = mv_arg+360; %make mv_arg give proper degree output
else
end
% compute peak firing rate and preferred angle
pref_angle = rad2deg(hd.mean+pi);
try
dir_south = hdDir > 150 & hdDir < 210;
spiketrain_south = spiketrain(dir_south);
speed_south = speed(dir_south);
dir_north = hdDir > 330 | hdDir < 30;
spiketrain_north = spiketrain(dir_north);
speed_north = speed(dir_north);
dir_east = hdDir < 120 & hdDir > 60;
spiketrain_east = spiketrain(dir_east);
speed_east = speed(dir_east);
dir_west = hdDir < 300 & hdDir > 240;
spiketrain_west = spiketrain(dir_west);
speed_west = speed(dir_west);
catch
    keyboard
end
north = 0; south = 0; east = 0; west = 0;
if sum(spiketrain_south) == 0
    disp('No south spikes?')
    south = 1;
end
if sum(spiketrain_north) == 0
    disp('No north spikes?')
   north = 1;
end
if sum(spiketrain_east) == 0
    disp('No east spikes?')
   east = 1;
end
if sum(spiketrain_west) == 0
    disp('No west spikes?')
    west = 1;
end

if north == 1 && south == 1 && east == 1 && west == 1
    keyboard
end

overall_id = find(spiketrain >=1);   % get the index of spikes and also the interspike interval (isi)
for y = 1:length(overall_id)-1
    isi_overall(y) = (overall_id(y+1) - overall_id(y))./20;
end
north_id = find(spiketrain_north >= 1);
for y = 1:length(north_id)-1
    isi_north(y) = (north_id(y+1) - north_id(y))./20;
end
south_id = find(spiketrain_south >= 1);
for y = 1:length(south_id)-1
    isi_south(y) = (south_id(y+1) - south_id(y))./20;
end
east_id = find(spiketrain_east >= 1);
for y = 1:length(east_id)-1
    isi_east(y) = (east_id(y+1) - east_id(y))./20;
end
west_id = find(spiketrain_west >= 1);
for y = 1:length(west_id)-1
    isi_west(y) = (west_id(y+1) - west_id(y))./20;
end

north_more_than_1(n,1) = sum(spiketrain_north > 1); % find the number of fast bursts (more than 1 spike/20ms bin)
south_more_than_1(n,1) = sum(spiketrain_south > 1);
west_more_than_1(n,1) = sum(spiketrain_west > 1);
east_more_than_1(n,1) = sum(spiketrain_east > 1);


[overall_burst,overall_num_spikes_300ms,spiketrain_burst_overall] = compute_burst(overall_id,spiketrain);
num_overall_spikes_300ms_episode{n} = overall_num_spikes_300ms;
num_overall_burst_spikes(n,1) = overall_burst;
percent_overall_burst(n,1) = overall_burst/nansum(spiketrain)*100;
overall_burst_index = find(spiketrain_burst_overall > 0);
spkDirburst = hdDir(overall_burst_index);
spkDirburst(isnan(spkDirburst)) = [];
spkDirburst_rads = deg2rad(spkDirburst);
[pref_dir_burst_circ conc_param_burst_burst] = circ_vmpar(spkDirburst_rads);
pref_dir_burst_circ = rad2deg(pref_dir_burst_circ + pi);
hdburst = hdstat(spkDirburst,hdDir,0.02, 1);
hdAxisburst = hdburst.axis; hd_frburst = hdburst.ratemap;
% compute mean vector length and argument
exp_dirburst = exp(-1i*hdAxisburst);
rayleigh_vectorburst = pi/(numel(hdAxisburst)*sin(pi/numel(hdAxisburst)))*sum(hd_frburst.*exp_dirburst)/sum(hd_frburst);
mvlburst = sqrt(real(rayleigh_vectorburst)^2+imag(rayleigh_vectorburst)^2);
mv_argburst = -atan2(imag(rayleigh_vectorburst),real(rayleigh_vectorburst))*180/pi;
if mv_argburst < 0
    mv_argburst = mv_argburst+360; %make mv_arg give proper degree output
else
end
% compute peak firing rate and preferred angle
pref_angleburst = rad2deg(hdburst.mean + pi);
[west_burst,west_num_spikes_300ms,spiketrain_burst_west] = compute_burst(west_id,spiketrain_west);
num_west_spikes_300ms_episode{n} = west_num_spikes_300ms;
num_west_burst_spikes(n,1) = west_burst;
percent_west_burst(n,1) = west_burst/nansum(spiketrain_west)*100;

[south_burst,south_num_spikes_300ms,spiketrain_burst_south] = compute_burst(south_id,spiketrain_south);
num_south_spikes_300ms_episode{n} = south_num_spikes_300ms;
num_south_burst_spikes(n,1) = south_burst;
percent_south_burst(n,1) = south_burst/nansum(spiketrain_south)*100;

[north_burst,north_num_spikes_300ms,spiketrain_burst_north] = compute_burst(north_id,spiketrain_north);
num_north_spikes_300ms_episode{n} = north_num_spikes_300ms;
num_north_burst_spikes(n,1) = north_burst;
percent_north_burst(n,1) = north_burst/nansum(spiketrain_north)*100;

[east_burst,east_num_spikes_300ms,spiketrain_burst_east] = compute_burst(east_id,spiketrain_east);
num_east_spikes_300ms_episode{n} = east_num_spikes_300ms;
num_east_burst_spikes(n,1) = east_burst;
percent_east_burst(n,1) = east_burst/nansum(spiketrain_east)*100;

spiketrain(isnan(spiketrain)) = 0;
spiketrain_burst_overall(isnan(spiketrain_burst_overall)) = 0;

fr = gauss_smoothing(spiketrain,20)*50;
fr_north = gauss_smoothing(spiketrain_north,20)*50;
fr_south = gauss_smoothing(spiketrain_south,20)*50;
fr_east = gauss_smoothing(spiketrain_east,20)*50;
fr_west = gauss_smoothing(spiketrain_west,20)*50;

speedScore_overall = corr(fr,speed,'rows','complete');
speedScore_north = corr(fr_north,speed_north,'rows','complete'); % do speedscore etc for all spikes
speedScore_south = corr(fr_south,speed_south,'rows','complete');
speedScore_east = corr(fr_east,speed_east,'rows','complete');
speedScore_west = corr(fr_west,speed_west,'rows','complete');

poly_ov = polyfit(fr,speed,1); slope_overall = poly_ov(1);
poly_n = polyfit(fr_north,speed_north,1); slope_north = poly_n(1);
poly_s = polyfit(fr_south,speed_south,1); slope_south = poly_s(1);
poly_e = polyfit(fr_east,speed_east,1); slope_east = poly_e(1);
poly_w = polyfit(fr_west,speed_west,1); slope_west = poly_w(1);

fr_burst_overall = gauss_smoothing(spiketrain_burst_overall,20)*50;
fr_burst_north = gauss_smoothing(spiketrain_burst_north,20)*50;
fr_burst_south = gauss_smoothing(spiketrain_burst_south,20)*50;
fr_burst_east = gauss_smoothing(spiketrain_burst_east,20)*50;
fr_burst_west = gauss_smoothing(spiketrain_burst_west,20)*50;

speedScore_burst_overall = corr(fr_burst_overall',speed(1:length(fr_burst_overall)),'rows','complete');
speedScore_burst_north = corr(fr_burst_north',speed_north(1:length(fr_burst_north)),'rows','complete'); % do speedscore etc for burst spikes
speedScore_burst_south = corr(fr_burst_south',speed_south(1:length(fr_burst_south)),'rows','complete');
speedScore_burst_east = corr(fr_burst_east',speed_east(1:length(fr_burst_east)),'rows','complete');
speedScore_burst_west = corr(fr_burst_west',speed_west(1:length(fr_burst_west)),'rows','complete');

poly_bov = polyfit(fr_burst_overall',speed(1:length(fr_burst_overall)),1); slope_burst_overall = poly_bov(1);
poly_bn = polyfit(fr_burst_north',speed_north(1:length(fr_burst_north)),1); slope_burst_north = poly_bn(1);
poly_bs = polyfit(fr_burst_south',speed_south(1:length(fr_burst_south)),1); slope_burst_south = poly_bs(1);
poly_be = polyfit(fr_burst_east',speed_east(1:length(fr_burst_east)),1); slope_burst_east = poly_be(1);
poly_bw = polyfit(fr_burst_west',speed_west(1:length(fr_burst_west)),1); slope_burst_west = poly_bw(1);

FF_overall(n,1) = nanvar(spiketrain)/nanmean(spiketrain); % Compute the Fano Factors (Poisson procress generates a FF of 1) and 95% CI of the respective gamma distributions
gamma_ci_overall(n,1:2) = gaminv([.025,.975],(size(spiketrain,1) - 1)/2,2/(size(spiketrain,1) - 1));

 if FF_overall(n,1) > gamma_ci_overall(n,2)
     Poisson_overall{n,1} = 'More Variable';
 elseif FF_overall(n,1) < gamma_ci_overall(n,1)
     Poisson_overall{n,1} = 'Less Variable';
 else
     Poisson_overall{n,1} = 'Poisson';
 end
     
FF_east(n,1) = nanvar(spiketrain_east)/nanmean(spiketrain_east);
gamma_ci_east(n,1:2) = gaminv([.025,.975],(size(spiketrain_east,1) - 1)/2,2/(size(spiketrain_east,1) - 1));

 if FF_east(n,1) > gamma_ci_east(n,2)
     Poisson_east{n,1} = 'More Variable';
 elseif FF_east(n,1) < gamma_ci_east(n,1)
     Poisson_east{n,1} = 'Less Variable';
 else
     Poisson_east{n,1} = 'Poisson';
 end

FF_west(n,1) = nanvar(spiketrain_west)/nanmean(spiketrain_west);
gamma_ci_west(n,1:2) = gaminv([.025,.975],(size(spiketrain_west,1) - 1)/2,2/(size(spiketrain_west,1) - 1));

 if FF_west(n,1) > gamma_ci_west(n,2)
     Poisson_west{n,1} = 'More Variable';
 elseif FF_west(n,1) < gamma_ci_west(n,1)
     Poisson_west{n,1} = 'Less Variable';
 else
     Poisson_west{n,1} = 'Poisson';
 end
 
FF_north(n,1) = nanvar(spiketrain_north)/nanmean(spiketrain_north);
gamma_ci_north(n,1:2) = gaminv([.025,.975],(size(spiketrain_north,1) - 1)/2,2/(size(spiketrain_north,1) - 1));

 if FF_north(n,1) > gamma_ci_north(n,2)
     Poisson_north{n,1} = 'More Variable';
 elseif FF_north(n,1) < gamma_ci_north(n,1)
     Poisson_north{n,1} = 'Less Variable';
 else
     Poisson_north{n,1} = 'Poisson';
 end

FF_south(n,1) = nanvar(spiketrain_south)/nanmean(spiketrain_south);
gamma_ci_south(n,1:2) = gaminv([.025,.975],(size(spiketrain_south,1) - 1)/2,2/(size(spiketrain_south,1) - 1));

 if FF_south(n,1) > gamma_ci_south(n,2)
     Poisson_south{n,1} = 'More Variable';
 elseif FF_south(n,1) < gamma_ci_south(n,1)
     Poisson_south{n,1} = 'Less Variable';
 else
     Poisson_south{n,1} = 'Poisson';
 end

auto_overall(1:21,n) = xcorr(spiketrain-nanmean(spiketrain),10,'coeff'); % Calculate the autocorrelation with a 100ms window (bins are 20ms)
CI_overall(1:2,n) = nanmean(spiketrain) + tinv([0.025  0.975],length(spiketrain)-1)*nanstd(spiketrain)/sqrt(length(spiketrain));   

auto_south(1:21,n) = xcorr(spiketrain_south-nanmean(spiketrain_south),10,'coeff');
CI_south(1:2,n) = nanmean(spiketrain_south) + tinv([0.025  0.975],length(spiketrain_south)-1)*nanstd(spiketrain_south)/sqrt(length(spiketrain_south));   

auto_north(1:21,n) = xcorr(spiketrain_north-nanmean(spiketrain_north),10,'coeff');
CI_north(1:2,n) = nanmean(spiketrain_north) + tinv([0.025  0.975],length(spiketrain_north)-1)*nanstd(spiketrain_north)/sqrt(length(spiketrain_north));   

auto_east(1:21,n) = xcorr(spiketrain_east-nanmean(spiketrain_east),10,'coeff');
CI_east(1:2,n) = nanmean(spiketrain_east) + tinv([0.025  0.975],length(spiketrain_east)-1)*nanstd(spiketrain_east)/sqrt(length(spiketrain_east));   

auto_west(1:21,n) = xcorr(spiketrain_west-nanmean(spiketrain_west),10,'coeff');
CI_west(1:2,n) = nanmean(spiketrain_west) + tinv([0.025  0.975],length(spiketrain_west)-1)*nanstd(spiketrain_west)/sqrt(length(spiketrain_west));   

h = hist(cellTS,post);  % compute x and y score and slopes using old TC method considering only posx or posy as opposed to broken spiketrains
h = gauss_smoothing(h,20)*50;
burstTS(:,1) = find(spiketrain_burst_overall >= 1);
burstTS = burstTS*0.02;
hb = hist(burstTS,post);
hb = gauss_smoothing(hb,20)*50;

if isempty(burstTS)
oldscore = NaN;
oldburstscore = NaN; 
oldscorex = NaN;
oldburstscorex = NaN; 
oldscorey = NaN;
oldburstscorey = NaN; 
oldslope = NaN;
oldslopeburst = NaN; 
oldslopey = NaN;
oldslopebursty = NaN; 
oldslopex = NaN;
oldslopeburstx = NaN;

return

end
speedx = inpaint_nans(speedx,1); % interpolate nans in directional speed array
speedy = inpaint_nans(speedy,1);
try
oldscore = corr(h',speed);
oldburstscore = corr(hb',speed);
oldscorex = corr(h',speedx,'rows','complete');
oldburstscorex = corr(hb',speedx,'rows','complete');
oldscorey = corr(h',speedy,'rows','complete');
oldburstscorey = corr(hb',speedy,'rows','complete');
oldslope = polyfit(h',speed,1); oldslope = oldslope(1);
oldslopeburst = polyfit(hb',speed,1); oldslopeburst = oldslopeburst(1);
oldslopey = polyfit(h',speedy,1); oldslopey = oldslopey(1);
oldslopebursty = polyfit(hb',speedy,1); oldslopebursty = oldslopebursty(1);
oldslopex = polyfit(h',speedx,1); oldslopex = oldslopex(1);
oldslopeburstx = polyfit(hb',speedx,1); oldslopeburstx = oldslopeburstx(1);
catch
    keyboard
end


end
