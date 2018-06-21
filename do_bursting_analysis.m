clear all
tab = uigetfile('*.xlsx','FEED ME A SPREADSHEET'); % Read spreadsheet that contains path to databasemaker files (*.pos,*.eeg,*.TXCX where X is tetrode (T) and Unit (C) number). Spreadsheet has columns Session (paths), Tetrode, and Unit.
A = readtable(tab);
numCell = numel(A.Tetrode);
spreadsheet = tab;
Burst_Results = table;
%% LOAD DATA %%
for n = 1:numCell % make arrays of the locations of all the relevant position (*_pos.mat), spike (*_TXCX.mat), and eeg (*_eeg.mat) files
posfile = strcat(A.Sessions{n},'_pos.mat');
spikefile = strcat(A.Sessions{n},'_T',num2str(A.Tetrode(n)),'C',num2str(A.Unit(n)),'.mat');
load(posfile)
load(spikefile)

if sum(cellTS) == 0
    continue
end

tic
[oldscore,oldscorex,oldscorey,oldslope,oldslopey,oldslopex,oldburstscore,oldburstscorex,oldburstscorey,oldslopeburst,oldslopeburstx,oldslopebursty,...
    pref_dir_circ,pref_dir_burst_circ,spike_directions,rotation,speedScore_overall,slope_overall,speedScore_burst_overall,slope_burst_overall,speedScore_north,slope_north,speedScore_south,slope_south,speedScore_east,slope_east,speedScore_west,slope_west,...
    speedScore_burst_north,slope_burst_north,speedScore_burst_south,slope_burst_south,speedScore_burst_east,slope_burst_east,speedScore_burst_west,slope_burst_west,...
    speedScore,hdAxis,hd_fr,hdAxisburst,hd_frburst,mvl,mvlburst,mv_arg,mv_argburst,pref_angle,pref_angleburst,...
    num_overall_burst_spikes,num_south_burst_spikes,num_north_burst_spikes,num_west_burst_spikes,num_east_burst_spikes,...
    speed,hdDir,spiketrain_burst_overall,isi_overall,FF_overall,gamma_ci_overall,Poisson_overall,percent_overall_burst,...
    spiketrain_burst_north,isi_north,FF_north,gamma_ci_north,Poisson_north,percent_north_burst,...
    spiketrain_burst_south,isi_south,FF_south,gamma_ci_south,Poisson_south,percent_south_burst,...
    spiketrain_burst_east,isi_east,FF_east,gamma_ci_east,Poisson_east,percent_east_burst,...
    spiketrain_burst_west,isi_west,FF_west,gamma_ci_west,Poisson_west,percent_west_burst] = bursting_analysis(posx,posx2,posy,posy2,post,cellTS);
toc

if rotation == 1
    disp('Open Field Rotation')
elseif rotation == 0
    disp('No Open Field Rotation')
end
clear posx posx2 posy posy2 post cellTS

posfile_sq = strcat(A.Sessions_sq{n},'_pos.mat');
spikefile_sq = strcat(A.Sessions_sq{n},'_T',num2str(A.Tetrode_sq(n)),'C',num2str(A.Unit_sq(n)),'.mat');
load(posfile_sq);
load(spikefile_sq);
posx_sq = posx; posx2_sq = posx2; posy_sq = posy; posy2_sq = posy2; post_sq = post; cellTS_sq = cellTS;

if sum(cellTS_sq) == 0
    continue
end

tic
[oldscore_sq,oldscorex_sq,oldscorey_sq,oldslope_sq,oldslopey_sq,oldslopex_sq,oldburstscore_sq,oldburstscorex_sq,oldburstscorey_sq,oldslopeburst_sq,oldslopeburstx_sq,oldslopebursty_sq,...
    pref_dir_circ_sq,pref_dir_burst_circ_sq,spike_directions_sq,rotation_sq,speedScore_overall_sq,slope_overall_sq,speedScore_burst_overall_sq,slope_burst_overall_sq,speedScore_north_sq,slope_north_sq,speedScore_south_sq,slope_south_sq,speedScore_east_sq,slope_east_sq,speedScore_west_sq,slope_west_sq,...
    speedScore_burst_north_sq,slope_burst_north_sq,speedScore_burst_south_sq,slope_burst_south_sq,speedScore_burst_east_sq,slope_burst_east_sq,speedScore_burst_west_sq,slope_burst_west_sq,...
    speedScore_sq,hdAxis_sq,hd_fr_sq,hdAxisburst_sq,hd_frburst_sq,mvl_sq,mvlburst_sq,mv_arg_sq,mv_argburst_sq,pref_angle_sq,pref_angleburst_sq,...
    num_overall_burst_spikes_sq,num_south_burst_spikes_sq,num_north_burst_spikes_sq,num_west_burst_spikes_sq,num_east_burst_spikes_sq,...
    speed_sq,hdDir_sq,spiketrain_burst_overall_sq,isi_overall_sq,FF_overall_sq,gamma_ci_overall_sq,Poisson_overall_sq,percent_overall_burst_sq,...
    spiketrain_burst_north_sq,isi_north_sq,FF_north_sq,gamma_ci_north_sq,Poisson_north_sq,percent_north_burst_sq,...
    spiketrain_burst_south_sq,isi_south_sq,FF_south_sq,gamma_ci_south_sq,Poisson_south_sq,percent_south_burst_sq,...
    spiketrain_burst_east_sq,isi_east_sq,FF_east_sq,gamma_ci_east_sq,Poisson_east_sq,percent_east_burst_sq,...
    spiketrain_burst_west_sq,isi_west_sq,FF_west_sq,gamma_ci_west_sq,Poisson_west_sq,percent_west_burst_sq] = bursting_analysis(posx_sq,posx2_sq,posy_sq,posy2_sq,post_sq,cellTS_sq);
toc

if rotation_sq == 1
    disp("Rotation")
elseif rotation_sq == 0
    disp("No Rotation")
else
    disp("Something done stuffed up somewhere")
end

Burst_Results.Spikefile{n,1} = spikefile;
Burst_Results.Percent_Burst_Overall(n,1) = percent_overall_burst;
Burst_Results.num_overall_burst_spikes(n,1) = num_overall_burst_spikes;
Burst_Results.spiketrain_burst_overall{n,1} = spiketrain_burst_overall;
Burst_Results.isi_overall{n,1} = isi_overall;
Burst_Results.FF_overall(n,1) = FF_overall;
Burst_Results.gamma_ci_overall{n,1} = gamma_ci_overall;
Burst_Results.Poisson_overall{n,1} = Poisson_overall;


Burst_Results.Percent_Burst_north(n,1) = percent_north_burst;
Burst_Results.num_north_burst_spikes(n,1) = num_north_burst_spikes;
Burst_Results.spiketrain_burst_north{n,1} = spiketrain_burst_north;
Burst_Results.isi_north{n,1} = isi_north;
Burst_Results.FF_north(n,1) = FF_north;
Burst_Results.gamma_ci_north{n,1} = gamma_ci_north;
Burst_Results.Poisson_north{n,1} = Poisson_north;

Burst_Results.Percent_Burst_south(n,1) = percent_south_burst;
Burst_Results.num_south_burst_spikes(n,1) = num_south_burst_spikes;
Burst_Results.spiketrain_burst_south{n,1} = spiketrain_burst_south;
Burst_Results.isi_south{n,1} = isi_south;
Burst_Results.FF_south(n,1) = FF_south;
Burst_Results.gamma_ci_south{n,1} = gamma_ci_south;
Burst_Results.Poisson_south{n,1} = Poisson_south;

Burst_Results.Percent_Burst_east(n,1) = percent_east_burst;
Burst_Results.spiketrain_burst_east{n,1} = spiketrain_burst_east;
Burst_Results.num_east_burst_spikes(n,1) = num_east_burst_spikes;
Burst_Results.isi_east{n,1} = isi_east;
Burst_Results.FF_east(n,1) = FF_east;
Burst_Results.gamma_ci_east{n,1} = gamma_ci_east;
Burst_Results.Poisson_east{n,1} = Poisson_east;

Burst_Results.Percent_Burst_west(n,1) = percent_west_burst;
Burst_Results.num_west_burst_spikes(n,1) = num_west_burst_spikes;
Burst_Results.spiketrain_burst_west{n,1} = spiketrain_burst_west;
Burst_Results.isi_west{n,1} = isi_west;
Burst_Results.FF_west(n,1) = FF_west;
Burst_Results.gamma_ci_west{n,1} = gamma_ci_west;
Burst_Results.Poisson_west{n,1} = Poisson_west;

% speed = speed(1:length(Burst_Results.spiketrain_burst_overall{n}));
% speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
% speed = gauss_smoothing(speed,10); % gaussian kernel smoothing
% % spiketrain_burst_overall = Burst_Results.spiketrain_burst_overall{n};
% % spiketrain_burst_overall = spiketrain_burst_overall';
% Burst_Results.speedScore_all(n,1) = speedScore;
% Burst_Results.speedScore_burst(n,1) =corr(speed,spiketrain_burst_overall,'rows','complete');

Burst_Results.mvl(n,1) = mvl;
Burst_Results.mvl_burst(n,1) = mvlburst;
Burst_Results.mv_arg(n,1) = mv_arg;
Burst_Results.mv_arg_burst(n,1) = mv_argburst;
Burst_Results.pref_angle_vonmises(n,1) = pref_dir_circ;
Burst_Results.pref_angle_burst_vonmises(n,1) = pref_dir_burst_circ;
Burst_Results.preferred_angle(n,1) = pref_angle;
Burst_Results.preferred_angle_burst(n,1) = pref_angleburst;
Burst_Results.HD_curve{n,1} = hdAxis;
Burst_Results.HD_curve{n,2} = hd_fr;
Burst_Results.HD_curve_Burst{n,1} = hdAxisburst;
Burst_Results.HD_curve_Burst{n,2} = hd_frburst;

Burst_Results.speedScore_overall(n,1) = speedScore_overall;
Burst_Results.speedScore_burst_overall(n,1) = speedScore_burst_overall;
Burst_Results.speedScore_north(n,1) = speedScore_north;
Burst_Results.speedScore_burst_north(n,1) = speedScore_burst_north;
Burst_Results.speedScore_south(n,1) = speedScore_south;
Burst_Results.speedScore_burst_south(n,1) = speedScore_burst_south;
Burst_Results.speedScore_east(n,1) = speedScore_east;
Burst_Results.speedScore_burst_east(n,1) = speedScore_burst_east;
Burst_Results.speedScore_west(n,1) = speedScore_west;
Burst_Results.speedScore_burst_west(n,1) = speedScore_burst_west;

Burst_Results.slope_overall(n,1) = slope_overall;
Burst_Results.slope_burst_overall(n,1) = slope_burst_overall;
Burst_Results.slope_north(n,1) = slope_north;
Burst_Results.slope_burst_north(n,1) = slope_burst_north;
Burst_Results.slope_south(n,1) = slope_south;
Burst_Results.slope_burst_south(n,1) = slope_burst_south;
Burst_Results.slope_east(n,1) = slope_east;
Burst_Results.slope_burst_east(n,1) = slope_burst_east;
Burst_Results.slope_west(n,1) = slope_west;
Burst_Results.slope_burst_west(n,1) = slope_burst_west;

Burst_Results.Spikefile_sq{n,1} = spikefile_sq; % Squish
Burst_Results.Percent_Burst_Overall_sq(n,1) = percent_overall_burst_sq;
Burst_Results.num_overall_burst_spikes_sq(n,1) = num_overall_burst_spikes_sq;
Burst_Results.spiketrain_burst_overall_sq{n,1} = spiketrain_burst_overall_sq;
Burst_Results.isi_overall_sq{n,1} = isi_overall_sq;
Burst_Results.FF_overall_sq(n,1) = FF_overall_sq;
Burst_Results.gamma_ci_overall_sq{n,1} = gamma_ci_overall_sq;
Burst_Results.Poisson_overall_sq{n,1} = Poisson_overall_sq;


Burst_Results.Percent_Burst_north_sq(n,1) = percent_north_burst_sq;
Burst_Results.num_north_burst_spikes_sq(n,1) = num_north_burst_spikes_sq;
Burst_Results.spiketrain_burst_north_sq{n,1} = spiketrain_burst_north_sq;
Burst_Results.isi_north_sq{n,1} = isi_north_sq;
Burst_Results.FF_north_sq(n,1) = FF_north_sq;
Burst_Results.gamma_ci_north_sq{n,1} = gamma_ci_north_sq;
Burst_Results.Poisson_north_sq{n,1} = Poisson_north_sq;

Burst_Results.Percent_Burst_south_sq(n,1) = percent_south_burst_sq;
Burst_Results.num_south_burst_spikes_sq(n,1) = num_south_burst_spikes_sq;
Burst_Results.spiketrain_burst_south_sq{n,1} = spiketrain_burst_south_sq;
Burst_Results.isi_south_sq{n,1} = isi_south_sq;
Burst_Results.FF_south_sq(n,1) = FF_south_sq;
Burst_Results.gamma_ci_south_sq{n,1} = gamma_ci_south_sq;
Burst_Results.Poisson_south_sq{n,1} = Poisson_south_sq;

Burst_Results.Percent_Burst_east_sq(n,1) = percent_east_burst_sq;
Burst_Results.spiketrain_burst_east_sq{n,1} = spiketrain_burst_east_sq;
Burst_Results.num_east_burst_spikes_sq(n,1) = num_east_burst_spikes_sq;
Burst_Results.isi_east_sq{n,1} = isi_east_sq;
Burst_Results.FF_east_sq(n,1) = FF_east_sq;
Burst_Results.gamma_ci_east_sq{n,1} = gamma_ci_east_sq;
Burst_Results.Poisson_east_sq{n,1} = Poisson_east_sq;

Burst_Results.Percent_Burst_west_sq(n,1) = percent_west_burst_sq;
Burst_Results.num_west_burst_spikes_sq(n,1) = num_west_burst_spikes_sq;
Burst_Results.spiketrain_burst_west_sq{n,1} = spiketrain_burst_west_sq;
Burst_Results.isi_west_sq{n,1} = isi_west_sq;
Burst_Results.FF_west_sq(n,1) = FF_west_sq;
Burst_Results.gamma_ci_west_sq{n,1} = gamma_ci_west_sq;
Burst_Results.Poisson_west_sq{n,1} = Poisson_west_sq;
% 
% speed_sq = speed(1:length(Burst_Results.spiketrain_burst_overall{n}));
% speed(isnan(speed)) = interp1(find(~isnan(speed)), speed(~isnan(speed)), find(isnan(speed)), 'pchip'); % interpolate NaNs
% speed = gauss_smoothing(speed,10); % gaussian kernel smoothing
% spiketrain_burst_overall = Burst_Results.spiketrain_burst_overall{n};
% % spiketrain_burst_overall = spiketrain_burst_overall';
% Burst_Results.speedScore_all_sq(n,1) = speedScore_sq;
% Burst_Results.speedScore_burst_sq(n,1) =corr(speed,spiketrain_burst_overall,'rows','complete');

Burst_Results.mvl_sq(n,1) = mvl_sq;
Burst_Results.mvl_burst_sq(n,1) = mvlburst_sq;
Burst_Results.mv_arg_sq(n,1) = mv_arg_sq;
Burst_Results.mv_arg_burst_sq(n,1) = mv_argburst_sq;
Burst_Results.preferred_angle_sq(n,1) = pref_angle_sq;
Burst_Results.preferred_angle_burst_sq(n,1) = pref_angleburst_sq;
Burst_Results.preferred_angle_vonmises_sq(n,1) = pref_dir_circ_sq;
Burst_Results.preferred_angle_burst_vonmises_sq(n,1) = pref_dir_burst_circ_sq;
Burst_Results.HD_curve_sq{n,1} = hdAxis_sq;
Burst_Results.HD_curve_sq{n,2} = hd_fr_sq;
Burst_Results.HD_curve_Burst_sq{n,1} = hdAxisburst_sq;
Burst_Results.HD_curve_Burst_sq{n,2} = hd_frburst_sq;

Burst_Results.speedScore_overall_sq(n,1) = speedScore_overall_sq;
Burst_Results.speedScore_burst_overall_sq(n,1) = speedScore_burst_overall_sq;
Burst_Results.speedScore_north_sq(n,1) = speedScore_north_sq;
Burst_Results.speedScore_burst_north_sq(n,1) = speedScore_burst_north_sq;
Burst_Results.speedScore_south_sq(n,1) = speedScore_south_sq;
Burst_Results.speedScore_burst_south_sq(n,1) = speedScore_burst_south_sq;
Burst_Results.speedScore_east_sq(n,1) = speedScore_east_sq;
Burst_Results.speedScore_burst_east_sq(n,1) = speedScore_burst_east_sq;
Burst_Results.speedScore_west_sq(n,1) = speedScore_west_sq;
Burst_Results.speedScore_burst_west_sq(n,1) = speedScore_burst_west_sq;

Burst_Results.slope_overall_sq(n,1) = slope_overall_sq;
Burst_Results.slope_burst_overall_sq(n,1) = slope_burst_overall_sq;
Burst_Results.slope_north_sq(n,1) = slope_north_sq;
Burst_Results.slope_burst_north_sq(n,1) = slope_burst_north_sq;
Burst_Results.slope_south_sq(n,1) = slope_south_sq;
Burst_Results.slope_burst_south_sq(n,1) = slope_burst_south_sq;
Burst_Results.slope_east_sq(n,1) = slope_east_sq;
Burst_Results.slope_burst_east_sq(n,1) = slope_burst_east_sq;
Burst_Results.slope_west_sq(n,1) = slope_west_sq;
Burst_Results.slope_burst_west_sq(n,1) = slope_burst_west_sq;

Burst_Results.old_score(n,1) = oldscore;
Burst_Results.old_score_sq(n,1) = oldscore_sq;
Burst_Results.old_scorey(n,1) = oldscorey;
Burst_Results.old_scorey_sq(n,1) = oldscorey_sq;
Burst_Results.old_scorex(n,1) = oldscorex;
Burst_Results.old_scorex_sq(n,1) = oldscorex_sq;
Burst_Results.old_slope(n,1) = oldslope;
Burst_Results.old_slope_sq(n,1) = oldslope_sq;
Burst_Results.old_slopey(n,1) = oldslopey;
Burst_Results.old_slopey_sq(n,1) = oldslopey_sq;
Burst_Results.old_slopex(n,1) = oldslopex;
Burst_Results.old_slopex_sq(n,1) = oldslopex_sq;

Burst_Results.old_score_burst(n,1) = oldburstscore;
Burst_Results.old_score_burst_sq(n,1) = oldburstscore_sq;
Burst_Results.old_scorey_burst(n,1) = oldburstscorey;
Burst_Results.old_scorey_burst_sq(n,1) = oldburstscorey_sq;
Burst_Results.old_scorex_burst(n,1) = oldburstscorex;
Burst_Results.old_scorex_burst_sq(n,1) = oldburstscorex_sq;
Burst_Results.old_slope_burst(n,1) = oldslopeburst;
Burst_Results.old_slope_burst_sq(n,1) = oldslopeburst_sq;
Burst_Results.old_slopey_burst(n,1) = oldslopebursty;
Burst_Results.old_slopey_burst_sq(n,1) = oldslopebursty_sq;
Burst_Results.old_slopex_burst(n,1) = oldslopeburstx;
Burst_Results.old_slopex_burst_sq(n,1) = oldslopeburstx_sq;
end

keyboard
