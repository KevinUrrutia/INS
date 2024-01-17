clear; close all; clc;

%% read satellite ephemeris data from rinex data
ephem_t = readtable("C:\Users\urrut\OneDrive\Documents\GitHub\INS\DDKU\dat\mlfp_2024-01-06T22_00_00.000000000\orbit_param.csv");
ephem_t.time = AtlasLib.UTC2gps(ephem_t.time)';
ephem_t = rmmissing(ephem_t, 1);

%% read satellite positions for validation
nav_t = readtable("C:\Users\urrut\OneDrive\Documents\GitHub\INS\DDKU\dat\mlfp_2024-01-06T22_00_00.000000000\orbit.csv");
nav_t.time = AtlasLib.UTC2gps(nav_t.time)';
nav_t = rmmissing(nav_t, 1);

%% read satellite observation data 
obs_t = readtable("C:\Users\urrut\OneDrive\Documents\GitHub\INS\DDKU\dat\mlfp_2024-01-06T22_00_00.000000000\obs.csv");
obs_t.time = AtlasLib.UTC2gps(obs_t.time)';

%% Construct Satellite Positions and Velocities
Atlas = AtlasLib;
[t, ia] = unique(ephem_t.time);
numBatch = length(ia);
r_S = struct([]);
for ii = 1:numBatch - 1
    prns = ephem_t.PRN(ia(ii) : ia(ii + 1) - 1);
    for jj = 0:length(prns) - 1
        if(~isfield(r_S, genvarname(prns{jj + 1})))
            r_S(1).(genvarname(prns{jj + 1})) = [];
            r_S(1).(genvarname(prns{jj + 1})).GPStime = [];
            r_S(1).(genvarname(prns{jj + 1})).r_sv = [];
            r_S(1).(genvarname(prns{jj + 1})).v_sv = [];
            r_S(1).(genvarname(prns{jj + 1})).SNR = [];

        end
        r_S(1).(genvarname(prns{jj + 1})).GPStime = [r_S(1).(genvarname(prns{jj + 1})).GPStime, t(ii)];
        ephem = struct('SVClockBias', ephem_t.SVClockBias(ii + jj), ...
                        'SVClockDrift', ephem_t.SVClockDrift(ii + jj), ...
                        'SVClockDriftRate', ephem_t.SVClockDriftRate(ii + jj), ...
                        'Crs', ephem_t.Crs(ii + jj), ...
                        'DeltaN', ephem_t.DeltaN(ii + jj), ...
                        'M0', ephem_t.M0(ii + jj), ...
                        'Cuc', ephem_t.Cuc(ii + jj), ...
                        'Eccentricity', ephem_t.Eccentricity(ii +  jj), ...
                        'Cus', ephem_t.Cus(ii + jj), ...
                        'sqrtA', ephem_t.sqrtA(ii + jj), ...
                        'Toe', ephem_t.Toe(ii + jj), ...
                        'Cic', ephem_t.Cic(ii + jj), ...
                        'Omega0', ephem_t.Omega0(ii + jj), ...
                        'Cis', ephem_t.Cis(ii + jj), ...
                        'Io', ephem_t.Io(ii + jj), ...
                        'Crc', ephem_t.Crc(ii + jj), ...
                        'omega', ephem_t.omega(ii + jj), ...
                        'OmegaDot', ephem_t.OmegaDot(ii + jj), ...
                        'IDOT', ephem_t.IDOT(ii + jj), ...
                        'TransTime', ephem_t.TransTime(ii + jj));
        [r_sv, v_sv] = AtlasLib.satellite_ephem(Atlas, ephem);
        r_S(1).(genvarname(prns{jj + 1})).r_sv = [r_S(1).(genvarname(prns{jj + 1})).r_sv, r_sv];
        r_S(1).(genvarname(prns{jj + 1})).v_sv = [r_S(1).(genvarname(prns{jj + 1})).v_sv, v_sv];
    end
end

%% Interpolate Satellite Positions to 10 Hz
fn = fieldnames(r_S);
for ii = 1:size(fn, 1)
    t_new = r_S(1).(fn{ii}).GPStime(1):1:r_S(1).(fn{ii}).GPStime(end);
    r_S(1).(fn{ii}).r_sv = interp1(double(r_S(1).(fn{ii}).GPStime), r_S(1).(fn{ii}).r_sv', double(t_new), 'spline');
    r_S(1).(fn{ii}).v_sv = interp1(double(r_S(1).(fn{ii}).GPStime), r_S(1).(fn{ii}).v_sv', double(t_new), 'spline')';
    r_S(1).(fn{ii}).GPStime = t_new;
end

%% Grab appropriate SNR and place them at the correct time instant, interpolate if needed
% for ii = 1:size(fn, 1)
%     r_S(1).(fn{ii}).SNR = [];
% 
%     %find the indices that correspond to current PRN
%     idx = any(strcmp(obs_t.PRN, fn{ii}));
% 
%     for jj = 1:length(r_S(1).(fn{ii}).GPStime)
%        time_idx = find(cell2mat(obs_t.time(idx)))
%     end
% end
