clear; close all; clc;

%% read satellite ephemeris data from rinex data
ephem_t = readtable("C:\Users\urrut\OneDrive\Documents\GitHub\INS\DDKU\dat\mlfp_2024-01-06T22_00_00.000000000\orbit_param.csv");
ephem_t.time = AtlasLib.UTC2gps(ephem_t.time)';

%% Construct Satellite Positions and Velocities
[t, ia] = unique(ephem_t.time);
numBatch = lenth(ia);
for ii = 1:numBatch - 1
    
end
