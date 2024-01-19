%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_dd: Tests the Double Difference Algorithm
%
% Dependencies: AtlasLib, config_diff
%
% Author: Kevin Urrutia 
%
% Revision History:
% v 1.0 Pre-Release Jan 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% load configurations
%Data set numbering: Austin_UTX: 0, 
data_set = 0;
p = config_diff(data_set); 

%% load data into super struct
S_s = loadData(p);