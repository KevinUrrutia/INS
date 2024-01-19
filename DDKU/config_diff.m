%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Config Diff: Configurations for double differences 
% 
% Inputs: 
%        data_num - selector for the data set to pull from
%
% Output: 
%        p - parameters for the double differencing algorithm
%
% Dependencies: AtlasLib
%
% Author: Kevin Urrutia 
%
% Revision History:
% v 1.0 Pre-Release Jan 17, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = config_diff(data_num)
    parent_dir = '/Data/';
    addpath(strcat('..', parent_dir));
   

    switch data_num
        case 0
            p.f_base_o = strcat(parent_dir, 'AustinUTX/Base/SEPT1290.19O'); %base observation file
            p.f_base_p = strcat(parent_dir, 'AustinUTX/Base/SEPT1290.19P'); %base mixed nav file
            p.base_LLA = [-97.736925, 30.290841, 191.724448]; %[deg, deg, m];
            p.base_ECEF = AtlasLib.LLA2ECEF(AtlasLib, p.base_LLA);
            
            p.f_rover_o = strcat(parent_dir, 'AustinUTX/Rover/SEPT1291.19O'); %rover observation file
            p.f_rover_p = strcat(parent_dir, 'AustinUTX/Rover/SEPT1290.19P'); %rover mixed nav file
            p.f_rover_gt = strcat(parent_dir, 'AustinUTX/Rover/ground_truth.log'); %rover ground truth
            p.f_rover_imu = strcat(parent_dir, 'AustinUTX/Rover/imu.log'); %rover imu file
            p.imu_type = 'ATLANS-C';

            p.f_type = 'rinex_v3';
           
        otherwise
            error("ERROR: Select a valid data path!");
    end
end