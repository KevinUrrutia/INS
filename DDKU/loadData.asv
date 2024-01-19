%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loadData: Loads in different types of navigation and observation files
%           and creates a uniform struct
% 
% Inputs: 
%        p - configuration struct that contains meta data like data
%        location
%
% Output: 
%        S_s - Super struct that maintains observations and satellite
%        positions
%
% Dependencies: AtlasLib
%
% Author: Kevin Urrutia 
%
% Revision History:
% v 1.0 Pre-Release Jan 17, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S_s = loadData(p)
    switch(p.f_type)
        case 'rinex_v3'
            obs_b = rinexread(p.f_base_o);
            nav_b = rinexread(p.f_base_p);

            obs_r = rinexread(p.f_rover_o);
            nav_r = rinexread(p.f_rover_p);
            imu_r = readtable(p.f_rover_imu, 'VariableNamingRule', 'preserve');
            gt_r = readtable(p.f_rover_gt, 'VariableNamingRule', 'preserve');

            %convert from UTC time to GPS time
            obs_b.GPS.time = AtlasLib.UTC2gps(obs_b.GPS.Time)';
            nav_b.GPS.time = AtlasLib.UTC2gps(nav_b.GPS.Time)';

            obs_r.GPS.time = AtlasLib.UTC2gps(obs_r.GPS.Time)';
            nav_r.GPS.time = AtlasLib.UTC2gps(nav_b.GPS.Time)';

            %construct satellite positions and velocities 
            
        otherwise
            error("ERROR! Invalid File type entered")
    end
end