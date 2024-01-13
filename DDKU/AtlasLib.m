%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atlas Library: Functions that deal with the space segment of navigation.
%                This includes creating satellite positions, atmospheric 
%                corrections, etc. 
%
% Dependencies: None
%
% References:
% 1. Calculation of Satellite Position from Ephemeris Data. Applies GPS for 
%    Engineers and Project Managers
%
% Author: Kevin Urrutia 
%
% Revision History:
% v 1.0 Pre-Release Jan 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef AtlasLib
    properties
        %constants
        GM = 3.986005e14; % [m^3/s^2] WSG-84 value for the product of graviational Constant and the Mass of the Earth
        omega_e = 7.292115e-5 %[rad/s] WSG-84 value of Earth's rotation rate
        half_week = 302400; %[s] number of seconds that represent half of GPS week
    end

    methods (Static)
        function [r, r_dot] = satellite_cart_coord(obj, ephem)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Inputs: 
            % ephem - satellite ephemeris
            %
            % Outputs:
            % r - ECEF Cartesian Coordinates of satellite positions
            % r_dot - ECEF Cartesian Coordinates of satellite velocities
            %
            % Author:
            % Kevin Urrutia 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            SVClockBias = ephem.SVClockBias;
            SVClockDrift = ephem.SVClockDrift;
            SVClockDriftRate = ephem.SVClockDriftRate;
            IODE = ephem.IODE;
            Crs = ephem.Crs;
            DeltaN = ephem.DeltaN;
            M0 = ephem.M0;
            Cuc = ephem.Cuc;
            Eccentricity = ephem.Eccentricity;
            Cus = ephem.Cus;
            sqrtA = ephem.sqrtA;
            Toe = ephem.Toe; 
            Cic = ephem.Cic;
            Omega0 = ephem.Omega0;
            Io = ephem.Io;
            Crc = ephem.Crc;
            omega = ephem.omega;
            OmegaDot = ephem.OmegaDot;
            IDOT = ephem.IDOT;
            GPSWeek = ephem.GPSWeek;
            TGD = ephem.TGD;
            TransTime = ephem.TransTime;

            r = zeros(3,1);
            r_dot = zeros(3,1);

            delta_t = obj.weekroll(TransTime - Toe); 
            delta_t_sv = SVClockBias + SVClockDrift*delta_t + SVClockDriftRate*delta_t^2;
            t_sv = TransTime - delta_t_sv; %true satellite time

            A = sqrtA^2;
            t_k = obj.weekroll(t_sv - Toe);
            n0 = sqrt(obj.GM / A^3); %Mean Motion
            n = n0 + DeltaN; % Corrected Mean Motion
            M_k = M0 + n*t_k; %Mean Anomoly

            E = zeros(1, 10); %kepler's equation of eccentricity to be solved through iteration
            E(1) = M_k;
            for i = 2:length(E)
                E(i) = M_k + Eccentricity*sin(E(i-1));
            end
            E_k = E(end);

            nu_k = 2*atan2(sqrt(1 + Eccentricity)*(sin(E_k)/2), sqrt(1+Eccentricity)*(cos(E_k)/2)); %obtain true anomoly
            Phi_k = nu_k + omega;  %Argument of Latitude

            delta_u_k = Cuc*cos(2*Phi_k) + Cus*sin(2*Phi_k); %Argument of latitude correction
            delta_r_k = Crc*cos(2*Phi_k) + Crs*sin(2*Phi_k); %Radius Correction
            delta_i_k = Cic*cost(2*Phi_k) + Cis*sin(2*Phi_k); %Inclinnation Correction
            u_k = Phi_k + delta_u_k; %Corrected Argument of Latitude
            r_k = A*(1 - Eccentricity*cos(E_k)) + delta_r_k; %Corrected Radius
            i_k = Io + delta_i_k + IDOT*t_k; %corrected incclination
            X_k_prime = r_k*cos(u_k); %Position in orbital plane
            Y_k_prime = r_k*sin(u_k); %Position in orbital Plane
            Omega_k = Omega0 + (OmegaDot - obj.omega_e)*t_k + obj.omega_e*t_k; %corrected longitude of ascending node

            r(1) = X_k_prime*cos(Omega_k) - Y_k_prime*sin(Omega_k)*cos(i_k); %x ECEF
            r(2) = X_k_prime*sin(Omega_k) + Y_k_prime*cos(Omega_k)*cos(i_k); %y ECEF
            r(3) = Y_k_prime*sin(i_k); %z ECEF

            

        end

        function t = weekroll(obj, time)
            t = time;
            if time > obj.half_week
                t = time - 2*obj.half_week;
            elseif time < -obj.half_week
                t = time + 2*obj.half_week;
            end

        end
    end
end