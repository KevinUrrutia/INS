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
        deg2rad = pi / 180; 
        rad2deg = 180 / pi;
        GM = 3.986005e14; % [m^3/s^2] WSG-84 value for the product of graviational Constant and the Mass of the Earth
        omega_e = 7.292115e-5 %[rad/s] WSG-84 value of Earth's rotation rate
        half_week = 302400; %[s] number of seconds that represent half of GPS week
        a = 6378137.0; %[m] WSG-84 Equitorial Radius/ semi-major axis
        f = 0.00335281; % WSG-84 flatness
        b = 6378137.0*(1-0.00335281); %[m] a*(1-f) WSG-84 Semi-minor axis
        e = sqrt(0.00335281*(2-0.00335281)); %sqrt(f*(2-f)) WSG-84 Eccentricity

    end

    methods (Static)
        function [r, r_dot] = satellite_ephem(obj, ephem)
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
            Cis = ephem.Cis;
            Io = ephem.Io;
            Crc = ephem.Crc;
            omega = ephem.omega;
            OmegaDot = ephem.OmegaDot;
            IDOT = ephem.IDOT;
            TransTime = ephem.TransTime;

            r = zeros(3,1);
            r_dot = zeros(3,1);

            %satellite positions in ECEF
            delta_t = obj.weekroll(obj, TransTime - Toe); 
            delta_t_sv = SVClockBias + SVClockDrift*delta_t + SVClockDriftRate*delta_t^2;
            t_sv = TransTime - delta_t_sv; %true satellite time

            A = sqrtA^2;
            t_k = obj.weekroll(obj, t_sv - Toe);
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
            delta_i_k = Cic*cos(2*Phi_k) + Cis*sin(2*Phi_k); %Inclinnation Correction
            u_k = Phi_k + delta_u_k; %Corrected Argument of Latitude
            r_k = A*(1 - Eccentricity*cos(E_k)) + delta_r_k; %Corrected Radius
            i_k = Io + delta_i_k + IDOT*t_k; %corrected incclination
            X_k_prime = r_k*cos(u_k); %Position in orbital plane
            Y_k_prime = r_k*sin(u_k); %Position in orbital Plane
            Omega_k = Omega0 + (OmegaDot - obj.omega_e)*t_k + obj.omega_e*t_k; %corrected longitude of ascending node

            r(1) = X_k_prime*cos(Omega_k) - Y_k_prime*sin(Omega_k)*cos(i_k); %x ECEF
            r(2) = X_k_prime*sin(Omega_k) + Y_k_prime*cos(Omega_k)*cos(i_k); %y ECEF
            r(3) = Y_k_prime*sin(i_k); %z ECEF

            %satellite velocities in ECEF
            E_k_dot = n / (1 - Eccentricity*cos(E_k)); %Eccentricity anomoly rate
            nu_k_dot = E_k_dot*sqrt(1 - Eccentricity^2) / (1 - Eccentricity*cos(E_k)); %True Anomoly Rate
            d_ik_dt = IDOT + 2*nu_k*(Cis*cos(2*Phi_k) - Cic*sin(2*Phi_k)); %Corrected Inclination Angle rate
            u_k_dot = nu_k_dot + 2*nu_k_dot*(Cus*cos(2*Phi_k) - Cuc*sin(2*Phi_k)); %corrected Argument of Latitude Rate
            r_k_dot = Eccentricity*A*E_k_dot*sin(E_k) + 2*nu_k*(Crs*cos(2*Phi_k) - Crc*sin(2*Phi_k)); %Corrected Radius rate
            Omega_k_dot = OmegaDot - obj.omega_e;
            X_k_prime_dot = r_k_dot*cos(u_k) - r_k*u_k_dot*sin(u_k);
            Y_k_prime_dot = r_k_dot*sin(u_k) + r_k*u_k_dot*cos(u_k);
            r_dot(1) = -X_k_prime*Omega_k_dot*sin(Omega_k) + X_k_prime_dot*cos(Omega_k) - Y_k_prime_dot*sin(Omega_k)*cos(i_k) - Y_k_prime*(Omega_k_dot*cos(Omega_k)*cos(i_k) - d_ik_dt*sin(Omega_k)*sin(i_k));
            r_dot(2) = X_k_prime*Omega_k_dot*cos(Omega_k) + X_k_prime_dot*sin(Omega_k) + Y_k_prime_dot*cos(Omega_k)*cos(i_k) - Y_k_prime*(Omega_k_dot*sin(Omega_k)*cos(i_k) + d_ik_dt*cos(Omega_k)*sin(i_k));
            r_dot(3) = Y_k_prime_dot*sin(i_k) + Y_k_prime*d_ik_dt*cos(i_k);

        end

        function t = weekroll(obj, time)
            t = time;
            if time > obj.half_week
                t = time - 2*obj.half_week;
            elseif time < -obj.half_week
                t = time + 2*obj.half_week;
            end

        end

        function leaps = getleaps()
            leaps = [46828800, 78364801, 109900802, 173059203, 252028804, 315187205, 346723206, 393984007, 425520008, 457056009, 504489610, 551750411, 599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017];
        end

        function isleap = checkLeapYear(gpsTime)
            leaps = AtlasLib.getleaps();
            isleap = ismember(gpsTime, leaps);
        end

        function nleaps = countLeaps(gpsTime, dirFlag)
            leaps = AtlasLib.getleaps();
            lenLeaps = length(leaps);
            nleaps = 0;
            for  i = 1:lenLeaps
                if strcmp(dirFlag, 'unix2gps')
                    if(gpsTime >= leaps(i) - i)
                        nleaps = nleaps + 1;
                    end
                elseif strcmp(dirFlag, 'gps2unix')
                    if(gpsTime >= leaps(i))
                        nleaps = nleaps + 1;
                    end
                else
                    error('Invalid Time Conversion')
                end
            end
        end

        function gpsTime = unix2gps(unixTime)
            if(mod(unixTime, 1) ~= 0)
                unixTime = unixTime - 0.5;
                isLeap = 1;
            else 
                isLeap = 0;
            end

            gpsTime = unixTime - 315964800;
            nleaps = AtlasLib.countLeaps(gpsTime, 'unix2gps');
            gpsTime = gpsTime + nleaps + isLeap;
        end

        function unixTime = gps2unix(gpsTime)
            unixTime = gpsTime + 315964800;
            nleaps = AtlasLib.countLeaps(gpsTime, 'gps2unix');
            unixTime = unixTime - nleaps;
            if(AtlasLib.isleap(gpsTime))
                unixTime = unixTime + 0.5;
            end
        end

        function unixTime = UTC2unix(UTC)
            unixTime = convertTo(UTC, 'epochtime', 'Epoch', '1970-01-01')';
        end

        function gpsTime = UTC2gps(UTC)
            unixTime = AtlasLib.UTC2unix(UTC);
            gpsTime = AtlasLib.unix2gps(unixTime);
        end

        function r_ECEF = LLA2ECEF(obj, LLA)
            r_ECEF = zeros(3, 1);
            lat = LLA(1) * obj.deg2rad;
            lon = LLA(2) * obj.deg2rad;
            alt = LLA(3);

            R_N = obj.a / (1 - obj.e^2*sin(lat)^2)^(1/2);

            r_ECEF(1) = (R_N + alt)*cos(lat)*cos(lon);
            r_ECEF(2) = (R_N + alt)*cos(lat)*sin(lon);
            r_ECEF(3) = (R_N*(1 - obj.e^2) + alt) * sin(lat);
        end

        function LLA = ECEF2LLA(obj, r_ECEF)
            x = r_ECEF(1);
            y = r_ECEF(2);
            z = r_ECEF(3);

            lon = atan2(y, x);
            
            %initialization
            h = 0;
            d_h = inf;
            R_N = obj.a;
            p = sqrt(x^2 + y^2);

            while(abs(d_h) > 1e-9) %wait for convergence
                denom = (1 - obj.e^2)*R_N + h;
                sin_phi = z / denom;
                lat = atan((z + obj.e^2 * R_N*sin_phi)/p);
                R_N = obj.a / sqrt(1 - obj.e^2*sin_phi^2);
                h_k = (p / cos(lat)) - R_N;
                d_h = h_k - h;
                h =h_k;
            end

            LLA(1) = lat * obj.rad2deg;
            LLA(2) = lon * obj.rad2deg;
            LLA(3) = h;
        end

    end
end