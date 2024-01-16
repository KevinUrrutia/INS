import georinex as gr
import numpy as np
import sys, getopt
import warnings
import os

def read_orbit(orbit_file, meas_dir):
    file_name = "orbit.csv"
    orbit = gr.load(orbit_file)
    os.chdir(meas_dir)

    with open(file_name, 'w') as f:
        f.write("time, PRN, X, Y, Z, Vx, Vy, Vz, clock, dclock\n")

        for i in range(len(orbit.time)):
            for j in range(len(orbit.sv)):
                sat = orbit.sv[j]
                sat_str = np.array2string(sat)
                sat_str = sat_str[1:len(sat_str) -1] #remove quotes

                time_str = np.datetime_as_string(orbit.time[i])
                X_str = str(orbit.sel(sv=sat, ECEF="x")['position'].values[i])
                Y_str = str(orbit.sel(sv=sat, ECEF="y")['position'].values[i])
                Z_str = str(orbit.sel(sv=sat, ECEF="z")['position'].values[i])
                Vx_str = str(orbit.sel(sv=sat, ECEF="x")['velocity'].values[i])
                Vy_str = str(orbit.sel(sv=sat, ECEF="y")['velocity'].values[i])
                Vz_str = str(orbit.sel(sv=sat, ECEF="z")['velocity'].values[i])
                clock = str(orbit.sel(sv=sat)['clock'].values[i])
                dclock = str(orbit.sel(sv=sat)['dclock'].values[i])


                if 'G' in sat_str:
                    idx = sat_str.index('G');
                    sat_str = "GPS_L1x_" + sat_str[idx + 1:idx + 3];
                elif 'E' in sat_str:
                    idx = sat_str.index('E')
                    sat_str = "GAL_E1a_" + sat_str[idx + 1:idx + 3];
                elif 'R' in sat_str:
                    idx = sat_str.index('R')
                    sat_str = "GLO_R1x_" + sat_str[idx + 1:idx + 3];


                f.write(time_str + ',' + sat_str + ',' + X_str + ',' + Y_str + ',' + Z_str + ",")
                f.write(Vx_str + "," + Vy_str + "," + Vz_str + "," + clock + "," + dclock + "\n")

        os.chdir("../")

def read_nav(navigation_file, meas_dir):
    file_name = "orbit_param.csv"
    nav = gr.load(navigation_file)

    os.chdir(meas_dir)

    with open(file_name, 'w') as f:
        f.write("time,PRN,SVClockBias,SVClockDrift,SVClockDriftRate,IODE,Crs,DeltaN,M0,Cuc,Eccentricity,Cus,sqrtA,")
        f.write("Toe,Cic,Omega0,Cis,Io,Crc,omega,OmegaDot,IDOT,CodesL2,GPSWeek,L2Pflag,")
        f.write("SVacc,health,TGD,IODC,TransTime,FitIntvl\n")

        for i in range(len(nav.time)):
            for j in range(len(nav.sv)):
                sat = nav.sv[j]
                sat_str = np.array2string(sat)
                sat_str = sat_str[1:len(sat_str) -1] #remove quotes

                time_str = np.datetime_as_string(nav.time[i])
                SVClockBias = str(nav.sel(sv=sat)['SVclockBias'].values[i])
                SVClockDrift = str(nav.sel(sv=sat)['SVclockDrift'].values[i])
                SVClockDriftRate = str(nav.sel(sv=sat)['SVclockDriftRate'].values[i])
                IODE = str(nav.sel(sv=sat)['IODE'].values[i])
                Crs = str(nav.sel(sv=sat)['Crs'].values[i])
                DeltaN = str(nav.sel(sv=sat)['DeltaN'].values[i])
                M0 = str(nav.sel(sv=sat)['M0'].values[i])
                Cuc = str(nav.sel(sv=sat)['Cuc'].values[i])
                Eccentricity = str(nav.sel(sv=sat)['Eccentricity'].values[i])
                Cus = str(nav.sel(sv=sat)['Cus'].values[i])
                sqrtA = str(nav.sel(sv=sat)['sqrtA'].values[i])
                Toe = str(nav.sel(sv=sat)['Toe'].values[i])
                Cic = str(nav.sel(sv=sat)['Cic'].values[i])
                Omega0 = str(nav.sel(sv=sat)['Omega0'].values[i])
                Cis = str(nav.sel(sv=sat)['Cis'].values[i])
                Io = str(nav.sel(sv=sat)['Io'].values[i])
                Crc = str(nav.sel(sv=sat)['Crc'].values[i])
                omega = str(nav.sel(sv=sat)['omega'].values[i])
                OmegaDot = str(nav.sel(sv=sat)['OmegaDot'].values[i])
                IDOT = str(nav.sel(sv=sat)['IDOT'].values[i])
                CodesL2 = str(nav.sel(sv=sat)['CodesL2'].values[i])
                GPSWeek = str(nav.sel(sv=sat)['GPSWeek'].values[i])
                L2Pflag = str(nav.sel(sv=sat)['L2Pflag'].values[i])
                SVacc = str(nav.sel(sv=sat)['SVacc'].values[i])
                health = str(nav.sel(sv=sat)['health'].values[i])
                TGD = str(nav.sel(sv=sat)['TGD'].values[i])
                IODC = str(nav.sel(sv=sat)['IODC'].values[i])
                TransTime = str(nav.sel(sv=sat)['TransTime'].values[i])
                FitIntvl = str(nav.sel(sv=sat)['FitIntvl'].values[i])

                if 'G' in sat_str:
                    idx = sat_str.index('G');
                    sat_str = "GPS_L1x_" + sat_str[idx + 1:idx + 3];
                elif 'E' in sat_str:
                    idx = sat_str.index('E')
                    sat_str = "GAL_E1a_" + sat_str[idx + 1:idx + 3];
                elif 'R' in sat_str:
                    idx = sat_str.index('R')
                    sat_str = "GLO_R1x_" + sat_str[idx + 1:idx + 3];

                f.write(time_str + "," + sat_str + "," + SVClockBias + "," + SVClockDrift + "," + SVClockDriftRate + ",")
                f.write(IODE + "," + Crs + "," + DeltaN + "," + M0 + "," + Cuc + "," + Eccentricity + "," + Cus + ",")
                f.write(sqrtA + "," + Toe + "," + Cic + "," + Omega0 + "," + Cis + "," + Io + "," + Crc + "," + omega + "," + OmegaDot + ",")
                f.write(IDOT + "," + CodesL2 + "," + GPSWeek + "," + L2Pflag + "," + SVacc + "," + health + "," + TGD + ",")
                f.write(IODC + "," + TransTime + "," + FitIntvl + '\n')

    file_name = "atmosphere_correction.csv"
    with open(file_name, 'w') as f:
        f.write("a0,a1,a2,a3,b0,b1,b2,b3\n")
        for i in range(len(nav.ionospheric_corr_GPS)):
            if i < len(nav.ionospheric_corr_GPS) - 1:
                f.write(str(nav.ionospheric_corr_GPS[i]) + ',')
            else:
                f.write(str(nav.ionospheric_corr_GPS[i]) + '\n')

    os.chdir("../")

def read_obs(obs_file, meas_dir):
    file_name = "obs.csv"

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        obs = gr.load(obs_file)

    os.chdir(meas_dir)

    with open(file_name, 'w') as f:
        f.write("time,PRN,rho,phi,SNR\n")

        for i in range(len(obs.time)):
            for j in range(len(obs.sv)):
                sat = obs.sv[j]
                sat_str = np.array2string(sat)
                sat_str = sat_str[1:len(sat_str) -1] #remove quotes

                time_str = np.datetime_as_string(obs.time[i])
                rho = str(obs.sel(sv=sat)['C1'].values[i])
                phi = str(obs.sel(sv=sat)['L1'].values[i])
                SNR = str(obs.sel(sv=sat)['S1'].values[i])

                if 'G' in sat_str:
                    idx = sat_str.index('G');
                    sat_str = "GPS_L1x_" + sat_str[idx + 1:idx + 3];
                elif 'E' in sat_str:
                    idx = sat_str.index('E')
                    sat_str = "GAL_E1a_" + sat_str[idx + 1:idx + 3];
                elif 'R' in sat_str:
                    idx = sat_str.index('R')
                    sat_str = "GLO_R1x_" + sat_str[idx + 1:idx + 3];

                f.write(time_str + "," + sat_str +  "," + rho + "," + phi  + "," + SNR + "\n")

    os.chdir("../")

def read_coord(coord_file, meas_dir):
    file_name = "coord.csv"
    str = "GEOID HEIGHT"

    with open(coord_file, 'r') as coord:
        for l_no, line in enumerate(coord):
            if str in line:
                X_str = coord.readline()
                idx_start = X_str.index("X  - ")
                idx_end = X_str.index(" (m")
                X_str = X_str[idx_start + 5:idx_end]

                Y_str = coord.readline()
                idx_start = Y_str.index("Y  - ")
                idx_end = Y_str.index(" (m")
                Y_str = Y_str[idx_start + 5:idx_end]

                Z_str = coord.readline()
                idx_start = Z_str.index("Z  - ")
                idx_end = Z_str.index(" (m")
                Z_str = Z_str[idx_start + 5:idx_end]
                break


    os.chdir(meas_dir)
    with open(file_name, 'w') as f:
        f.write("X,Y,Z\n")
        f.write(X_str + "," + Y_str + "," + Z_str + "\n")

    os.chdir("../")

#####################################################################################################################################################
# Rinex Converter: Converts Rinex V 2.1 Files gathered from CORS servers to csv
# Inputs: Rinex V 2.1 Observation, Navigation, Orbit and Coordinate Files
#         in the form of SSSDDDH.YYo, SSSDDDH.YYn, .sp3, SSSS.ds respectively
#         where:
#         SSSS represents the site id
#         DDD represents the day of the year starting with 001 as Jan 1st
#         H represents the hour of the day
#         YY represents the 2-digit year ex. 2024 -> 24
#
# Output: atmosphere_correction.csv - contains GPS atmosphere correction params a0 - a1, b0 - b1
#         coord.csv - contains ECEF coordinates in meters for base station
#         obs.csv - contains psuedorange, carrier pahse and SNR per satellite at UTC time intervals
#         orbit_param.csv - contains satellite keplarian, clock, health parameters at UTC time intervals
#         orbit.csv - contains satellite positions and velocity in ECEF coordinates in meters and clock correcrtion paramters at UTC time intervals
#
# Dependencies: georinex
#
# Author: Kevin Urrutia
#
# Revision History:
# v 1.0 Pre-Release Jan 11, 2024
 #####################################################################################################################################################



def main(argv):

    #read file names from command line arguments
   observation_file = ""
   navigation_file = ""
   orbit_file = ""
   coord_file = ""
   opts, args = getopt.getopt(argv,"ho:n:r:c:",["obsfile=", "navfile=", "orbitfile=", "coordfile="])

   for opt, arg in opts:
       if opt == '-h':
           print('rinex_converter.py -o <observation file> -n <navigation file> -r <orbit file> -c <coordinate file>')
           sys.exit()
       elif opt in ("-o", "--obsfile"):
           observation_file = arg
       elif opt in ("-n", "--navfile"):
           navigation_file = arg
       elif opt in ('-r', "--orbitfile"):
           orbit_file = arg
       elif opt in('-c', '--coordfile'):
           coord_file = arg

   print("Observation file: ", observation_file)
   print("Navigation file: ", navigation_file)
   print("Orbit file: ", orbit_file)
   print("Coordinate file: ", coord_file)

   #create a directory to store the converted rinex data into
   indx = observation_file.index('.')
   file_dir = observation_file[0:indx - 1] + '/'
   base_loc = observation_file[0:4]

   observation_file = file_dir + observation_file
   navigation_file = file_dir + navigation_file
   orbit_file = file_dir + orbit_file
   coord_file = file_dir + coord_file


   nav = gr.load(navigation_file)
   meas_dir = base_loc + "_" +  np.datetime_as_string(nav.time[0])
   meas_dir = meas_dir.replace(":", "_")
   if not os.path.exists(meas_dir):
    os.mkdir(meas_dir)

   read_orbit(orbit_file, meas_dir)
   print("Finished Reading Orbit File")
   read_nav(navigation_file, meas_dir)
   print("Finished Reading Navigation File")
   read_obs(observation_file, meas_dir)
   print("Finished Reading Observation File")
   read_coord(coord_file, meas_dir)
   print("Finished Reading Coordinate File")


if __name__ == "__main__":
    main(sys.argv[1:])
