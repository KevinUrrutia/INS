import georinex as gr
import numpy as np
import sys, getopt
import os

def read_orbit(orbit_file, meas_dir):
    file_name = "orbit.csv"
    orbit = gr.load(orbit_file)
    os.chdir(meas_dir)

    with open(file_name, 'w') as f:
        f.write("time, sat, X, Y, Z, Vx, Vy, Vz, clock, dclock\n")

        for i in range(len(orbit.time)):
            for j in range(len(orbit.sv)):
                sat = orbit.sv[j]
                sat_str = np.array2string(sat)

                time_str = np.datetime_as_string(orbit.time[i])
                X_str = str(orbit.sel(sv=sat, ECEF="x")['position'].values[i])
                Y_str = str(orbit.sel(sv=sat, ECEF="y")['position'].values[i])
                Z_str = str(orbit.sel(sv=sat, ECEF="z")['position'].values[i])
                Vx_str = str(orbit.sel(sv=sat, ECEF="x")['velocity'].values[i])
                Vy_str = str(orbit.sel(sv=sat, ECEF="y")['velocity'].values[i])
                Vz_str = str(orbit.sel(sv=sat, ECEF="z")['velocity'].values[i])
                clock = str(orbit.sel(sv=sat)['clock'].values[i])
                dclock = str(orbit.sel(sv=sat)['dclock'].values[i])

                f.write(time_str + ',' + sat_str + ',' + X_str + ',' + Y_str + ',' + Z_str + ",")
                f.write(Vx_str + "," + Vy_str + "," + Vz_str + "," + clock + "," + dclock + "\n")

        os.chdir("../")




def main(argv):

    #read file names from command line arguments
   observation_file = ""
   navigation_file = ""
   orbit_file = ""
   opts, args = getopt.getopt(argv,"ho:n:r:",["obsfile=", "navfile=", "orbitfile="])

   for opt, arg in opts:
       if opt == '-h':
           print('rinex_converter.py -o <observation file> -n <navigation file> -r <orbit file>')
           sys.exit()
       elif opt in ("-o", "--obsfile"):
           observation_file = arg
       elif opt in ("-n", "--navfile"):
           navigation_file = arg
       elif opt in ('-r', "--orbitfile"):
           orbit_file = arg

   print("Observation file: ", observation_file)
   print("Navigation file: ", navigation_file)
   print("Orbit file:", orbit_file)

   #create a directory to store the converted rinex data into
   indx = observation_file.index('.')
   file_dir = observation_file[0:indx - 1] + '/'
   base_loc = observation_file[0:4]

   observation_file = file_dir + observation_file
   navigation_file = file_dir + navigation_file
   orbit_file = file_dir + orbit_file


   nav = gr.load(navigation_file)
   meas_dir = base_loc + "_" +  np.datetime_as_string(nav.time[0])
   if not os.path.exists(meas_dir):
    os.mkdir(meas_dir)

   read_orbit(orbit_file, meas_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
