import math
import numpy as np
from satellite import Satellite
from receiver import Receiver
from prettytable import PrettyTable

###Task 1###
## Compute the coordinates of the 7 GPS satellites.

## Reads data from file on path, and returns a list of satellites
def read_data_from_file(path):
    satellites = []
    try:
        file = open(path, "r")
        lines = file.readlines()
        for line in lines:
            rawData = line.split(" ")
            name = rawData[0]
            year = rawData[1]
            month = rawData[2]
            day = rawData[3]
            hour = rawData[4]
            minute = rawData[5]
            second = rawData[6]
            SVClockBias = rawData[7]
            SVClockDrift = rawData[8]
            SVClockDriftRate = rawData[9]
            IODE = rawData[10]
            Crs = rawData[11]
            DeltaN = rawData[12]
            M0 = rawData[13]
            Cuc = rawData[14]
            e = rawData[15]
            Cus = rawData[16]
            sqrtA = rawData[17]
            Toe = rawData[18]
            Cic = rawData[19]
            Omega0 = rawData[20]
            Cis = rawData[21]
            i0 = rawData[22]
            Crc = rawData[23]
            omega = rawData[24]
            OmegaDot = rawData[25]
            IDOT = rawData[26]
            CodesL2 = rawData[27]
            GPSWeek = rawData[28]
            L2Pflag = rawData[29]
            SVacc = rawData[30]
            SVhealth = rawData[31]
            TGD = rawData[32]
            IODC = rawData[33]
            Ttom = rawData[34]
            P = rawData[35]
            DT = rawData[36]
            DION = rawData[37]
            DTROP = rawData[38]
            satellites.append(Satellite(name, year, month, day, hour, minute, second, SVClockBias, SVClockDrift, SVClockDriftRate, IODE, Crs, DeltaN, M0, Cuc, e, Cus, sqrtA, Toe, Cic, Omega0, Cis, i0, Crc, omega, OmegaDot, IDOT, CodesL2, GPSWeek, L2Pflag, SVacc, SVhealth, TGD, IODC, Ttom, P, DT, DION, DTROP))
        file.close()
    except:
        print("file or game not found")
    return satellites

def create_table_satellites(coordinates, title):
    t = PrettyTable()
    t.field_names = (["Sat Name", "X", "Y", "Z"])
    t.align = "l"
    t.title = title
    for key in coordinates:
        sat_name = key
        x = coordinates[key][0]
        y = coordinates[key][1]
        z = coordinates[key][2]
        t.add_row([sat_name, x, y, z])

    # for sat in satellites:
    #     name = sat.name
    #     x = sat.x
    #     y = sat.y
    #     z = sat.z
    #     t.add_row([name, x, y, z])
    #     # print(sat.name)
    #     # print(calculated_coordinates[0])
    return t

def calculate_difference_between_coordinates(coordinates, coordinates0, title):
    t = PrettyTable()
    t.field_names = (["Sat Name", "diff X", "diff Y", "diff Z"])
    t.align = "l"
    t.title = title
    for key in coordinates:
        sat_name = key
        x = round(abs(coordinates[key][0] - coordinates0[key][0]),4)
        y = round(abs(coordinates[key][1] - coordinates0[key][1]),4)
        z = round(abs(coordinates[key][2] - coordinates0[key][2]),4)
        t.add_row([sat_name, x, y, z])
    return t

if "__main__" == __name__:
    path = "project1/data.txt"
    ## Read data from file
    print("Reading data from file...")
    satellites = read_data_from_file(path)
    print("Done reading data from file.\n")


    ###Task 1
    print("-------Task 1---------")
    print("Computing coordinates for each satellite...")
    coordinates = {}
    for sat in satellites:
        coordinates[sat.name] = sat.calculate_coordinates()
    tb = create_table_satellites(coordinates, "Coordinates for every satellite")
    print(tb)

    ##Fasit for SV03
    X = 23098418.6200952
    Y = -12669381.2111253
    Z = 2686141.36110425
    hr= 176.356 
    for sat in satellites:
        if sat.name == "G03":
            x = sat.x
            y = sat.y
            z = sat.z
            print("Diff X-value: ", abs(X-x))
            print("Diff Y-value: ", abs(Y-y))
            print("Diff Z-value: ", abs(Z-z))


    # ###Task 2
    # print("\n-------Task 2---------")
    # print("Setting the 9 correction terms to 0...")
    # coordinates0 = {}
    # for sat in satellites:
    #     sat.set_coorection_terms_to_zero()
    #     coordinates0[sat.name] = sat.calculate_coordinates()
    #     sat.reset_correction_terms()

    # tb0 = create_table_satellites(coordinates0, "Coordinates for every satellite with correction terms=0")
    # print(tb0)
    # print("\n")

    # print("Difference between the coordinates (abs value):")
    # tbDiff = calculate_difference_between_coordinates(coordinates, coordinates0, "Difference between coordinates (abs value)")
    # print(tbDiff)


    # ###Task 3
    # print("\n-------Task 3---------")
    # print("Transform the approximate receiver coordinates to approximate Cartesian coordinates")
    # # Approximate geodetic receiver coordinates
    # phi_r0 = 63.2
    # lambda_r0 = 10.2
    # h_r0 = 100
    # r = Receiver(phi_r0, lambda_r0, h_r0)
    # print("Calculate cartesian coordinates...")
    # r.transform_coordinates_from_ellipsodial_to_cartesion_coordinates()
    # print("Done calculating cartesian coordinates.")
    # print("Approximate receiver coordinates:")
    # tbReceiver = r.create_position_table()
    # print(tbReceiver)


    # ###Task 4
    # print("\n-------Task 4---------")
    # print("Estimated receiver position:")
    # table = r.calculate_estimated_receiver_position(satellites)
    # print(table)

    # ###Task 5
    # print("-------Task 5---------")
    # print("Calculate variance-covariance matrix...")
    # q_values, table = r.calculate_variance_covariance(satellites)
    # print("Done calculating variance-covariance matrix.")
    # print("Variance-covariance matrix:")
    # print(table)

    # PDOP = r.calculate_PDOP(q_values)
    # print("PDOP: ", PDOP)

