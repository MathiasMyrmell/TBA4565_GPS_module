# Import External libraries
from prettytable import PrettyTable
# Import classes
from satellite import Satellite
from receiver import Receiver

## Reads data from file in path, and returns a list of satellite objects
def read_data_from_file(path):
    satellites = []
    try:
        file = open(path, "r")
        lines = file.readlines()
        for line in lines:
            data = line.split(" ")
            satellites.append(Satellite(data))
        file.close()
    except:
        print("file or game not found")
    return satellites

## Creates table to visualize the coordinates for each satellite
# For task 1 and 2
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
    return t

## Creates table to visualize the difference between the coordinates for each satellite
## with or without correction terms. Task 2
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
    ## Path to data file
    path = "project1/data.txt"

    ## Read data from file
    print("-Reading data from file...")
    satellites = read_data_from_file(path)
    print("-Done reading data from file.\n")


    ###Task 1
    print("-------Task 1---------")
    print("-Computing coordinates for each satellite...")
    coordinates = {}
    for sat in satellites:
        coordinates[sat.name] = sat.calculate_coordinates()
    tb = create_table_satellites(coordinates, "Coordinates for every satellite")
    print(tb)

  
    ###Task 2
    print("\n-------Task 2---------")
    print("-Setting the 9 correction terms to 0...")
    coordinates0 = {}
    for sat in satellites:
        sat.set_coorection_terms_to_zero()
        coordinates0[sat.name] = sat.calculate_coordinates()
        sat.reset_correction_terms()

    tb0 = create_table_satellites(coordinates0, "Coordinates for every satellite with correction terms=0")
    print(tb0)
    print("\n")

    print("-Difference between the coordinates (abs value):")
    tbDiff = calculate_difference_between_coordinates(coordinates, coordinates0, "Difference between coordinates (abs value)")
    print(tbDiff)


    ###Task 3
    print("\n-------Task 3---------")
    print("-Transform the approximate receiver coordinates to approximate Cartesian coordinates")
    # Approximate geodetic receiver coordinates
    phi_r0 = 63.2
    lambda_r0 = 10.2
    h_r0 = 100
    r = Receiver(phi_r0, lambda_r0, h_r0)
    print("-Calculate cartesian coordinates...")
    r.transform_coordinates_from_ellipsodial_to_cartesion_coordinates()
    print("-Done calculating cartesian coordinates.")
    print("-Approximate receiver coordinates:")
    tbReceiver = r.create_position_table()
    print(tbReceiver)


    ###Task 4
    print("\n-------Task 4---------")
    print("-Estimated receiver position:")
    r.calculate_estimated_receiver_position(satellites)


    ###Task 5
    print("-------Task 5---------")
    print("-Calculate variance-covariance matrix...")
    q_values, table = r.calculate_variance_covariance(satellites)
    print("-Done calculating variance-covariance matrix.")
    print("-Variance-covariance matrix:")
    print(table)
    print("\n")

    print("-Calculate PDOP...")
    PDOP = r.calculate_PDOP(q_values)
    print("-PDOP: ", PDOP)


    ###Task 6
    print("Calculating height from estimated coordinates...")
    r.calculate_estimated_height()
    print("Done calculating height from estimated coordinates.")

    ###Task 7
# 63.41547533840656 | 10.405195308850008 | 119.88603231683373 | -1.5860074209543168e-08
# 63.415463330158545 | 10.405973021493388 | 176.37842045351863 | -8.306456802888356e-08
