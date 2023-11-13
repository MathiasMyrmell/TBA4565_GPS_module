from receiver import Receiver
from prettytable import PrettyTable
from satellite import Satellite
import math
import numpy as np

def load_data(path):
    satellites = []
    try:
        file = open(path, "r")
        lines = file.readlines()
        for line in lines:
            rawData = line.split(" ")
            name = rawData[0]
            XA1 = rawData[1]
            YA1 = rawData[2]
            ZA1 = rawData[3]
            CPOA1 = rawData[4]
            XB1 = rawData[5]
            YB1 = rawData[6]
            ZB1 = rawData[7]
            CPOB1 = rawData[8]
            XA2 = rawData[9]
            YA2 = rawData[10]
            ZA2 = rawData[11]
            CPOA2 = rawData[12]
            XB2 = rawData[13]
            YB2 = rawData[14]
            ZB2 = rawData[15]
            CPOB2 = rawData[16]
            satellites.append(Satellite(name, XA1, YA1, ZA1, CPOA1, XB1, YB1, ZB1, CPOB1, XA2, YA2, ZA2, CPOA2, XB2, YB2, ZB2, CPOB2))
        file.close()
    except:
        print("file or game not found")
    return satellites

## Print help
def print_ellipsoidal_and_cartesian_coordinates(recA, recB):
    table = PrettyTable()
    table.title = "Ellipsoidal coordinates | Cartesian coordinates         "
    table.field_names = ["Receiver", "phi","lambda","height", "X", "Y", "Z"]
    table.add_row(["A", math.degrees(recA.phi), math.degrees(recA.lambda_), recA.h, recA.approx_X, recA.approx_Y, recA.approx_Z])
    table.add_row(["B", round(math.degrees(recB.phi),2), math.degrees(recB.lambda_), recB.h, recB.approx_X, recB.approx_Y, recB.approx_Z])
    table.align = "l"
    print(table)

def print_variance_covariance_matrix(VCV_matrix):
    table = PrettyTable()
    table.title = "Variance Covariance Matrix"
    for r in VCV_matrix:
        table.add_row(r)
    print(table)

def print_L_matrix(L_matrix):
    table = PrettyTable()
    table.title = "L matrix"
    table.field_names = ["L"]
    for r in L_matrix:
        table.add_row(r)
    table.align = "l"
    print(table)

def print_A_matrix(A_matrix):
    table = PrettyTable()
    table.title = "A matrix"
    table.field_names = ["AX", "AY", "AZ", "l1", "l2", "l3", "l4"]
    for r in A_matrix:
        table.add_row(r)
    table.align = "l"
    print(table)

def print_X_matrix(X_matrix):
    table = PrettyTable()
    table.field_names = ["X matrix"]
    for r in X_matrix:
        table.add_row(r)
    table.align = "l"
    print(table)

def print_variance_covariance_matrix(VCV_matrix):
    table = PrettyTable()
    table.title = "Variance Covariance Matrix"
    for r in VCV_matrix:
        table.add_row(r)
    print(table)




###Task 2
def calculate_estimated_receiver_position(unknown_receiver, known_receiver, satellites, t1, t2):
    print("-Calculating estimated receiver position-")
    unknown_receiver.estimated_X = unknown_receiver.approx_X
    unknown_receiver.estimated_Y = unknown_receiver.approx_Y
    unknown_receiver.estimated_Z = unknown_receiver.approx_Z

    delta_xr = 1
    delta_yr = 1
    delta_zr = 1
    error = 8*10**-9
    i=0
    while abs(delta_xr) >error or abs(delta_yr) >error or abs(delta_zr) >error:
        x_matrix = X_matrix(satellites, [known_receiver, unknown_receiver], t1, t2)
        # print_X_matrix(x_matrix)
        delta_xr = x_matrix[0][0]
        delta_yr = x_matrix[1][0]
        delta_zr = x_matrix[2][0]
        unknown_receiver.estimated_X += delta_xr
        unknown_receiver.estimated_Y += delta_yr
        unknown_receiver.estimated_Z += delta_zr
        i+=1

    print("Iterations:", i) 

    positions = PrettyTable()
    positions.field_names = (["","X", "Y", "Z"])
    positions.add_row(["Estimated (New)",unknown_receiver.estimated_X, unknown_receiver.estimated_Y, unknown_receiver.estimated_Z])
    positions.add_row(["Approximated (Old)",unknown_receiver.approx_X, unknown_receiver.approx_Y, unknown_receiver.approx_Z])
    positions.add_row(["Difference", unknown_receiver.estimated_X-unknown_receiver.approx_X, unknown_receiver.estimated_Y-unknown_receiver.approx_Y, unknown_receiver.estimated_Z-unknown_receiver.approx_Z])
    positions.align = "l"
    positions.title = "Receiver position (m)"
    print(positions)

#X matrix
def X_matrix(satellites, receivers, t1, t2):
    L = L_matrix(satellites, receivers, t1, t2)
    # print_L_matrix(L)
    A = A_matrix(satellites, receivers[1], t1, t2)
    # print_A_matrix(A)
    X = np.linalg.inv(A.T.dot(A)).dot(A.T).dot(L)
    # print_X_matrix(X)
    return X

#L matrix
def L_matrix(satellites, receivers, t1, t2):
    SI = satellites[0]
    times = [t1, t2]
    delta_L = [[], []]
    for i in range(0, len(times)):
        t = times[i]
        for j in range(1, len(satellites)):
            S = satellites[j]
            L = calc_L([SI, S], receivers, t)
            delta_L[i].append(L)
    
    L = np.array([
        [delta_L[0][0]],
        [delta_L[0][1]],
        [delta_L[0][2]],
        [delta_L[0][3]],
        [delta_L[1][0]],
        [delta_L[1][1]],
        [delta_L[1][2]],
        [delta_L[1][3]]
    ])
    return L

def calc_L(satellites, receivers, t):
    RA = receivers[0]
    RB = receivers[1]
    SI = satellites[0]
    SJ = satellites[1]

    phi = _calc_phi(satellites, receivers, t)
    p_JB0 = _calc_p0(SJ, RB, t)
    p_IB0 = _calc_p0(SI, RB, t)
    p_JA0 =_calc_p0(SJ, RA, t)
    p_IA0 = _calc_p0(SI, RA, t)

    L = phi - p_JB0 + p_IB0 + p_JA0 - p_IA0
    return L

def _calc_phi(satellites, receivers, t):
    RA = receivers[0]
    RB = receivers[1]
    SI = satellites[0]
    SJ = satellites[1]
    phi_AI = SI.get_carrier_phase_observation(RA.name, t)
    phi_AJ = SJ.get_carrier_phase_observation(RA.name, t)
    phi_BI = SI.get_carrier_phase_observation(RB.name, t)
    phi_BJ = SJ.get_carrier_phase_observation(RB.name, t)

    phi = phi_BJ-phi_BI-phi_AJ+phi_AI
    return phi

def _calc_p0(satellite, receiver, t):
    #X^j = saetllite
    #X_i0 = receiver
    satC = satellite.get_coordinates(receiver.name, t)
    recC = receiver.get_estimated_coordinates_cartesian()
    p0 = math.sqrt(((satC[0]-recC[0])**2) + ((satC[1]-recC[1])**2) + ((satC[2]-recC[2])**2))
    return p0

#A matrix
def A_matrix(satellites, receiver, t1, t2):
    l = math.radians(receiver.lambda_)
    print(receiver.lambda_)
    times = [t1, t2]
    #Calculate AX
    SI = satellites[0]
    AX = [[], []]
    for i in range(0,len(times)):
        t = times[i]
        for j in range(1, len(satellites)):
            S = satellites[j]
            A = calc_AX([SI, S], receiver, t)
            AX[i].append(A)
    
    #Calculate AY
    AY = [[], []]
    for i in range(0,len(times)):
        t = times[i]
        for j in range(1, len(satellites)):
            S = satellites[j]
            A = calc_AY([SI, S], receiver, t)
            AY[i].append(A)

    #Calculate AZ
    AZ = [[], []]
    for i in range(0,len(times)):
        t = times[i]
        for j in range(1, len(satellites)):
            S = satellites[j]
            A = calc_AZ([SI, S], receiver, t)
            AZ[i].append(A)

    #A matrix
    A = np.array([
        [AX[0][0], AY[0][0], AZ[0][0], l, 0, 0, 0],
        [AX[0][1], AY[0][1], AZ[0][1], 0, l, 0, 0],
        [AX[0][2], AY[0][2], AZ[0][2], 0, 0, l, 0],
        [AX[0][3], AY[0][3], AZ[0][3], 0, 0, 0, l],
        [AX[1][0], AY[1][0], AZ[1][0], l, 0, 0, 0],
        [AX[1][1], AY[1][1], AZ[1][1], 0, l, 0, 0],
        [AX[1][2], AY[1][2], AZ[1][2], 0, 0, l, 0],
        [AX[1][3], AY[1][3], AZ[1][3], 0, 0, 0, l]
    ])

    return A

def calc_AX(satellites, receiver, time):
    SI = satellites[0]
    SJ = satellites[1]
    R = receiver

    XJ = SJ.get_coordinates(R.name, time)[0]
    PJR = _calc_p0(SJ, R, time)

    XI = SI.get_coordinates(R.name, time)[0]
    PIR = _calc_p0(SI, R, time)

    RX = R.get_estimated_coordinates_cartesian()[0]

    AX = -(XJ-RX)/PJR + (XI-RX)/PIR
    return AX

def calc_AY(satellites, receiver, time):
    SI = satellites[0]
    SJ = satellites[1]
    R = receiver

    YJ = SJ.get_coordinates(R.name, time)[1]
    PJR = _calc_p0(SJ, R, time)

    YI = SI.get_coordinates(R.name, time)[1]
    PIR = _calc_p0(SI, R, time)

    RY = R.get_estimated_coordinates_cartesian()[1]
    AY = -(YJ-RY)/PJR + (YI-RY)/PIR
    return AY

def calc_AZ(satellites, receiver, time):
    SI = satellites[0]
    SJ = satellites[1]
    R = receiver
    
    ZJ = SJ.get_coordinates(R.name, time)[2]
    PJR = _calc_p0(SJ, R, time)

    ZI = SI.get_coordinates(R.name, time)[2]
    PIR = _calc_p0(SI, R, time)

    RZ = R.get_estimated_coordinates_cartesian()[2]

    AZ = -(ZJ-RZ)/PJR + (ZI-RZ)/PIR
    return AZ

#VarianceCovarianceMatrix
def calculate_variance_covariance_matrix(satellites, receiver, t1, t2):
    A = A_matrix(satellites, receiver, t1, t2)
    Q = np.linalg.inv(A.T.dot(A))
    return Q

def calculate_PDOP(VCVM):
    qxx = VCVM[0][0]
    qyy = VCVM[1][1]
    qzz = VCVM[2][2]

    PDOP = math.sqrt(qxx+qyy+qzz)
    return PDOP
###Task 3


###Task 4
def calculate_geodetic_coordinates(receiver):
    height, phi = _calculate_geodetic_height(receiver)
    _lambda = math.atan(receiver.estimated_Y/receiver.estimated_X)

    geodetic = PrettyTable()
    geodetic.field_names = (["phi","lambda", "h"])
    geodetic.add_row([math.degrees(phi), math.degrees(_lambda), height])
    geodetic.align = "l"
    geodetic.title = "Estimated receiver position, geodetic coordinates"
    print(geodetic)

def _calculate_geodetic_height(receiver):
    X = receiver.estimated_X
    Y = receiver.estimated_Y
    Z = receiver.estimated_Z

    #Step 1: Calc p
    p = math.sqrt(X**2+Y**2)

    #Step 2: Calc phi0
    phi0 = math.atan((Z/p)*(1-((receiver.b**2)/(receiver.a**2))))

    i = 0
    while True:
        #Step 3: Calc N0
        N0 = receiver.a**2/(math.sqrt(receiver.a**2*(math.cos(phi0)**2)+receiver.b**2*(math.sin(phi0)**2)))

        #Step 4: Calc h
        h = (p/math.cos(phi0))-N0

        #Step 5: Improved phi0
        e = (receiver.a**2-receiver.b**2)/receiver.a**2
        improved_phi0 = math.atan((Z/p)*(1-e*(N0/(N0+h)))**-1)
        

        #Step 6: Compare approximate and improved phi
        diff = abs(phi0-improved_phi0)
        phi0 = improved_phi0
        
        i+=1
        if(diff==0):
            break
    

    print("Iterations:", i)

    return h, phi0


if __name__ == "__main__":
    path = "project2/data.txt"
    satellites = load_data(path)

    #ReceiverA (Known station)
    phiA = -32.003884648
    lambdaA = 115.894802001
    heightA = 23.983
    known = True
    RA = Receiver("A",phiA, lambdaA, heightA, known)

    #ReceiverB (Unknown station)
    phiB = -31.9
    lambdaB = 115.75
    heightB = 50
    known = False
    RB = Receiver("B",phiB, lambdaB, heightB, known)

    ###Task 1
    #Transform receiver coordinates from ellipsoidal to cartesian coordinates
    print("------Task 1------")
    print("-Transiforming receiver coordinates from ellipsoidal to cartesian coordinates-")
    RA.calculate_cartesian_coordinates()
    RB.calculate_cartesian_coordinates()
    print_ellipsoidal_and_cartesian_coordinates(RA, RB)
    #Correct

    ###Task 2
    #Estimate position of receiver B
    print("------Task 2------")
    print("-Estimate position of receiver B-")
    t1 = 172800
    t2 = 175020
    calculate_estimated_receiver_position(RB, RA, satellites, t1, t2) #Unknown receiver, known receiver, satellites, time 1, time 2

    print("-b ) Estimate the variance-covariance matrix-")
    Q = calculate_variance_covariance_matrix(satellites, RB, t1, t2)
    print_variance_covariance_matrix(Q)

    print("- Calculating PDOP")
    PDOP = calculate_PDOP(Q)
    print("PDOP:", PDOP)

    # ###Task 3
    # print("------Task 3------")
    # print("-a: Do nothing and use the real ambiguity in the next step-")


    # print("-b: Fix the real ambiguity to nearest integer value-")


    # print("-c: Use the standard deviation of ambiguities computed from LS in step 2 and fix the ambiguities considering the standard deviations.")


    ###Task 4
    print("------Task 4------")
    calculate_geodetic_coordinates(RB) #Receiver. Only needs X, Y, Z

