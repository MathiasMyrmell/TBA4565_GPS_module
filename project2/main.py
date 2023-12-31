from receiver import Receiver
from prettytable import PrettyTable
from satellite import Satellite
import math
import numpy as np
import copy
#Global variables
L1 = 1575.42*10**6
c = 299792458
_lambda = c/L1

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
    phase_ambiguity =[_lambda]*(len(satellites)-1)
    
    unknown_receiver.estimated_X = unknown_receiver.approx_X
    unknown_receiver.estimated_Y = unknown_receiver.approx_Y
    unknown_receiver.estimated_Z = unknown_receiver.approx_Z

    delta_xr = 1
    delta_yr = 1
    delta_zr = 1

    error = 8*10**-9
    i=0
    while abs(delta_xr) >error or abs(delta_yr) >error or abs(delta_zr) >error:
        L = L_matrix(satellites, [known_receiver, unknown_receiver], t1, t2)
        A = A_matrix(satellites, unknown_receiver, t1, t2, phase_ambiguity)
        P = P_matrix()
        X = np.linalg.inv(A.T.dot(P).dot(A)).dot(A.T).dot(P).dot(L)
        delta_xr = X[0][0]
        delta_yr = X[1][0]
        delta_zr = X[2][0]
        unknown_receiver.estimated_X += delta_xr
        unknown_receiver.estimated_Y += delta_yr
        unknown_receiver.estimated_Z += delta_zr
        unknown_receiver.matrices = [X, L, A, P]
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

    Q = calculate_variance_covariance_matrix(satellites, unknown_receiver, t1, t2, phase_ambiguity)
    print_variance_covariance_matrix(Q)
    PDOP = calculate_PDOP(Q)
    unknown_receiver.PDOP = PDOP
    print("PDOP:", PDOP)
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

    phi = (_lambda)*(phi_BJ-phi_BI-phi_AJ+phi_AI)
    return phi

def _calc_p0(satellite, receiver, t):
    satC = satellite.get_coordinates(receiver.name, t)
    recC = receiver.get_estimated_coordinates_cartesian()
    p0 = math.sqrt(((satC[0]-recC[0])**2) + ((satC[1]-recC[1])**2) + ((satC[2]-recC[2])**2))
    return p0

#A matrix
def A_matrix(satellites, receiver, t1, t2, phase_ambiguity):
    times = [t1, t2]
    SI = satellites[0]
    #Calculate AX
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

    ## A matrix
    A = np.array([
        [AX[0][0], AY[0][0], AZ[0][0], phase_ambiguity[0],  0,                  0,                  0                   ],
        [AX[0][1], AY[0][1], AZ[0][1], 0,                   phase_ambiguity[1], 0,                  0                   ],
        [AX[0][2], AY[0][2], AZ[0][2], 0,                   0,                  phase_ambiguity[2], 0                   ],
        [AX[0][3], AY[0][3], AZ[0][3], 0,                   0,                  0,                  phase_ambiguity[3]  ],
        [AX[1][0], AY[1][0], AZ[1][0], phase_ambiguity[0],  0,                  0,                  0                   ],
        [AX[1][1], AY[1][1], AZ[1][1], 0,                   phase_ambiguity[1], 0,                  0                   ],
        [AX[1][2], AY[1][2], AZ[1][2], 0,                   0,                  phase_ambiguity[2], 0                   ],
        [AX[1][3], AY[1][3], AZ[1][3], 0,                   0,                  0,                  phase_ambiguity[3]  ]
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

#P matrix
def P_matrix():
    sigma = 0.002
    factor = (1/(2*(sigma**2)))*(1/5)
    matrix = np.array([
        [4, -1, -1, -1, 0, 0, 0, 0],
        [-1, 4, -1, -1, 0, 0, 0, 0],
        [-1, -1, 4, -1, 0, 0, 0, 0],
        [-1, -1, -1, 4, 0, 0, 0, 0],
        [0, 0, 0, 0, 4, -1, -1, -1],
        [0, 0, 0, 0, -1, 4, -1, -1],
        [0, 0, 0, 0, -1, -1, 4, -1],
        [0, 0, 0, 0, -1, -1, -1, 4]
    ])
    return factor*matrix

#VarianceCovarianceMatrix
def calculate_variance_covariance_matrix(satellites, receiver, t1, t2, phase_ambiguity):
    A = A_matrix(satellites, receiver, t1, t2, phase_ambiguity)
    P = P_matrix()
    Q = np.linalg.inv(A.T.dot(P).dot(A))
    return Q

def calculate_PDOP(VCVM):
    qxx = VCVM[0][0]
    qyy = VCVM[1][1]
    qzz = VCVM[2][2]
    PDOP = math.sqrt(qxx+qyy+qzz)
    return PDOP

###Task 3
def estimate_receiver_coordinates_real_ambiguity(receiver):
    X = receiver.matrices[0]
    L = receiver.matrices[1]
    A = receiver.matrices[2]
    P = receiver.matrices[3]

    n_values = X[3:7]
    n_values = n_values*_lambda
    for i in range(0, len(n_values)):
        L[i][0]-=(n_values[i])
        L[i+4][0]-=(n_values[i])
    A = np.delete(A, [3,4,5,6], 1)
    X = np.linalg.inv(A.T.dot(P).dot(A)).dot(A.T).dot(P).dot(L)
    receiver.estimated_X+=X[0][0]
    receiver.estimated_Y+=X[1][0]
    receiver.estimated_Z+=X[2][0]
    receiver.matrices = [X, L, A, P]

    table = PrettyTable()
    table.title = "Estimated Coordinates, Real ambiguity"
    table.field_names = ["X", "Y", "Z"]
    table.add_row([receiver.estimated_X, receiver.estimated_Y, receiver.estimated_Z])
    table.align = "l"
    print(table)
def estimate_receiver_coordinates_integer_ambiguity(receiver):
    X = receiver.matrices[0]
    L = receiver.matrices[1]
    A = receiver.matrices[2]
    P = receiver.matrices[3]

    n_values = X[3:7]
    for i in range(0, len(n_values)):
        n_values[i][0] = round(n_values[i][0])
    n_values = n_values*_lambda
    for i in range(0, len(n_values)):
        L[i][0]-=(n_values[i])
        L[i+4][0]-=(n_values[i])
 
    A = np.delete(A, [3,4,5,6], 1)
    X = np.linalg.inv(A.T.dot(P).dot(A)).dot(A.T).dot(P).dot(L)
    receiver.estimated_X+=X[0][0]
    receiver.estimated_Y+=X[1][0]
    receiver.estimated_Z+=X[2][0]
    receiver.matrices = [X, L, A, P]

    table = PrettyTable()
    table.title = "Estimated Coordinates, Integer ambiguity"
    table.field_names = ["X", "Y", "Z"]
    table.add_row([receiver.estimated_X, receiver.estimated_Y, receiver.estimated_Z])
    table.align = "l"
    print(table)
 
def estimate_receiver_coordinates_std_ambiguity(receiver):
    X = receiver.matrices[0]
    A = receiver.matrices[2]
    P = receiver.matrices[3]
    CX = np.linalg.inv(A.T.dot(P).dot(A))
    n_values1 = X[3:7]
    n_values2 = X[3:7]
    for i in range(0, len(n_values1)):
        n_values1[i][0] = round(n_values1[i][0]-CX[i][i])
    for i in range(0, len(n_values2)):
        n_values2[i][0] = round(n_values2[i][0]+CX[i][i])
    
    if n_values1[0] == n_values2[0] and n_values1[1] == n_values2[1] and n_values1[2] == n_values2[2] and n_values1[3] == n_values2[3]:
        print("Results will be the same as for integer ambiguity")
    else:
        print("Results will be different than for integer ambiguity")

def calculate_variance_covariance(receiver):
    A = receiver.matrices[2]
    P = receiver.matrices[3]
    CX = np.linalg.inv(A.T.dot(P).dot(A))
    print_variance_covariance_matrix(CX)
    PDOP = calculate_PDOP(CX)
    print("PDOP:", PDOP)

if __name__ == "__main__":
    path        = "project2/data.txt"
    satellites  = load_data(path)

    #ReceiverA (Known station)
    phiA    = -32.003884648
    lambdaA = 115.894802001
    heightA = 23.983
    known   = True
    RA      = Receiver("A", phiA, lambdaA, heightA, known)

    #ReceiverB (Unknown station)
    phiB    = -31.9
    lambdaB = 115.75
    heightB = 50
    known   = False
    RB      = Receiver("B", phiB, lambdaB, heightB, known)

    ###Task 1
    print("------Task 1------")
    print("-Transiforming receiver coordinates from ellipsoidal to cartesian coordinates-")
    RA.calculate_cartesian_coordinates()
    RB.calculate_cartesian_coordinates()
    print_ellipsoidal_and_cartesian_coordinates(RA, RB)

    ###Task 2
    print("\n------Task 2------")
    print("-Estimate position of receiver B-")
    t1 = 172800
    t2 = 175020
    calculate_estimated_receiver_position(RB, RA, satellites, t1, t2) 

    ###Task 3
    print("\n------Task 3------")
    print("---------a--------")
    print("-Do nothing and use the real ambiguity-")
    receiver_float = copy.deepcopy(RB)
    estimate_receiver_coordinates_real_ambiguity(receiver_float)
    calculate_variance_covariance(receiver_float)

    print("\n--------b--------")
    print("-Fix the real ambiguity to nearest integer value-")
    receiver_fix = copy.deepcopy(RB)
    estimate_receiver_coordinates_integer_ambiguity(receiver_fix)
    calculate_variance_covariance(receiver_fix)

    print("\n------c------")
    print("-Use the standard deviation of ambiguities computed from LS in step 2 and fix the ambiguities considering the standard deviations.")
    receiver_std = copy.deepcopy(RB)
    estimate_receiver_coordinates_std_ambiguity(receiver_std)

    # ##Task 4
    print("\n------Task 4------")
    chosen_receiver = receiver_fix
    chosen_receiver.calculate_estimated_height()