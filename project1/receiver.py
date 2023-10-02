import math
import numpy as np
from prettytable import PrettyTable

class Receiver:

    def __init__(self, phi_r0, lambda_r0, h_r0):
        self.phi_r0 = self.degree_to_radians(phi_r0)
        self.lambda_r0 = self.degree_to_radians(lambda_r0)
        self.h_r0 = h_r0
        ##From internettet
        self.a = 6378137
        self.b = 6356752.314245
        self.c = 299792458
        ##Calculated
        self.approxX = None
        self.approxY = None
        self.approxZ = None
        self.apporxN = None
        self.p = None
        self.phi0 =None
        self.N0 = None
        self.h = None
        self.improved_phi0 = None
        #Estimated
        self.estimatedX = None
        self.estimatedY = None
        self.estimatedZ = None
        self.receiver_clock_correction = None
        self.estimatedHeight = None

    def degree_to_radians(self, degree):
        return math.radians(degree)
    
    def radians_to_degree(self, radians):
        return math.degrees(radians)
    
    ###Task 3
    ## Transform coordinates from ellipsodial to cartesion coordinates
    def transform_coordinates_from_ellipsodial_to_cartesion_coordinates(self):
        # Calculate N
        N = (self.a**2)/(math.sqrt((self.a**2)*math.cos(self.phi_r0)**2 + (self.b**2)*math.sin(self.phi_r0)**2))
        # # Transform x
        x = (N+self.h_r0)*math.cos(self.phi_r0)*math.cos(self.lambda_r0)
        # # Transform y
        y = (N+self.h_r0)*math.cos(self.phi_r0)*math.sin(self.lambda_r0)
        # # Transform z
        z = ((self.b**2/self.a**2)*N+self.h_r0)*math.sin(self.phi_r0)

        self.approxX = x
        self.approxY = y
        self.approxZ = z
        self.apporxN = N
        return N,x,y,z
    
    ## Display result for Task 3
    def create_position_table(self):
        t = PrettyTable()
        t.field_names = (["X", "Y", "Z"])
        t.add_row([self.approxX, self.approxY, self.approxZ])
        t.align = "l"
        t.title = "Approximate receiver coordinates"
        return t
    




    ###Task 4
    def calculate_estimated_receiver_position(self, satellites):
        x_matrix = self.X(satellites)
        delta_x = x_matrix[0][0]
        delta_y = x_matrix[1][0]
        delta_z = x_matrix[2][0]
        receiver_clock_correction = x_matrix[3][0]

        self.estimatedX = self.approxX + delta_x
        self.estimatedY = self.approxY + delta_y
        self.estimatedZ = self.approxZ + delta_z
        self.receiver_clock_correction = receiver_clock_correction

        self.calculate_estimated_height()
        self.estimatedHeight = self.h
        print("phi", self.radians_to_degree(self.phi0))
        print("lambda", self.radians_to_degree(math.atan(-self.estimatedY/self.estimatedX)))

        table = PrettyTable()
        table.field_names = (["Variable", "Value"])
        table.add_row(["X", self.estimatedX])
        table.add_row(["Y", self.estimatedY])
        table.add_row(["Z", self.estimatedZ])
        table.add_row(["T", self.receiver_clock_correction])
        table.add_row(["Height", self.estimatedHeight])

        return table

    ###Observation Equation in the Least Squares Method
    ## L matrix
    def L(self, satellites):
        rows = []
        for sat in satellites:
            rows.append([sat.P-self.p_i0(sat)-self.c*sat.DT-sat.DION - sat.DTROP])
        matrix = np.array(rows)
        return matrix

    ## A matrix
    def A(self, satellites):
        rows = []
        for sat in satellites:
            rows.append([self._ax(sat), self._ay(sat), self._az(sat),-self.c])
        matrix = np.array(rows)
        return matrix
    
    # Help functions for A matrix
    def _ax(self, satellite):
        return -(satellite.x-self.approxX)/self.p_i0(satellite)

    def _ay(self, satellite):
        return -(satellite.y-self.approxY)/self.p_i0(satellite)
    
    def _az(self, satellite):
        return -(satellite.z-self.approxZ)/self.p_i0(satellite)

    ## X matrix
    def X(self, satellites):
        L = self.L(satellites)
        A = self.A(satellites)
        firstPart = np.linalg.inv(np.matmul(A.T,A))
        secondPart = np.matmul(A.T,L)
        X = np.matmul(firstPart,secondPart)
        # X = np.matmul(np.linalg.inv(A),L)
        # X = np.linalg.inv(A.T@A)@A.T@L
        return X


    def p_i0(self, satellite):
        satX = satellite.x
        satY = satellite.y
        satZ = satellite.z

        recX = self.approxX
        recY = self.approxY
        recZ = self.approxZ
        x = (satX-recX)**2
        y = (satY-recY)**2
        z = (satZ-recZ)**2

        pi0 = math.sqrt(x+y+z)
        return pi0


    ## Transformation from cartesion to ellipsodial
    def calculate_estimated_height(self):
        #Step 1: Calc p
        self._calc_p()
        
        #Step 2: Calc phi0
        self._calc_phi0()

        i = 0
        self.improved_phi0 = self.phi0
        while self.phi0 == self.improved_phi0 and i < 1000:
            #Step 3: Calc N0
            self._calc_N0()

            #Step 4: Calc h
            self._calc_h()
            #Step 5: Improved phi0
            self._calc_improved_phi0()
            
            if(self.phi0 != self.improved_phi0):
                self.phi0 = self.improved_phi0
            else:
                break
            i+=1
        print("i", i)
    
    # Help functions for transformation from cartesion to ellipsodial
    def _calc_p(self):
        self.p = math.sqrt(self.estimatedX**2+self.estimatedY**2)

    def _calc_phi0(self):
        self.phi0 = math.atan((self.estimatedZ/self.p)*(self.b**2/self.a**2)**-1)

    def _calc_N0(self):
        self.N0 = (self.a**2)/(math.sqrt(self.a**2*math.cos(self.phi0)**2+self.b**2*math.sin(self.phi0)**2))
    
    def _calc_h(self):
        self.h = (self.p/math.cos(self.phi0))-self.N0

    def _calc_improved_phi0(self):
        self.improved_phi0 = math.atan((self.estimatedZ/self.p)*((self.b**2/self.a**2)*(self.N0/(self.N0+self.h))))


    ###Task 5
    ## Calculate variance-covariance matrix
    def calculate_variance_covariance(self, satellites):
        A = self.A(satellites)
        A_T = A.T
        A_T_A = np.matmul(A_T, A)
        A_T_A_inv = np.linalg.inv(A_T_A)
        qxx = A_T_A_inv[0][0]
        qyy = A_T_A_inv[1][1]
        qzz = A_T_A_inv[2][2]

        table = PrettyTable()
        table.title = "Variance-covariance matrix"
        table.header = False
        table.add_row(A_T_A_inv[0])
        table.add_row(A_T_A_inv[1])
        table.add_row(A_T_A_inv[2])
        table.add_row(A_T_A_inv[3])
        table.align = "l"

        return [qxx, qyy, qzz], table
    
    ## Calculate PDOP
    def calculate_PDOP(self, q_values):
        PDOP = math.sqrt(q_values[0]+q_values[1]+q_values[2])
        return PDOP

