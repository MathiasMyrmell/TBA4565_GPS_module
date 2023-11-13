import numpy as np
from prettytable import PrettyTable
import math

class Receiver:
    def __init__(self, name, phi_0, lambda_0, h_0, known):
        self.name = name
        self.known = known
        #Ellipsiodal coordinates
        self.phi = math.radians(phi_0)
        self.lambda_ = math.radians(lambda_0)
        self.h = h_0

        #Approximate coordinates 
        self.approx_X = None
        self.approx_Y = None
        self.approx_Z = None

        #Cartesian coordinates 
        self.estimated_X = None
        self.estimated_Y = None
        self.estimated_Z = None
        
        #Other values
        self.a = 6378137
        self.b = 6356752.3142


    ##Task 1
    def calculate_cartesian_coordinates(self):
        N = self._calculate_N()
        X = (N + self.h) * math.cos(self.phi) * math.cos(self.lambda_)
        Y = (N + self.h) * math.cos(self.phi) * math.sin(self.lambda_)
        Z = (N * (self.b**2 / self.a**2) + self.h) * math.sin(self.phi)

        self.approx_X = X
        self.approx_Y = Y
        self.approx_Z = Z
        
        if(self.known == True):
            self.estimated_X = X
            self.estimated_Y = Y
            self.estimated_Z = Z


    def _calculate_N(self):
        a = self.a
        b = self.b
        phi = self.phi
        N = a**2 / (math.sqrt(a**2 * math.cos(phi)**2 + b**2 * math.sin(phi)**2))
        return N


    ##Task 2

    def calc_p(self, satellit, t):
        #p_j^i(t)
        # x = (satellite_X-receiver_X)^2
        satellite_coordinates = satellit.get_coordinates(self.name, t)
        satellite_X = satellite_coordinates[0]
        satellite_Y = satellite_coordinates[1]
        satellite_Z = satellite_coordinates[2]
        
        return math.sqrt(((satellite_X -self.x)**2) + ((satellite_Y -self.y)**2) + ((satellite_Z -self.z)**2))

    def get_estimated_coordinates_cartesian(self):
        return [self.estimated_X, self.estimated_Y, self.estimated_Z]
    




    ###Task 4
        ## Transformation from cartesion to ellipsodial
    def calculate_estimated_height(self):
        #Step 1: Calc p
        p = self._calc_p()
        #Step 2: Calc phi0
        phi0 = self._calc_phi0(p)
        improved_phi0 = self.phi0
        diff = None
        while True:
            #Step 3: Calc N0
            N0 = self._calc_N0(phi0)
            #Step 4: Calc h
            h = self._calc_h(phi0, N0, p)
            #Step 5: Improved phi0
            improved_phi0 = self._calc_phi0(p, N0, h)
            #Step 6: Compare approximate and improved phi
            diff = abs(phi0-improved_phi0)
            phi0 = improved_phi0
            if(diff==0):
                break

        self.estimatedHeight = self.h

        estimated = PrettyTable()
        estimated.field_names = (["phi","lambda", "h","T"])
        estimated.add_row([self.radians_to_degree(self.phi0), self.radians_to_degree(math.atan(self.approxY/self.approxX)), self.estimatedHeight, self.receiver_clock_correction])
        estimated.align = "l"
        estimated.title = "Estimated receiver position, geodetic coordinates"
        print(estimated)


    # Help functions for transformation from cartesion to ellipsodial
    def _calc_p(self):
        p = math.sqrt(self.approxX**2+self.approxY**2)
        self.p = p
        return p

    def _calc_phi0(self, p, N0 = 1, h = 0):
        Z = self.approxZ
        e = (self.a**2-self.b**2)/self.a**2
        NNH = N0/(N0+h)
        phi0 = math.atan((Z/p)*(1-e*NNH)**-1)
        self.phi0 = phi0
        return phi0
        
    def _calc_N0(self, phi0):
        N0 = (self.a**2)/(math.sqrt(self.a**2*(math.cos(phi0)**2)+self.b**2*(math.sin(phi0)**2)))
        self.N0 = N0
        return N0
    
    def _calc_h(self, phi0, N0, p):
        h = (p/math.cos(phi0))-N0
        self.h = h
        return h


    


