import math
import numpy as np
from prettytable import PrettyTable


class Satellite:
    def __init__(self, name, year, month, day, hour, minute, second, SVClockBias, SVClockDrift, SVClockDriftRate, IODE, Crs, DeltaN, M0, Cuc, e, Cus, sqrtA, Toe, Cic, Omega0, Cis, i0, Crc, omega, OmegaDot, IDOT, CodesL2, GPSWeek, L2Pflag, SVacc, SVhealth, TGD, IODC, Ttom, P, DT, DION, DTROP):
        ## From file
        # SV/EPOCH/SV CLK
        self.name = name
        self.year = int(year)
        self.month = int(month)
        self.day = int(day)
        self.hour = int(hour)
        self.minute = int(minute)
        self.second = int(second)
        self.SVClockBias = float(SVClockBias)
        self.SVClockDrift = float(SVClockDrift)
        self.SVClockDriftRate = float(SVClockDriftRate)
        # BROADCAST ORBIT - 1
        self.IODE = float(IODE)
        self.Crs = float(Crs)
        self.DeltaN =float( DeltaN)
        self.M0 = float(M0)
        # BROADCAST ORBIT - 2
        self.Cuc = float(Cuc)
        self.e = float(e)
        self.Cus = float(Cus)
        self.sqrtA = float(sqrtA)
        # BROADCAST ORBIT - 3
        self.Toe = float(Toe)
        self.Cic = float(Cic)
        self.Omega0 = float(Omega0) #egentlig lambda 0
        self.Cis = float(Cis)
        # BROADCAST ORBIT - 4
        self.i0 = float(i0)
        self.Crc = float(Crc)
        self.omega = float(omega)
        self.OmegaDot = float(OmegaDot)
        # BROADCAST ORBIT - 5
        self.IDOT = float(IDOT)
        self.CodesL2 = float(CodesL2)
        self.GPSWeek = float(GPSWeek)
        self.L2Pflag = float(L2Pflag)
        # BROADCAST ORBIT - 6
        self.SVacc = float(SVacc)
        self.SVhealth = float(SVhealth)
        self.TGD = float(TGD)
        self.IODC = float(IODC)
        # BROADCAST ORBIT - 7
        self.Ttom = float(Ttom)
        # self.fitInterval = float(fitInterval)
        # self.spare = float(spare)

        # EXTRA
        self.P = float(P)
        self.DT = float(DT)
        self.DION = float(DION)
        self.DTROP = float(DTROP)

        ##Given
        self.t = 558000
        self.GM = 3.986005*(10**14)
        self.omega_e = 7.2921151467*(10**-5)

        ## Calculated
        self.a = self.sqrtA**2
        self.t_k = self.calc_tk()
        self.m_k = self.calc_mk()
        self.E_k = self.calc_Ek()
        self.f_k = self.calc_fk()
        self.r_k = self.calc_rk()
        self.i_k = self.calc_ik()
        self.lambda_k = self.calc_lambda_k()
        self.u_k = self.calc_uk()

        self.x = None
        self.y = None
        self.z = None

        ##Saved values
        saved = []

    def __str__(self):
        table = PrettyTable()
        column_names = ["Variable","Value"]
        table.add_column(column_names[0],["SV/EPOCH/SV","Name","Year","Month","Day","Hour","Minute","Second","SVClockBias","SVClockDrift","SVClockDriftRate","BROADCAST ORBIT - 1","IODE","Crs","DeltaN","M0","BROADCAST ORBIT - 2","Cuc","e","Cus","sqrtA","BROADCAST ORBIT - 3","Toe","Cic","Omega0","Cis","BROADCAST ORBIT - 4","i0","Crc","omega","OmegaDot","BROADCAST ORBIT - 5","IDOT","CodesL2","GPSWeek","L2Pflag","BROADCAST ORBIT - 6","SVacc","SVhealth","TGD","IODC","BROADCAST ORBIT - 7","Ttom", "EXTRA", "P", "DT", "DION", "DTROP"])
        table.add_column(column_names[1],["----------------------",self.name,self.year,self.month,self.day,self.hour,self.minute,self.second,self.SVClockBias,self.SVClockDrift,self.SVClockDriftRate,"----------------------",self.IODE,self.Crs,self.DeltaN,self.M0,"----------------------",self.Cuc,self.e,self.Cus,self.sqrtA,"----------------------",self.Toe,self.Cic,self.Omega0,self.Cis,"----------------------",self.i0,self.Crc,self.omega,self.OmegaDot,"----------------------",self.IDOT,self.CodesL2,self.GPSWeek,self.L2Pflag,"----------------------",self.SVacc,self.SVhealth,self.TGD,self.IODC,"----------------------",self.Ttom, "----------------------", self.P, self.DT, self.DION, self.DTROP])
        table.align = "l"
        return table.get_string(title="Satellite: "+self.name)

    def __repr__(self):
        return self.__str__()
    
    ##Calculations
    def calc_tk(self):
        tk = self.t-self.Toe
        if(tk>302400):
            tk = tk-604800
        elif(tk<-302400):
            tk = tk+604800
        return tk
    
    def calc_mk(self):
        return self.M0 + (math.sqrt(self.GM/self.a**3) + self.DeltaN)*self.t_k

    def calc_Ek(self):
        e_old = self.m_k
        for i in range(1,4):
            e_new = e_old + (self.m_k-e_old+self.e*math.sin(e_old))/(1-self.e*math.cos(e_old))
            e_old = e_new
        return e_old

    def calc_fk(self):
        return 2*math.atan(math.sqrt((1+self.e)/(1-self.e))*math.tan(self.E_k/2))

    def cos2(self,x):
        cos_2 = math.cos(x)**2
        sin_2 = math.sin(x)**2
        return cos_2-sin_2

    def sin2(self, x):
        return 2*math.sin(x)*math.cos(x)

    def calc_uk(self):
        return self.omega + self.f_k + self.Cuc*self.cos2(self.omega+self.f_k)+self.Cus*self.sin2(self.omega+self.f_k)

    def calc_rk(self):
        return self.a*(1-self.e*math.cos(self.E_k))+self.Crc*self.cos2(self.omega+self.f_k)+self.Crs*self.sin2(self.omega+self.f_k)

    def calc_ik(self):
        # return self.i0 + self.IDOT*self.t_k + self.Cic*self.cos2(self.omega+self.f_k)+self.Cis*self.sin2(self.omega+self.f_k)
        return self.i0 + self.IDOT*self.t_k + self.Cic*math.cos(2*(self.omega+self.f_k))+self.Cis*math.sin(2*(self.omega+self.f_k))

    def calc_lambda_k(self):
        return self.Omega0 + (self.OmegaDot-self.omega_e)*self.t_k-self.omega_e*self.Toe

    ##Matrices
    def r1(self, x):
        array = np.array([
            [1,0,0],
            [0,math.cos(x),-math.sin(x)],
            [0,math.sin(x),math.cos(x)],
        ])
        # print("r1: \n", array)
        return array
    
    def r2(self, x):
        array = np.array([
            [math.cos(x),0,math.sin(x)],
            [0,1,0],
            [-math.sin(x),0,math.cos(x)],
        ])
        # print("r2: \n", array)
        return array
    
    def r3(self, x):
        array = np.array([
            [math.cos(x),-math.sin(x),0],
            [math.sin(x),math.cos(x),0],
            [0,0,1],
        ])
        # print("r3: \n", array)
        return array

    def calculate_coordinates(self):#Funker (trur æ)
        rk = np.array([
            [self.r_k],
            [0],
            [0]
        ])

        # print("-lambda_k: ", -self.lambda_k)
        # print("-i_k: ", -self.i_k)
        # print("-u_k: ", -self.u_k)
        lamk = self.r3(-self.lambda_k)
        ik = self.r1(-self.i_k)
        uk = self.r3(-self.u_k)
        # print("R3(-lamk): \n", lamk)
        # print("\n")
        # print("R1(-ik): \n", ik)
        # print("\n")
        # print("R3(-uk): \n", uk)
        # print("\n")
        # print("rk: \n", rk)
        # print("\n")
        coordinates = np.matmul(np.matmul(np.matmul(lamk,ik),uk),rk)
        self.x = coordinates[0][0]
        self.y = coordinates[1][0]
        self.z = coordinates[2][0]


        table = self.coordinate_formula_table(lamk, ik, uk, rk)
        return [self.x, self.y, self.z]

    def coordinate_formula_table(self, lamk, ik, uk, rk):
        table = PrettyTable()
        coordinates = np.array([
            [self.x],
            [self.y],
            [self.z]
        ])
        table.field_names = (["R3(-lamk)","R1(-ik)", "R3(-uk)", "rk", "Coordinates"])
        table.add_row([lamk,ik, uk, rk, coordinates])
        return table
    
    def set_coorection_terms_to_zero(self):
        self.saved = [self.DeltaN, self.IDOT, self.OmegaDot, self.Cuc, self.Cus, self.Crc, self.Crs, self.Cic, self.Cis]
        self.DeltaN = 0
        self.IDOT = 0
        self.OmegaDot = 0
        self.Cuc = 0
        self.Cus = 0
        self.Crc = 0
        self.Crs = 0
        self.Cic = 0
        self.Cis = 0

        ##Calculate with new values
        self.t_k = self.calc_tk()
        self.m_k = self.calc_mk()
        self.E_k = self.calc_Ek()
        self.f_k = self.calc_fk()
        self.r_k = self.calc_rk()
        self.i_k = self.calc_ik()
        self.lambda_k = self.calc_lambda_k()
        self.u_k = self.calc_uk()

    def reset_correction_terms(self):
        self.DeltaN = self.saved[0]
        self.IDOT = self.saved[1]
        self.OmegaDot = self.saved[2]
        self.Cuc = self.saved[3]
        self.Cus = self.saved[4]
        self.Crc = self.saved[5]
        self.Crs = self.saved[6]
        self.Cic = self.saved[7]
        self.Cis = self.saved[8]
        ##Calculate with new values
        self.t_k = self.calc_tk()
        self.m_k = self.calc_mk()
        self.E_k = self.calc_Ek()
        self.f_k = self.calc_fk()
        self.r_k = self.calc_rk()
        self.i_k = self.calc_ik()
        self.lambda_k = self.calc_lambda_k()
        self.u_k = self.calc_uk()

        self.calculate_coordinates()

if "__main__" == __name__:
    name = "G03"
    year = 2023
    month = 7
    day = 22
    hour = 12
    minute = 0
    second = 0
    SVClockBias = -0.0002218186855316
    SVClockDrift = 1.648459146963e-11
    SVClockDriftRate = 0.0
    IODE = 21.0
    Crs = -12.75
    DeltaN = 4.415183910132e-09
    M0 = -0.4148605833742
    Cuc = -6.780028343201e-07
    e = 0.00505669566337
    Cus = 3.432855010033e-06
    sqrtA = 5153.62795639
    Toe = 561600.0
    Cic = -4.470348358154e-08
    Omega0 = 2.420706131693
    Cis = 4.470348358154e-08
    i0 = 0.9798129497365
    Crc = 320.5
    omega = 1.070620602228
    OmegaDot = -8.298202796312e-09
    IDOT = -1.157191058795e-10
    CodesL2 = 0.0
    GPSWeek = 2271.0
    L2Pflag = 0.0
    SVacc = 2.0
    SVhealth = 0.0
    TGD = 1.862645149231e-09
    IODC = 21
    Ttom = 554340
    fitInterval = 0
    spare = 0
    s = Satellite(name, year, month, day, hour, minute, second, SVClockBias, SVClockDrift, SVClockDriftRate, IODE, Crs, DeltaN, M0, Cuc, e, Cus, sqrtA, Toe, Cic, Omega0, Cis, i0, Crc, omega, OmegaDot, IDOT, CodesL2, GPSWeek, L2Pflag, SVacc, SVhealth, TGD, IODC, Ttom, fitInterval, spare)

    print("t_k: ", s.t_k)
    print("m_k: ", s.m_k)
    print("E_k: ", s.E_k)
    print("f_k: ", s.f_k)
    print("u_k: ", s.u_k)
    print("r_k: ", s.r_k)
    print("i_k: ", s.i_k)
    print("lambda_k: ", s.lambda_k)

    calculated_coordinates, table = s.calculate_coordinates()
    print("table: \n", table)

    s.set_coorection_terms_to_zero()
    calculated_coordinates, table = s.calculate_coordinates()
    print("table: \n", table)

