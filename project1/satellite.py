import math
import numpy as np
from prettytable import PrettyTable
class Satellite:
    def __init__(self, data):
        # SV/EPOCH/SV CLK
        self.name = data[0]
        self.year = int(data[1])
        self.month = int(data[2])
        self.day = int(data[3])
        self.hour = int(data[4])
        self.minute = int(data[5])
        self.second = int(data[6])
        self.SVClockBias = float(data[7])
        self.SVClockDrift = float(data[8])
        self.SVClockDriftRate = float(data[9])
        # BROADCAST ORBIT - 1
        self.IODE = float(data[10])
        self.Crs = float(data[11])
        self.DeltaN =float(data[12])
        self.M0 = float(data[13])
        # BROADCAST ORBIT - 2
        self.Cuc = float(data[14])
        self.e = float(data[15])
        self.Cus = float(data[16])
        self.sqrtA = float(data[17])
        # BROADCAST ORBIT - 3
        self.Toe = float(data[18])
        self.Cic = float(data[19])
        self.Omega0 = float(data[20])
        self.Cis = float(data[21])
        # BROADCAST ORBIT - 4
        self.i0 = float(data[22])
        self.Crc = float(data[23])
        self.omega = float(data[24])
        self.OmegaDot = float(data[25])
        # BROADCAST ORBIT - 5
        self.IDOT = float(data[26])
        self.CodesL2 = float(data[27])
        self.GPSWeek = float(data[28])
        self.L2Pflag = float(data[29])
        # BROADCAST ORBIT - 6
        self.SVacc = float(data[30])
        self.SVhealth = float(data[31])
        self.TGD =float(data[32])
        self.IODC = float(data[33])
        # BROADCAST ORBIT - 7
        self.Ttom = float(data[34])

        # EXTRA
        self.P = float(data[35])#float(P)
        self.DT = float(data[36])#float(DT)
        self.DION = float(data[37])#float(DION)
        self.DTROP = float(data[38])#float(DTROP)

        ##Given
        self.GM = 3.986005*(10**14)
        self.omega_e = 7.2921151467*(10**-5)
        self.c = 299792458
        self.t = 558000#self.calc_transmission_time()

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
        self.saved = []
    def __str__(self):
        table = PrettyTable()
        column_names = ["Variable","Value"]
        table.add_column(column_names[0],["SV/EPOCH/SV","Name","Year","Month","Day","Hour","Minute","Second","SVClockBias","SVClockDrift","SVClockDriftRate","BROADCAST ORBIT - 1","IODE","Crs","DeltaN","M0","BROADCAST ORBIT - 2","Cuc","e","Cus","sqrtA","BROADCAST ORBIT - 3","Toe","Cic","Omega0","Cis","BROADCAST ORBIT - 4","i0","Crc","omega","OmegaDot","BROADCAST ORBIT - 5","IDOT","CodesL2","GPSWeek","L2Pflag","BROADCAST ORBIT - 6","SVacc","SVhealth","TGD","IODC","BROADCAST ORBIT - 7","Ttom", "EXTRA", "P", "DT", "DION", "DTROP"])
        table.add_column(column_names[1],["----------------------",self.name,self.year,self.month,self.day,self.hour,self.minute,self.second,self.SVClockBias,self.SVClockDrift,self.SVClockDriftRate,"----------------------",self.IODE,self.Crs,self.DeltaN,self.M0,"----------------------",self.Cuc,self.e,self.Cus,self.sqrtA,"----------------------",self.Toe,self.Cic,self.Omega0,self.Cis,"----------------------",self.i0,self.Crc,self.omega,self.OmegaDot,"----------------------",self.IDOT,self.CodesL2,self.GPSWeek,self.L2Pflag,"----------------------",self.SVacc,self.SVhealth,self.TGD,self.IODC,"----------------------",self.Ttom, "----------------------", self.P, self.DT, self.DION, self.DTROP])
        table.align = "l"
        return table.get_string(title="Satellite: "+self.name)

    def __repr__(self):
        return self.__str__()
    
    def calc_transmission_time(self):
        transmission_time = 558000- (self.P/self.c) + self.DT
        return transmission_time

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
        return array
    
    def r2(self, x):
        array = np.array([
            [math.cos(x),0,math.sin(x)],
            [0,1,0],
            [-math.sin(x),0,math.cos(x)],
        ])
        return array
    
    def r3(self, x):
        array = np.array([
            [math.cos(x),-math.sin(x),0],
            [math.sin(x),math.cos(x),0],
            [0,0,1],
        ])
        return array

    def calculate_coordinates(self):
        rk = np.array([
            [self.r_k],
            [0],
            [0]
        ])
        lamk = self.r3(self.lambda_k)
        ik = self.r1(self.i_k)
        uk = self.r3(self.u_k)#
        coordinates = np.matmul(np.matmul(np.matmul(lamk,ik),uk),rk)
        self.x = coordinates[0][0]
        self.y = coordinates[1][0]
        self.z = coordinates[2][0]
        return [self.x, self.y, self.z]
   
    ##Change correction terms
    #From original to zero
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

    #From zero to original
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