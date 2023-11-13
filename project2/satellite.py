
from prettytable import PrettyTable

class Satellite:
    def __init__(self,name, XA1, YA1, ZA1, CPOA1, XB1, YB1, ZB1, CPOB1, XA2, YA2, ZA2, CPOA2, XB2, YB2, ZB2, CPOB2):
        self.name = name
        self.t1 = 172800
        self.t2 = 175020
        ## Time t1: 172800
        # Seen from A
        self.XA1 = float(XA1)
        self.YA1 = float(YA1)
        self.ZA1 = float(ZA1)
        self.CPOA1 = float(CPOA1)
        # Seen from B
        self.XB1 = float(XB1)
        self.YB1 = float(YB1)
        self.ZB1 = float(ZB1)
        self.CPOB1 = float(CPOB1)

        ## Time t2: 175020
        # Seen from A
        self.XA2 = float(XA2)
        self.YA2 = float(YA2)
        self.ZA2 = float(ZA2)
        self.CPOA2 = float(CPOA2)
        # Seen from B
        self.XB2 = float(XB2)
        self.YB2 = float(YB2)
        self.ZB2 = float(ZB2)
        self.CPOB2 = float(CPOB2)

    def __str__(self):
        table = PrettyTable()
        column_names = ["Variable","Value"]
        table.add_column(column_names[0],["T1", "xA", "yA", "zA", "CPOA", "xB", "yB", "zB", "CPOB", "t2", "xA", "yA", "zA", "CPOA", "xB", "yB", "zB", "CPOB"])
        table.add_column(column_names[1],[self.t1, self.XA1, self.YA1, self.ZA1, self.CPOA1, self.XB1, self.YB1, self.ZB1, self.CPOB1, self.t2, self.XA2, self.YA2, self.ZA2, self.CPOA2, self.XB2, self.YB2, self.ZB2, self.CPOB2])
        table.align = "r"
        return table.get_string(title="Satellite: "+self.name)

    def __repr__(self):
        return self.__str__()


    def get_coordinates(self, receiver, time):
        #T1
        if(time == self.t1):
            if(receiver == "A"):
                return [self.XA1, self.YA1, self.ZA1]
            elif(receiver == "B"):
                return [self.XB1, self.YB1, self.ZB1]
        #T2
        elif(time == self.t2):
            if(receiver == "A"):
                return [self.XA2, self.YA2, self.ZA2]
            elif(receiver == "B"):
                return [self.XB2, self.YB2, self.ZB2]
        else:
            return None
        
    def get_carrier_phase_observation(self, receiver, time):
        #T1
        if(time == self.t1):
            if(receiver == "A"):
                return self.CPOA1
            elif(receiver == "B"):
                return self.CPOB1
        #T2
        elif(time == self.t2):
            if(receiver == "A"):
                return self.CPOA2
            elif(receiver == "B"):
                return self.CPOB2
        else:
            return None
        
     
        

    