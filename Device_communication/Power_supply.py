import pyvisa
import sys
import time
class Power_Supply():
    
    def __init__(self,addr) :
        self.address = addr
        resources = pyvisa.ResourceManager()

        self.rs = resources.open_resource(addr)
        self.rs.timeout(10000)

    def set_voltage(self,V):
        try:
            # Setting Voltage
            if V<0:
                raise Exception("Voltage must lager than 0")
            if V>750:
                raise Exception("Voltage can not larger than 750")
            else:
                voltage = "VOLT" + " " + +str(V)
                self.rs.write(voltage)
                time.sleep(1)
        except Exception as e:
            print("Wrong Voltage",e)

    def set_current(self,A):
        try:
            # Setting Voltage
            if A<0:
                raise Exception("Current must lager than 0")
            if A>60:
                raise Exception("Current can not larger than 60")
            else:
                Current = "CURR" + " " + str(A)
                self.rs.write(Current)
                time.sleep(1)
        except Exception as e:
            print("Wrong Currnet",e)


    def Power_setting(self, V, A, M="FIXed"):
        
        # Setting Mode  FIXed / STEP
        try:
            if (M != "FIXed") or (M != "STEP"):
                Exception("it is not included Order")
            else:
                Mode = "MODE" + " "+ M
                self.rs.write(Mode)
                time.sleep(1)

        except Exception as e:
            print("Wrong_Mode",e)


    def Power_on(self):
        self.rs.write("OUTput 1")
        time.sleep(1)
    
    def Power_of(self):
        self.rs.write("OUtput 0")
        time.sleep(1)



    