import pyvisa
import time

class Function_Generator():
    
    def __init__(self,addr) :
        self.address = addr
        resources = pyvisa.ResourceManager()

        self.rs = resources.open_resource(addr)
        self.rs.timeout(10000)


    def arb(self, V, fun_path, Sample, V_off=0):
        
        self.rs.write("FUNCtion ARB")
        # Function select
        Fun= "FUNCtion:ARBitrary" + " " + fun_path
        self.re.write(Fun)
        time.sleep(1)

        # Setting Voltage
        voltage = "VOLTage" + " " + "+"+str(V)
        self.rs.write(voltage)
        time.sleep(1)

        # Setting Offset
        offset = "VOLTage:OFFSet" + " " + "+"+str(V_off)
        self.rs.write(offset)
        time.sleep(1)
        
        # Sample Rate
        Sample_rate = "FUNC:ARB:SRAT" + " " + Sample

        self.rs.write(Sample_rate)
        time.sleep(1)

    def run_trigger(self,ch_num):
        trigger = "TRIGger" + str(ch_num)
        self.rs.write(trigger)
        time.sleep(1)

    def complete(self):
        self.rs.write("OPC?")
