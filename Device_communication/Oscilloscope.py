import pyvisa

class Oscilo():
    def __init__(self,addr) :
        self.address = addr
        resources = pyvisa.ResourceManager()

        self.rs = resources.open_resource(addr)
        self.rs.write_termination ='\n'
        self.rs.read_termination = '\n'

    def set_file_name(self, name):
