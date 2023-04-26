import matplotlib.pyplot as plt
import csv
## Parameter_setting

## Base Unit micro
Input_time = 10 # fixed
Charging_time = list(range(10,65,5)) # Maximu Value = 60*10^(-6)
Waiting_time = 50 # Fixed
Second_pulse_Time = [20] # Maximu Value = 20*10^(-6)
End_TIme = 50 # Fixed

# Base Unit Mega
Sample_rate=10

## Calcultor
Sample_time=1/Sample_rate

for ch in Charging_time:
    for Sp in Second_pulse_Time:
        Input = ["0"] * int(Input_time/Sample_time)
        Charging = ["5"] * int(ch/Sample_time)
        Waiting = ["0"] * int(Waiting_time/Sample_time)
        Second_pulse = ["5"] * int(Sp/Sample_time)
        End = ["0"] * int(End_TIme/Sample_time)

        Result= Input + Charging + Waiting + Second_pulse + End


        # with open('Function_Generator/Charging_{0}_Double_pulse_Test_Time_{1}.csv'.format(ch,Sp), 'w', newline='') as f:
        #     writer = csv.writer(f)
        #     writer.writerow(Result)
        with open('Function_Generator/Charging_{0}_Double_pulse_Test_Time_{1}.csv'.format(ch,Sp), 'w', newline='') as f:
            writer = csv.writer(f, delimiter="\n")
            for item in Result:
                writer.writerow(item)

