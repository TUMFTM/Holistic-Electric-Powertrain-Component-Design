# File based on Transmissionpy.py in Gloabl Drive Project
import os
import re
import argparse
import time

def transmission(TP,T,wheel_torque,wheel_rpm):
    with open('./models/Transmission/WTplus/WTplus_2.4.2/work/Beispiele/TP_final/TP_final.wte', 'r', encoding='utf-8') as f:
        fileContents = f.readlines()
        
    if TP == 1:
        fileContents[48] = 'TEMP_EINSPRITZ = %f' % T + '\n'
        fileContents[84] = 'DREHZAHL   = %f' % wheel_rpm + '\n'
        fileContents[90] = 'UK = 180  DA = 62    TUA = %f' % (-wheel_torque) + '\n'
    elif TP == 2:
        fileContents[48] = 'TEMP_EINSPRITZ = %f' % T + '\n'
        fileContents[84] = 'DREHZAHL   = %f' % wheel_rpm + '\n'
        fileContents[90] = 'UK = 180  DA = 62    TUA = %f' % (-wheel_torque) + '\n'
    elif TP == 3:
        fileContents[37] = 'TEMP_EINSPRITZ = %f' % T + '\n'
        fileContents[111] = 'DREHZAHL   = %f' % wheel_rpm + '\n'
        fileContents[117] = 'UK = 180  DA = 62    TUA = %f' % (-wheel_torque) + '\n'
    elif TP == 4:
        fileContents[36] = 'TEMP_EINSPRITZ = %f' % T + '\n'
        fileContents[69] = 'DREHZAHL   = %f' % wheel_rpm + '\n'
        fileContents[74] = 'UK = 20  DA = 62    TUA = %f' % (-wheel_torque) + '\n'
    elif TP == 5:
        fileContents[46] = 'TEMP_EINSPRITZ = %f' % T + '\n'
        fileContents[76] = 'DREHZAHL   = %f' % wheel_rpm + '\n'
        fileContents[81] = 'UK = 180  DA = 59.07    TUA = %f' % (-wheel_torque) + '\n'

    with open('./models/Transmission/WTplus/WTplus_2.4.2/work/Beispiele/TP_final/TP_final.wte', 'w', encoding='utf-8') as file:
        file.writelines(fileContents)

    os.chdir('./models/Transmission/WTplus/WTplus_2.4.2/bin')
    os.system('WTplus.exe > NUL');
    os.chdir('../../../../../')
    
    with open('./models/Transmission/WTplus/WTplus_2.4.2/work/Beispiele/TP_final/TP_final.wta', 'rb') as f:
        fileContents = f.readlines()
        
    fileContents = [i.decode('ascii', errors='replace') for i in fileContents]

    for i in range(len(fileContents)):
        if 'Welle  1:' in fileContents[i]:
            gearboxrpmStr = fileContents[i+2]
            gearbox_rpm = float(re.findall(r'-?\d+\.?\d*', gearboxrpmStr)[0])
            gearboxtorqueStr = fileContents[i+3]
            gearbox_torque = float(re.findall(r'-?\d+\.?\d*', gearboxtorqueStr)[0])
            break

    for i in range(len(fileContents)):
        if 'Wirkungsgrad    .' in fileContents[i]:
            etaStr = str(fileContents[i])
            eta = re.findall(r'-?\d+\.?\d*', etaStr)
            eta = float(eta[0])
            break

    for i in range(len(fileContents)):
        if 'Gesamtverlustleistung' in fileContents[i]:
            powerlossStr = str(fileContents[i])
            powerloss = re.findall(r'-?\d+\.?\d*', powerlossStr)
            powerloss = float(powerloss[0])
            break

        
    return gearbox_rpm, gearbox_torque, powerloss, eta


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Compute transmission")
    parser.add_argument("--tp", type=int)
    parser.add_argument("--temp", type=float)
    parser.add_argument("--torque", type=float)
    parser.add_argument("--rpm", type=float)
    opt = parser.parse_args()
    gr, gt, p, et = transmission(opt.tp, opt.temp, opt.torque, opt.rpm)
    # print("gearbox rpm:", gr)
    # print("gearbox torque:", gt)
    # print("power loss:", p)
    with open('./transmission_output.txt', 'w') as f:
        #f.write(f"{gr}\n")
        #f.write(f"{gt}\n")
        f.write(f"{p}\n")
        #f.write(f"{et}\n")

