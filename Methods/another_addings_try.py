import wlmData
import wlmConst
import ctypes
import sys
import WLM_methods
import pyvisa
import ThorlabsPM100
import matplotlib.pyplot as plt
import time
from multiprocessing import Process, Lock
#########################################################
# Set the DLL_PATH variable according to your environment
#########################################################
DLL_PATH = "wlmData.dll"

def resonator_info_maker():
    mode = True
    down_reference = 580.038905
    d_freq = wlmData.dll.ConvertUnit(down_reference, wlmConst.cReturnWavelengthVac,
                                     wlmConst.cReturnFrequency)
    upper_reference = 580.03750
    u_freq = wlmData.dll.ConvertUnit(upper_reference, wlmConst.cReturnWavelengthVac,
                                     wlmConst.cReturnFrequency)
    stabilisation_time = 10
    # negative to make higher frequency
    PID_step_mV = -0.25
    time_limit = 240
    channel_to_measure = 2

    d = WLM_methods.stepping_PID_course(mode, d_freq, u_freq, stabilisation_time,
                                        PID_step_mV, time_limit, channel_to_measure)

    x_absc = [d[0][i][0] for i in range(d[1])]
    y_ord = [d[0][i][1] for i in range(d[1])]
    x = x_absc.copy()
    y = y_ord.copy()

    cort = WLM_methods.find_max_mod(x,y)
    max_mod_power = cort[0]
    max_mod_freq = cort[1]
    max_mod_index = cort[2]
    max_mod_breadth = WLM_methods.find_breadth_mod(x, y, max_mod_index)

    WLM_methods.del_mod(x, y, max_mod_index, max_mod_breadth)

    cort2 = WLM_methods.find_max_mod(x,y)
    max_mod_power2 = cort2[0]
    max_mod_freq2 = cort2[1]
    max_mod_index2 = cort2[2]
    max_mod_breadth2 = WLM_methods.find_breadth_mod(x, y, max_mod_index2)

    Delta_mod = max_mod_freq - max_mod_freq2
    title = "1st mod breadth: " + str(max_mod_breadth) + "\n2nd mod breadth: " + str(max_mod_breadth2) + "\nDelta mod: " + str(Delta_mod)

    # plot find_k returning list
    fig, ax = plt.subplots(figsize=(5, 3), layout='constrained')
    ax.set_xlabel('Delta frequency [THz]')
    ax.set_ylabel('Power [W]')
    ax.set_title('Power frequency dependency')
    ax.grid()
    ax.plot(x_absc, y_ord, label=title)
    plt.ylim(min(y_ord), max(y_ord)*1.05)
    plt.show()

# Load DLL from DLL_PATH
try:
    wlmData.LoadDLL(DLL_PATH)
except:
    sys.exit("Error: Couldn't find DLL on path %s. Please check the DLL_PATH variable!" % DLL_PATH)

# Checks the number of WLM server instance(s)
if wlmData.dll.GetWLMCount(0) == 0:
    print("There is no running wlmServer instance(s).")
else:
    # Read Type, Version, Revision and Build number
    Version_type = wlmData.dll.GetWLMVersion(0)
    Version_ver = wlmData.dll.GetWLMVersion(1)
    Version_rev = wlmData.dll.GetWLMVersion(2)
    Version_build = wlmData.dll.GetWLMVersion(3)
    print("WLM Version: [%s.%s.%s.%s]" % (Version_type, Version_ver, Version_rev, Version_build))

    # for Powermeter include from resonator
    # PM400_connect = "USB0::0x1313::0x8075::P5004922::INSTR"
    # PM100_connect = "USB0::0x1313::0x8078::P0008894::INSTR"
    # rm = pyvisa.ResourceManager()
    #
    # inst = rm.open_resource(PM400_connect)
    # power_meter = ThorlabsPM100.ThorlabsPM100(inst=inst)
    # print(power_meter.read)

    # check for channel in charge
    # it's sometimes a problem to set the channel, so it's needed to
    # make setting more than 1 time until set function returns 0

    # port = ctypes.c_long(0)
    # channel = ctypes.c_long(0)
    # active_chnl = wlmData.dll.GetActiveChannel(1, ctypes.byref(port), 0)
    # print(active_chnl, port.value)
    # x = wlmData.dll.SetActiveChannel(3, 1, 2,  0)
    # active_chnl1 = wlmData.dll.GetActiveChannel(3, ctypes.byref(port), 0)
    # print(x, active_chnl1, port.value)
    # print(round(wlmData.dll.GetWavelengthNum(2,0),6))
    # 580.0415 - 4096
    # 580.03530 - 0

    wave = 580.04
    time_pause = wlmData.dll.GetExposureNum(2, 1, 0)
    start_PID = wlmData.dll.GetDeviationSignalNum(1, 0)
    print(WLM_methods.reference_const_PID_stabilisator_with_timing_version_for_test(False, wave,
                                                           WLM_methods.cDependFrequencyPID, 4096,
                                                           time_pause*3, 240, start_PID, 2))

    # up = wlmData.dll.ConvertUnit(580.03530, wlmConst.cReturnWavelengthVac,
    #                                                wlmConst.cReturnFrequency)
    # down = wlmData.dll.ConvertUnit(580.04145, wlmConst.cReturnWavelengthVac,
    #                                                   wlmConst.cReturnFrequency)
    #
    # print((up - down) / 4096)


    # different meanings of koefficient in frequency - PID dependency
    # - is for we use frequency
    k = -1.3960736881714986e-06
    k2 = -3.76914227226224e-08
    k3 = -1.1114092415877479e-05
    k4 = -2.23e-06 # better last
    k5 = -7.598893318173907e-07
    k6 = k4 * 1.5225

    # points = 500
    # reference_wl = 580.038015
    # delta = round((reference_wl - wlmData.dll.GetWavelengthNum(1,0)), 7)
    # print(delta / k)
    # d = WLM_methods.wavelength_PID_bond(points, 0.1, wlmData.dll.GetDeviationSignalNum(1,0), exposition).copy()
    # koef = (d[points-1][1] - d[1][1]) / (d[points-1][0] - d[1][0])
    # print(koef)

    # find k

    # exposition = wlmData.dll.GetExposureNum(1, 1, 0)
    # start = wlmData.dll.GetDeviationSignalNum(1, 0)
    # cur_freq = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
    #                                    wlmConst.cReturnFrequency)
    # wlmData.dll.SetDeviationSignalNum(1, 1850)
    # WLM_methods.time_counter(cur_freq)
    # print(WLM_methods.find_k2(1000, 0.125, 1850))
    # wlmData.dll.SetDeviationSignalNum(1, start)

    # ref = 580.039705
    # mode = False
    # timer = 50
    # exposition = wlmData.dll.GetExposureNum(1, 1, 0)
    # koef = k4
    # max_dev = 2000
    # max_PID_val = 4096
    # start_PID_point = wlmData.dll.GetDeviationSignalNum(1, 0)
    # WLM_methods.reference_const_PID_stabilisator(mode, ref, koef, max_PID_val,
    #                                              wlmData.dll.GetExposureNum(1, 1, 0)*1.2, timer, start_PID_point)

    # plots 2 graphs near each other. good to see the difference
    # fig, (ax1, ax2) = plt.subplots(1,2,figsize=(5, 3), layout='constrained')
    # ax1.set_xlabel('Delta frequency [THz]')
    # ax1.set_ylabel('Power [W]')
    # ax1.set_title('Power frequency dependency')
    # ax1.grid()
    # ax1.plot(x, y)
    # ax2.set_xlabel('Delta frequency [THz]')
    # ax2.set_ylabel('Power [W]')
    # ax2.set_title('Power frequency dependency')
    # ax2.grid()
    # ax2.plot(x_absc, y_ord)
    # plt.ylim(min(y_ord), max(y_ord)+0.05*max(y_ord))
    # plt.show()


    # changing of wavelength in 3 digit from 9 to 1 after comma returns about 3184.714925527756 mV of delta PID
    # 398.0845614748144 - 9 to 8
    # 796.1704955501652 - 9 to 7
    # 1194.2578022768425 - 9 to 6
    # 1592.3464816548462 - 9 to 5
    # 1990.4365336333865 - 9 to 4
    # 2388.527958263253 - 9 to 3
    # 2786.620755544446 - 9 to 2
    # changes that we make in frequency are nearly the same as we've made in wavelength




