import ctypes

import numpy as np

import wlmData
import wlmConst
import time
import matplotlib.pyplot as plt
import pyvisa
import ThorlabsPM100
import math
def set_wl(amount):
    '''
    The function can set PID_course of every type, not only constant
    we can place here some law, e.g. sin or cos but there is some style of doing
    this that is said in manual to WLM
    :param amount: the set of PID course
    :return:
    '''
    new_b = bytes('= ' + str(amount), encoding='utf-8')
    new_PIDC = ctypes.create_string_buffer(new_b)
    print(new_PIDC.value.decode())
    if (wlmData.dll.SetPIDCourseNum(1, new_PIDC) == wlmConst.ResERR_NoErr):
        print("Successful write-in %s" % new_PIDC.value.decode())

def set_wl2(amount):
    '''
    The same as set_wl but without feedback
    :param amount: the set of PID course
    :return:
    '''
    new_b = bytes('= ' + str(amount), encoding='utf-8')
    new_PIDC = ctypes.create_string_buffer(new_b)
    wlmData.dll.SetPIDCourseNum(1, new_PIDC)

def wavelength_regulation(cycle_steps, initial_wave, delta_wave, time_pause_set, time_pause_getpwr, wl_precision):
    ''' Definition:
        _________________________________________________________________
        The function makes needed steps of 'set - measure' wavelength
        in WLM with explicit interruptions and wavelength step setting
        _________________________________________________________________
        Parameters:
        ________________________________________________________________
        cycle_steps - how many cycles of 'set - measure' do we need
        initial_wave - current wavelength, set in WLM
        delta_wave - step of wavelength to change WLM initial wavelength
        time_pause_set - time to wait after setting wavelength of WLM in ms
        time_pause_getpwr - time to wait after getting power of measurement shot in WLM of al-m in ms
        wl_precision - precision of wavelength in number of digits
        _________________________________________________________________
        Return:
        _________________________________________________________________
        A dictionary of tuples that consist of 3 parts.
        1 - time between set_wl and the moment after last pause including all pauses
        2 - wavelength that was set
        3 - wavelength that was got after set to check whether it was set right
        4 - power of the measurement shot
    '''
    d = {}
    i = 1
    while(i <= cycle_steps):
        wl_set = round(initial_wave + i*delta_wave, wl_precision)
        time1 = time.time()
        set_wl2(wl_set)
        time.sleep(time_pause_set*0.001)
        wavelength1 = round(wlmData.dll.GetWavelengthNum(1, 0), wl_precision)
        power1 = round(wlmData.dll.GetPowerNum(1,0),2)
        time.sleep(time_pause_getpwr * 0.001)
        time2 = time.time()
        d[i] = (round((time2-time1)*1000,2), wl_set, wavelength1, power1)
        i+=1
    return d

def wavelength_PID_bond(points, PID_step, PID_start, expo_time):
    i = 0
    d = {}
    while (i < points):
        PID = PID_start + i*PID_step
        wlmData.dll.SetDeviationSignalNum(1, PID)
        time.sleep((expo_time)/1000)
        wave = wlmData.dll.GetWavelengthNum(1, 0)
        d[i] = (PID, wave)
        i+=1
    wlmData.dll.SetDeviationSignalNum(1, PID_start)
    return d

def plot_wavelength_PID_bond(points, PID_step):
    expo_time =  wlmData.dll.GetExposureNum(1,1,0)
    start_pid = wlmData.dll.GetDeviationSignalNum(1,0)
    d = wavelength_PID_bond(points, PID_step, start_pid, expo_time).copy()
    x = [round(d[i][0],2) for i in range(points)]
    y = [round(d[i][1],7) for i in range(points)]
    fig, ax = plt.subplots(figsize=(5,3), layout='constrained')
    ax.plot(x,y)
    plt.ylim(min(y), max(y))
    plt.show()

def wl_stabilisation_after_set(wave, time_pause, precision):
    time1 = time.time()
    set_wl2(wave)
    delta = round(abs(wave - wlmData.dll.GetWavelengthNum(1,0)),6)
    print(delta)
    while(delta > precision):
        time.sleep(time_pause*0.001)
        delta = round(abs(wave - wlmData.dll.GetWavelengthNum(1,0)),6)
    time2 = time.time()
    return round((time2-time1)*1000,2)

def reference_const_PID_stabilisator(mode: bool, reference_wl, koef, max_PID_val, time_pause, timer, start_PID_point = 1860):
    # if bool mode true than frequency other way - wavelength
    '''
        The function stabilises the reference value of frequency for "timer" seconds;
        2nd version of algorithm

        :param mode: True - reference_wl should be frequency
        :param reference_wl: frequency or wavelength (no matter in fact, unless you didn't set the right mode. We use frequencies inside algorithm anyway)
        :param koef: koef of dependency between PID mV and frequency
        :param max_PID_val: max val in mV for PID
        :param time_pause: pause after setting value needed
        :param timer: time of algorithm execution
        :param start_PID_point: starting point of PID setting
        :return:
    '''
    stabilised = False
    PID_step = 0
    PID_current = start_PID_point
    reference = reference_wl
    if(not mode):
        reference = wlmData.dll.ConvertUnit(reference, wlmConst.cReturnWavelengthVac,
                                wlmConst.cReturnFrequency)
    time1 = time.time()
    while(True):
        if(not stabilised):
            PID_current = PID_current + PID_step
            wlmData.dll.SetDeviationSignalNum(1, PID_current)
            time.sleep(time_pause / 1000)
        wave_current = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                       wlmConst.cReturnFrequency)
        delta = reference - wave_current
        abs_dev = abs(delta)
        if(abs_dev <= 10e-07):
            if (delta > 0):
                PID_step = -0.125
            elif (delta < 0):
                PID_step = 0.125
        elif(abs_dev > 10e-07 and abs_dev <= 10e-06):
            if (delta > 0):
                PID_step = -1.25
            elif (delta < 0):
                PID_step = 1.25
        elif(abs_dev > 10e-06 ):
            PID_step = delta / koef
        if(( PID_current + PID_step ) > max_PID_val or ( PID_current + PID_step ) < 0):
            print("Out of bounds PID setting")
            break
        if(round(delta,8) == 0.00000000):
            stabilised = True
        elif(stabilised):
            stabilised = False
        if(time.time() - time1 > timer):
            print("Time limit exceeded")
            break

def wl_stabilisation_through_PID_const(reference_wl, koef, max_dev, time_pause, timer, start_PID_point = 1860):
    '''
    Function stabilises the wavelength on reference_val
    1st version of algorithm
    :param reference_wl: reference value in nm
    :param koef: koefficient of bond between PID in mV and wavelength as frequency
    :param max_dev: max deviation in mV
    :param time_pause: pause after setting PID inside alg-m
    :param timer: lifetime of function in sec
    :param start_PID_point: start point for al-m to change PID
    :return:
    '''
    flag = False
    flag2 = True # regulates PID_STEP if laser swims out of reference
    first_stabilisation_flag = True # regulates PID_STEP if laser isn't stabilised yet
    counter = 0
    PID_step = 0
    time2 = 0
    i = 0
    k = koef
    PID_prev = 0
    reference_wl_Thz = wlmData.dll.ConvertUnit(reference_wl, wlmConst.cReturnWavelengthVac,
                                wlmConst.cReturnFrequency)
    wlmData.dll.SetDeviationSignalNum(1, start_PID_point)
    time.sleep(time_pause/1000)
    time1 = time.time()
    wave = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                   wlmConst.cReturnFrequency)
    while(True):
        current = wlmData.dll.GetDeviationSignalNum(1,0)
        delta = reference_wl_Thz - wave
        if(flag2):
            if (delta > 0):
                PID_step = -0.125
            elif (delta < 0):
                PID_step = 0.125
        if(not flag and first_stabilisation_flag):
            PID_step = delta / k
        if(abs(PID_step - PID_prev) > max_dev):
            print("!")
            if(delta > 0):
                PID_step = -100
            elif(delta < 0):
                PID_step = 100
        if(round(delta,7) == 0.00000000):
            if(not flag):
                print("stabilized!")
                flag = True
                flag2 = False
                first_stabilisation_flag = False
        elif(flag):
            flag = False
            flag2 = True
        if(not flag):
            wlmData.dll.SetDeviationSignalNum(1, current + PID_step)
        time.sleep(time_pause/1000)
        wave = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                       wlmConst.cReturnFrequency)
        PID_prev = PID_step
        time2 = time.time()
        if((time2-time1)>timer):
            break

def find_k(points, PID_step, PID_start, time_pause):
    '''
      wave - before set
      wave2 - after set
      :param points:
      :param PID_step:
      :param PID_start:
      :param time_pause:
      :param delta_time: wait time until we get another wavelength is our step time_pause was false
      :return:
      '''
    assert points > 2, "Error: 2 or more points needed"
    i = 0
    d = {}
    while (i < points):
        PID = PID_start + i * PID_step
        # wave = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
        #                                wlmConst.cReturnFrequency)
        wlmData.dll.SetDeviationSignalNum(1, PID)
        time.sleep(time_pause/1000)
        wave2 = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                       wlmConst.cReturnFrequency)
        d[i] = (PID, wave2)
        i += 1
    wlmData.dll.SetDeviationSignalNum(1, PID_start)
    k = 1
    koefs_list = []
    while(k < points):
        koefs_list.append((d[k][1]-d[k-1][1])/PID_step)
        k+=1
    return np.average(koefs_list)

def find_k2(points, PID_step, PID_start):
    '''
      wave - before set
      wave2 - after set
      :param points:
      :param PID_step:
      :param PID_start:
      :param time_pause:
      :param delta_time: wait time until we get another wavelength is our step time_pause was false
      :return:
      '''
    assert points > 2, "Error: 2 or more points needed"
    assert (points%2) == 0, "Error: even points number needed"
    i = 0
    d = {}
    while (i < points):
        PID = PID_start + i * PID_step
        wave = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                       wlmConst.cReturnFrequency)
        wlmData.dll.SetDeviationSignalNum(1, PID)
        d[i] = (PID, time_counter(wave))
        i += 1
    wlmData.dll.SetDeviationSignalNum(1, PID_start)
    k = 0
    koefs_list = []
    index_max = points-1
    while(k < points/2):
        koefs_list.append((d[index_max - k][1]-d[k][1])/(d[index_max - k][0]-d[k][0]))
        k+=1
    return np.average(koefs_list)

def plot_wavelength_PID_find_k(points, PID_step):
    '''
    To use the function you need to return dictionary d in find_k function
    It will show you the graph of dependency between PID in mV and frequency measured in THz
    :param points: number of different measurements
    :param PID_step: step of PID in mV to use in find_k / or find_k2.
    :param time_pause: time to pause after setting the current PID in mV
    Remark: highly recommended to set delay more than exposition time to get an adequate result
    :return:
    '''
    start_pid = wlmData.dll.GetDeviationSignalNum(1,0)
    d = find_k2(points, PID_step, start_pid).copy()
    x = [round(d[i][0],2) for i in range(points)]
    y = [round(d[i][1],7) for i in range(points)]
    fig, ax = plt.subplots(figsize=(5,3), layout='constrained')
    ax.scatter(x,y)
    plt.ylim(min(y), max(y))
    plt.show()

def time_counter(ref_frequency):
    '''
    Counts time until we get new measure
    We can use this function to get new frequency of measurement with more accurate time difference
    Remark: we can get the previous frequency even if the delay is exposition time

    :param ref_frequency: it's previous frequency meaning
    :return: new frequency
    '''
    cur_freq = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                       wlmConst.cReturnFrequency)
    while((cur_freq - ref_frequency) == 0.00000000):
        cur_freq = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                        wlmConst.cReturnFrequency)
    return cur_freq

def triangle_PID_course(mode: bool, down_reference, upper_reference, stabilisation_time,
                        PID_step_mV, time_limit):
    assert stabilisation_time < time_limit, "Error: too short time limit"
    koef = -2.23e-06*1.5225
    d_reference = down_reference
    u_reference = upper_reference
    delta_freq = koef*PID_step_mV
    d = {}
    i = 0
    if(not mode):
        d_reference = wlmData.dll.ConvertUnit(d_reference, wlmConst.cReturnWavelengthVac,
                                wlmConst.cReturnFrequency)
        u_reference = wlmData.dll.ConvertUnit(upper_reference, wlmConst.cReturnWavelengthVac,
                                wlmConst.cReturnFrequency)
    max_PID_val = 4096
    time1 = time.time()
    while(True):
        start_PID_point = wlmData.dll.GetDeviationSignalNum(1, 0)
        reference_const_PID_stabilisator(True, d_reference, koef, max_PID_val,
                                                     wlmData.dll.GetExposureNum(1, 1, 0) * 1.2, stabilisation_time, start_PID_point)
        PID_current = wlmData.dll.GetDeviationSignalNum(1, 0)
        while(True):
            PID_current = PID_current + PID_step_mV
            cur_freq = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1,0), wlmConst.cReturnWavelengthVac,
                                               wlmConst.cReturnFrequency)
            if(cur_freq + delta_freq > u_reference):
                break
            wlmData.dll.SetDeviationSignalNum(1, PID_current)
            time_counter(cur_freq)
            if(time.time()-time1 > time_limit):
                break
        if (time.time() - time1 > time_limit):
            break


def initialise_power_meter(connection_row):
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource(connection_row)
    power_meter = ThorlabsPM100.ThorlabsPM100(inst=inst)
    return power_meter

def find_max_mod(absc_list, ord_list):
    max_mod = max(ord_list)
    freq_of_max_mod = 0.
    index = -1
    for k in range(len(absc_list)):
        if (ord_list[k] == max_mod):
            freq_of_max_mod = absc_list[k]
            index = k
            break
    return max_mod,freq_of_max_mod,index

def find_breadth_mod(absc_list, ord_list, index_mod):
    size = len(absc_list)
    assert index_mod >= 0 and index_mod < size, "Error: index of mod is out of bounds. Check lists of values and index"
    k1 = index_mod
    k2 = index_mod
    barier = ord_list[index_mod] / 2
    while(ord_list[k1] >= barier):
        k1+=1
    while(ord_list[k2] >= barier):
        k2-=1
    assert k1 >= 0 and k2 < size, "Error: index is out of bounds. The mod can be cut off"
    return absc_list[k1]-absc_list[k2]
def del_mod(absc_list, ord_list, index_mod, breadth):
    size = len(absc_list)
    assert index_mod >= 0 and index_mod < size, "Error: index of mod is out of bounds. Check lists of values and index"
    k1 = index_mod
    k2 = index_mod
    left = absc_list[index_mod] - breadth/2
    right = absc_list[index_mod] + breadth/2
    while (absc_list[k1] <= right):
        k1+=1
    while (absc_list[k2] >= left):
        k2 -= 1
    for _ in range(k1 - k2 + 1):
        del absc_list[k2]
        del ord_list[k2]
def stepping_PID_course(mode: bool, down_reference, upper_reference, stabilisation_time,
                        PID_step_mV, time_limit):
    assert stabilisation_time < time_limit, "Error: too short time limit"
    koef = -2.23e-06 * 1.5225
    flag = False
    d_reference = down_reference
    u_reference = upper_reference
    delta_freq = koef * PID_step_mV
    d = {}
    i = 0
    if (not mode):
        d_reference = wlmData.dll.ConvertUnit(d_reference, wlmConst.cReturnWavelengthVac,
                                              wlmConst.cReturnFrequency)
        u_reference = wlmData.dll.ConvertUnit(upper_reference, wlmConst.cReturnWavelengthVac,
                                              wlmConst.cReturnFrequency)
    max_PID_val = 4096

    power_meter = initialise_power_meter("USB0::0x1313::0x8078::P0008894::INSTR")

    time1 = time.time()
    start_PID_point = wlmData.dll.GetDeviationSignalNum(1, 0)
    reference_const_PID_stabilisator(True, d_reference, koef, max_PID_val,
                                     wlmData.dll.GetExposureNum(1, 1, 0) * 1.2, stabilisation_time, start_PID_point)
    PID_current = wlmData.dll.GetDeviationSignalNum(1, 0)

    d[i] = (abs(wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                        wlmConst.cReturnFrequency) - d_reference), power_meter.read)
    i += 1
    while (True):
        flag = False
        PID_current = PID_current + PID_step_mV
        cur_freq = wlmData.dll.ConvertUnit(wlmData.dll.GetWavelengthNum(1, 0), wlmConst.cReturnWavelengthVac,
                                           wlmConst.cReturnFrequency)
        if (cur_freq + delta_freq > u_reference):
            percent = (u_reference-cur_freq) / delta_freq
            if( percent > 0.1 and percent <= 1. and abs(PID_step_mV*percent) >= 0.125):
                PID_current = PID_current - PID_step_mV
                PID_current = PID_current + PID_step_mV*percent
                wlmData.dll.SetDeviationSignalNum(1, PID_current)
                freq = time_counter(cur_freq)
                d[i] = (abs(freq - d_reference), power_meter.read)
                if (d[i][0] / d[i - 1][0] > 1000):
                    del d[i]
                    flag = True
                if (not flag):
                    i += 1
            break
        wlmData.dll.SetDeviationSignalNum(1, PID_current)
        freq = time_counter(cur_freq)
        d[i] = (abs(freq - d_reference), power_meter.read)
        if(d[i][0]/d[i-1][0] > 1000):
            del d[i]
            flag = True
        if(not flag):
            i += 1
        if (time.time() - time1 > time_limit):
            break
    return d, i