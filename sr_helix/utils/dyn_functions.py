#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
from datetime import datetime
from math import cos, sin
import numpy as np
import scipy.io
from utils.constants import *

def calc_error(p_d, p_0):
    err = math.sqrt((p_d[0]-p_0[0])**2 + (p_d[1]-p_0[1])**2 + (p_d[2]-p_0[2])**2)
    return err

def to_motor_steps(th0):
    """
    Input: the theta needed (in radians)
    Returns: new motor step command to be inputted into set_goal_position method
    """
    if len(th0) == 1:
        steps = th0 * (4096/(2 *math.pi))
        return int(steps)
    else:
        steps = th0 * (4096/(2 *math.pi))
        return steps.astype(int).tolist()
    
def to_radians(steps):
    """
    Input: takes in dynamixel motor steps pos
    Output: radian value of dynamixel
    """
    rad_val = (steps * ((2 * math.pi)/4096))
    return rad_val

def to_rad_secs(revs):
    """
    Input:  Present velocity of dynamixel in rev/min
    Output: Velocity in rad/sec
    """
    rad_sec = revs * (2 * math.pi)/60
    return rad_sec
def diff(q, q_old, dt):
    """
    Input: current config and old config
    Output: calculates the derivative over some dt
    """
    q_dot = (q - q_old)/dt
    return q_dot

def get_time():
    """
    Output: returns time object you can use to calculate dt
    To read the time stamp, you do "t_obj.time()"
    """
    t = datetime.now().strftime("%H:%M:%S.%f")
    t_obj = datetime.strptime(t, "%H:%M:%S.%f")
    return t_obj

def get_qindex(mod_clock, tvec):
    """
    input: takes in the current time we're getting and makes sure we stay at proper
    time element 
    output: outputs the proper time vector index we need for our pos, vel, and accel vectors
    Notes:
    - recall tvec is a numpy matrix object that is (1,56) dim 
    """
    qindex = 0
    for t in range(1, tvec.shape[1]):
        if mod_clock >= tvec[0, t-1] and mod_clock < tvec[0, t]:
            qindex = t - 1
            return qindex
        elif t==(tvec.shape[1] - 1) and mod_clock >= tvec[0, t]:
            qindex = t
            return qindex
    return qindex

def mat2np(fname, typ):
    """
    Function that converts mat file to numpy matrix
    Parameters
    ----------------------------
        **both inputs are in string format
        fname = name of the file, i.e 'q.mat'
        typ = name of matrix you're trying to pull from dictionary, i.e 'q'
    """
    mat = scipy.io.loadmat(fname)[typ]
    # # for the 56 samples you created grab the gen coords
    return mat

def grab_current(tau_cables, min_torque, max_torque):
    """
    Parameters
    ----------------------------
    tau_cables: the raw currents calculated from controller (in mA?)
    min_torque: the minimum torque required to actuate module
    max_torque: the maximum torque module can be actuated without 
                damaging the module (basially point where plates touch)
    l: current cable lengths
        can do something like if current makes no
    Returns
    -----------------------------
    Ensures we only pass safe currents that the module can physically handle
    """
    mod_input = [0, 0, 0]
    for i in range(len(tau_cables)):
        if tau_cables[i][0].item() < 0:
            # negative current case
            if tau_cables[i][0].item() < -max_torque:
                mod_input[i] = -max_torque
            elif tau_cables[i][0].item() > -min_torque:
                mod_input[i] = 0
            else:
                mod_input[i] = int(tau_cables[i][0].item())
        else:    
            # positive current case
            if tau_cables[i][0].item() > max_torque:
                mod_input[i] = max_torque
            elif tau_cables[i][0].item() < min_torque:
                mod_input[i] = 0
            else:
                mod_input[i] = int(tau_cables[i][0].item())   
    m1 = mod_input[0]
    m2 = mod_input[1]
    m3 = mod_input[2]
    mod_input = [m1, m2, m3]
    # print(f"Our new mod input: {mod_input}\n")
    return mod_input

def grab_arm_current(tau, min_torque, max_torque):
    """
    Parameters
    ----------------------------
    tau_cables: [tau1, tau2, tau3]
    tau1 = [tau for motor1, tau for motor2, tau for motor3]
    min_torque: the minimum torque required to actuate module
    max_torque: the maximum torque module can be actuated without 
                damaging the module (basially point where plates touch)
    Returns
    -----------------------------
    Ensures we only pass safe currents that the module can physically handle
    """
    arm_input = [0] * 16
    for i in range(len(tau)):
        if(i != 0 and i != 4 and i != 8 and i != 12):
            if tau[i][0].item() < 0:
                # negative current case
                if tau[i][0].item() < -max_torque:
                    arm_input[i] = -max_torque
                elif tau[i][0].item() > -min_torque:
                    arm_input[i] = 0
                else:
                    arm_input[i] = int(tau[i][0].item())
            else:    
                # positive current case
                if tau[i][0].item() > max_torque:
                    arm_input[i] = max_torque
                elif tau[i][0].item() < min_torque:
                    arm_input[i] = 0
                else:
                    arm_input[i] = int(tau[i][0].item())   
        else:
            if tau[i][0].item() < 0:
                # negative current case
                if tau[i][0].item() < -xm_max_torque:
                    arm_input[i] = -xm_max_torque
                elif tau[i][0].item() > -xm_min_torque:
                    arm_input[i] = 0
                else:
                    arm_input[i] = int(tau[i][0].item())
            else:    
                # positive current case
                if tau[i][0].item() > xm_max_torque:
                    arm_input[i] = xm_max_torque
                elif tau[i][0].item() < xm_min_torque:
                    arm_input[i] = 0
                else:
                    arm_input[i] = int(tau[i][0].item())   
    # print(f"Our new mod input: {mod_input}\n")
    return arm_input

def grab_cable_lens(Mod1, l1, l1_0, th0, r):
    """
    Reads current motor angles from a module to calculate the current
    cable lengths.
    """
    dtheta1 = [0, 0, 0]
    for i in range(len(Mod1)):
        th1 = to_radians(Mod1[i].get_present_pos())
        dtheta1[i] = th1 - th0[i]
    for i in range(len(l1)):
        l1[i] = l1_0[i] - (dtheta1[i] * r)
    return l1

def  grab_arm_cable_lens(pos, l, l0, th0, r):
    """
    Reads current motor angles from a module to calculate the current
    cable lengths.
    Arm = Mod class instance
    l = [l1, l2, l3]
    l1  = [l11, l21, l31]
    l0 = [l1_0, l2_0, l3_0]
    th0 = [th01, th02, th03]
    th01 = [m10, m20, m30]
    """
    pos = [pos[:3], pos[3:6], pos[6:9], pos[9:]]
    for i in range(len(pos)):
        # grab cable lengths for each module
        dth = [0, 0, 0, 0]
        for j in range(len(pos[0])):
            # into 3 motors
            th = pos[i][j]
            dth[j] = th - th0[i][j]
        for k in range(len(l[i])):
            l[i][k] = l0[i][k] - (dth[k] * r)
    return l

def  grab_helix_cable_lens(pos, l, l0, th0, r):
    """
    Reads current motor angles from a module to calculate the current
    cable lengths.
    Arm = Mod class instance
    """
    pos = [pos[:3], pos[3:6], pos[6:9]]
    print(f"pos: {pos}")
    for i in range(len(pos)):
        print(f"i: {i}")
        # grab cable lengths for each module
        dth = [0, 0, 0]
        for j in range(len(pos[0])):
            # into 3 motors
            th = pos[i][j]
            dth[j] = th - th0[i][j]
        for k in range(len(l[i])):
            l[i][k] = l0[i][k] - (dth[k] * r)
    return l

def grab_helicoid_q(l1, l2, l3, mj0, s, d):
    """
    Params: 
    @l1: cable lengths for mod 1
    @l2: cable lengths for mod 2
    @l3: cable lengths for mod 3
    @s:  arc length 
    @d:  radius of module

    Returns:
    @q:  generalized coordinates [dx, dy, dL]
    """
    phi1 = math.atan2(((math.sqrt(3)/3) * (l1[1] + l1[0] - 2 * l1[2])),(l1[0] - l1[1]))
    k1 = 2 * math.sqrt(l1[2]**2 + l1[0]**2 + l1[1]**2 - (l1[2]*l1[0]) - (l1[0] * l1[1]) - (l1[2]*l1[1]))/(d* (l1[2] + l1[0] + l1[1] + 3 * Lm))
    phi2 = math.atan2(((math.sqrt(3)/3) * (l2[1] + l2[0] - 2 * l2[2])),(l2[0] - l2[1]))
    k2 = 2 * math.sqrt(l2[2]**2 + l2[0]**2 + l2[1]**2 - (l2[2]*l2[0]) - (l2[0] * l2[1]) - (l2[2]*l2[1]))/(d* (l2[2] + l2[0] + l2[1]+ 3 * Lm))
    phi3 = math.atan2(((math.sqrt(3)/3) * (l3[1] + l3[0] - 2 * l3[2])),(l3[0] - l3[1]))
    k3 = 2 * math.sqrt(l3[2]**2 + l3[0]**2 + l3[1]**2 - (l3[2]*l3[0]) - (l3[0] * l3[1]) - (l3[2]*l3[1]))/(d* (l3[2] + l3[0] + l3[1]+ 3 * Lm))

    s_curr1 = (l1[2] + l1[0] + l1[1])/3
    s_curr2 = (l2[2] + l2[0] + l2[1])/3
    s_curr3 = (l3[2] + l3[0] + l3[1])/3

    dL1 = s_curr1 - s
    dx1 = k1 * s_curr1 * d * cos(phi1)
    dy1 = k1 * s_curr1 * d * sin(phi1)
    dL2 = s_curr2 - s
    dx2 = k2 * s_curr2 * d * cos(phi2)
    dy2 = k2 * s_curr2 * d * sin(phi2)
    dL3 = s_curr3 - s 
    dx3 = k3 * s_curr3 * d * cos(phi3)
    dy3 = k3 * s_curr3 * d * sin(phi3)

    q = np.array([mj0, dx1, dy1, dL1, dx2, dy2, dL2, dx3, dy3, dL3]).reshape(-1,1)
    return q

def grab_qd(qd_str):
    """
    Returns the generalized coordinates of your desired state based off the following:
    for more detailed understanding of what these variables do, 
    refer to Fig 2. of Casimo paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8961972&tag=1

    ki   : curvature of module i
    phii : angle of rotation of curvature
    dLi  : change in arc length (current arc length - neutral arc length)
    ex. phi = 0 radian has curvature k entirely in direction of x-axis
        phi = 1.5708 (90 degrees) has curvature k in direction of y-axis
    """
    k1 = qd_str['k1']
    k2 = qd_str['k2']
    k3 = qd_str['k3']
    phi1 = qd_str['phi1']
    phi2 = qd_str['phi2']
    phi3 = qd_str['phi3']
    dL1 = qd_str['dL1']
    dL2 = qd_str['dL2']
    dL3 = qd_str['dL3']
    mj0 = qd_str['mj0']
    dx1 = k1 * d * cos(phi1)
    dy1 = k1 * d * sin(phi1)
    dx2 = k2 * d * cos(phi2)
    dy2 = k2 * d * sin(phi2)
    dx3 = k3 * d * cos(phi3)
    dy3 = k3 * d * sin(phi3)
    q = np.array([mj0, dx1, dy1, dL1, dx2, dy2, dL2, dx3, dy3, dL3]).reshape(-1,1)
    return q
