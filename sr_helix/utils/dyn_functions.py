#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
import json
import scipy.io
import numpy as np
from datetime import datetime
from math import cos, sin
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
    arm_input = [0] * 10
    for i in range(len(tau)):
        if(i != 0):
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

def torque_to_current(tau,l):
    tau_clip = np.concatenate((np.zeros((1,1)), tau[1:]))
    tau_cables = np.maximum(tau_clip,-30 * np.ones((10,1)))
    # put back og joint torque values since we didn't want to clip those
    tau_cables[0,0] = tau[0,0]
    arm_input = grab_arm_current(tau_cables, min_torque, max_torque)
    mod_clip =  arm_input[1:]
    for mod in range(len(limits)):
        for cable in range(len(limits[0])):
            idx = (mod * 3) + cable
            idx2 = (mod * 4) + cable
            if mod_clip[idx] < 0:
                if l[mod][cable] >= limits[mod][cable]:
                    mod_clip[idx] = 0
    mod_cmds = mod_clip
    return arm_input, mod_cmds

def grab_helix_cable_lens(pos, l, l0, th0, r):
    """
    Reads current motor angles from a module to calculate the current
    cable lengths.
    Arm = Mod class instance
    """
    # print("[DEBUG] in cable method")
    pos = [pos[:3], pos[3:6], pos[6:9]]
    # print(f"[DEBUG] pos: {pos}")
    # print(f"[DEBUG] th0: {th0}")
    # print(f"[DEBUG] l0: {l0}")

    for i in range(len(pos)):
        # print(f"i: {i}")
        # grab cable lengths for each module
        dth = [0, 0, 0]
        for j in range(len(pos[0])):
            # into 3 motors
            th = pos[i][j]
            dth[j] = th - th0[i][j]
        for k in range(len(l[i])):
            l[i][k] = l0[i][k] - (dth[k] * r)
    return l

def grab_helix_q(l1, l2, l3, mj0, s, d):
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
    # need to flip direction for mod 1 and mod 3

    l1_1 = l1[0]
    l1_2 = l1[1]
    l1_3 = l1[2]
    phi1 = math.atan2(((math.sqrt(3)/3) * (l1_3 + l1_2 - 2 * l1_1)),(l1_2 - l1_3)) 
    k1 = 2 * math.sqrt(l1_1**2 + l1_2**2 + l1_3**2 - (l1_1*l1_2) - (l1_2 * l1_3) - (l1_1*l1_3))/(d* (l1_1 + l1_2 + l1_3+ 3 * Lm))
    
    l2_1 = l2[0]
    l2_2 = l2[2]
    l2_3 = l2[1]
    phi2 = math.atan2(((math.sqrt(3)/3) * (l2_3 + l2_2 - 2 * l2_1)),(l2_2 - l2_3)) + np.pi
    k2 = 2 * math.sqrt(l2_1**2 + l2_2**2 + l2_3**2 - (l2_1*l2_2) - (l2_2 * l2_3) - (l2_1*l2_3))/(d* (l2_1 + l2_2 + l2_3+ 3 * Lm))
    
    l3_1 = l3[1]
    l3_2 = l3[2]
    l3_3 = l3[0]
    phi3 = math.atan2(((math.sqrt(3)/3) * (l3_3 + l3_2 - 2 * l3_1)),(l3_2 - l3_3))
    k3 = 2 * math.sqrt(l3_1**2 + l3_2**2 + l3_3**2 - (l3_1*l3_2) - (l3_2 * l3_3) - (l3_1*l3_3))/(d* (l3_1 + l3_2 + l3_3+ 3 * Lm))

    s_curr1 = (l1_1 + l1_2 + l1_3)/3
    s_curr2 = (l2_1 + l2_2 + l2_3)/3
    s_curr3 = (l3_1 + l3_2 + l3_3)/3

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

def grab_helix_qd(qd_str):
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
