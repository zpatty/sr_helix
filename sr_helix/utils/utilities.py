#!/usr/bin/env python
# -*- coding: utf-8 -*-
from datetime import datetime
import os
import numpy as np
import scipy
import json
import time
import copy
from utils.constants import *
from scipy.ndimage import uniform_filter1d
from ctypes import *

file_dir = os.path.dirname(os.path.realpath('__file__')) + '/sr_helix/'

def load_helix_controller():
    file_name = os.path.join(file_dir, '../../helix_control_codegen/codegen/dll/helix_controller/helix_controller.so')
    file_name = os.path.abspath(os.path.realpath(file_name))
    ctrl = CDLL(file_name)
    helix_controller = ctrl.helix_controller
    helix_controller.restype = None
    helix_controller.argtypes = np_mat_type(10), np_mat_type(10), np_mat_type(10), np_mat_type(10), np_mat_type(10),\
    c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double,\
    np_mat_type(100), np_mat_type(100), c_double, c_double, np_mat_type(3), np_mat_type(3), np_mat_type(3),\
    c_double, c_double, np_mat_type(10), np_mat_type(10), np_mat_type(3), np_mat_type(100), np_mat_type(10), \
    np_mat_type(100), np_mat_type(3)
    
    return helix_controller

def load_robot_config():
    """
    Loads general robotic module/dynamixel config contents 
    """
    file_name = os.path.join(file_dir, '../conf/robot_config.json')
    file_name = os.path.abspath(os.path.realpath(file_name))
    with open(file_name) as config:
        mod_params = json.load(config)
    return mod_params

def load_controller_config():
    """
    Loads the controller config file contents
    """
    file_name = os.path.join(file_dir, '../conf/controller_config.json')
    file_name = os.path.abspath(os.path.realpath(file_name))
    with open(file_name) as config:
        cntrl_params = json.load(config)
    return cntrl_params

def get_pos_cmds_from_cables():
    """
    Reads from cable lengths json file, where 0-1 refers to 0-100% compression
    We then scale this to motor steps and return a pos_cmd array for the set_pos_cmd defined in Mod.py
    Motor indices for each module
    mod1: 2, 5, 8
    mod2: 1, 4, 7, 
    mod3: 0, 3, 6

    msteps = [m3, m2, m1, m3, m2, m1, m3, m2, m1]
    """
    file_name = os.path.join(file_dir, '../conf/cable_lengths.json')
    file_name = os.path.abspath(os.path.realpath(file_name))

    with open(file_name) as config:
        cable_lens = json.load(config)
    mods = np.array((cable_lens["mod1"], cable_lens["mod2"], cable_lens["mod3"])).reshape(9)
    print(f"mods: {mods}\n")
    print(f"mods shape: {mods.shape}")
    dL = np.clip(mods, 0, 1)               # clip to between 0-100*
    delta_degs = dL * MAX_DELTA            # scale percentage to degrees  
    delta_steps = (delta_degs/DEG_PER_PULSE).astype(int).tolist()
    mod1 = delta_steps[:3]
    mod2 = delta_steps[3:6]
    mod3 = delta_steps[6:]
    # convert into 
    pos_cmd = [mod3[0], mod2[0], mod1[0],
               mod3[1], mod2[1], mod1[1],
               mod3[2], mod2[2], mod1[2]]
    print(f"new pos from config: {pos_cmd}\n")
    print(f"check type: {type(pos_cmd)}\n")
    print(f"check element type: {type(pos_cmd[0])}")
    return pos_cmd

def parse_setpoint():
    file_name = os.path.join(file_dir, '../conf/q.json')
    file_name = os.path.abspath(os.path.realpath(file_name))
    with open(file_name) as q_json:
        qd_str = json.load(q_json)
    return qd_str

def save_data(q_data, qd, tau_data, input_data, c_data, t_0, timestamps, config_params, qd_params, dt_loop, traj = False, x_data = None, F_fb_data = None):

    print(f"time since: {time.time() - t_0}\n")
    t = datetime.now().strftime("%m_%d_%Y_%H_%M_%S")
    folder_name =  "data/" + t
    os.makedirs(folder_name, exist_ok=True)
    if F_fb_data is None:
        scipy.io.savemat(folder_name + "/data.mat", {'q_data': q_data.T,'tau_data': tau_data.T,'time_data': timestamps,'q_desired': qd,'c_data': c_data.T,'input_data':input_data.T})
    else:
        scipy.io.savemat(folder_name + "/data.mat", {'q_data': q_data.T,'tau_data': tau_data.T,'time_data': timestamps,'q_desired': qd,'c_data': c_data.T,'input_data':input_data.T,'F_fb_data':F_fb_data.T})
    new_config = folder_name + "/config.json"
    with open(new_config, "w") as outfile:
        outfile.write(config_params)
    new_config = folder_name + "/config.json"
    with open(new_config, "w") as outfile:
        outfile.write(config_params)
                # Writing to new config.json
    if not traj:
        qd_config = folder_name + "/qd.json"
        with open(qd_config, "w") as outfile:
            outfile.write(qd_params)

def np_mat_type(dim, element_type=np.float64):
    return np.ctypeslib.ndpointer(dtype=element_type, shape=dim, ndim = 1, flags="C_CONTIGUOUS")

def helix_controller_wrapper(q,dq,qd,dqd,ddqd,xd,dxd,dxr,d,r,param,helix_controller, Lm=Lm):
    zero =  np.zeros((len(q),1))
    kb = param['kb']
    ks = param['ks']
    bb = param['bb']
    bs = param['bs']
    bm = param['bm']
    kp = param['kp']
    kl = param['kl']
    k_base = param['k_base']
    KD_mod = param['KD_mod']
    KD_m = param['KD_m']
    Kpx = param['Kpx']
    KDx = param['KDx']
    conv_pcc = param['conv_pcc']
    conv_motor = param["conv_motor"]
    Kp = np.diag([k_base,kp,kp,kl,kp,kp,kl,kp,kp,kl])
    KD = np.diag([KD_m,KD_mod,KD_mod,KD_mod,KD_mod,KD_mod,KD_mod,KD_mod,KD_mod,KD_mod])
    Kpvec = np.reshape(Kp.astype(float),(100,1))
    KDvec = np.reshape(KD.astype(float),(100,1))
    tau = np.zeros((10,1))
    tau_r = np.zeros((10,1))
    x = np.zeros((3,1))
    cont = np.zeros((4,1))
    A = np.zeros((100,1))
    M = np.zeros((100,1))
    C = np.zeros((10,1))
    cq = np.zeros((3,1))
    N = 4
    helix_controller(np.squeeze(q), np.squeeze(dq), np.squeeze(qd),  np.squeeze(zero), np.squeeze(zero), d, m, \
                        r, kb, ks, bb, bs, bm, L0, np.squeeze(Kpvec), np.squeeze(KDvec), Kpx, KDx, np.squeeze(xd), np.squeeze(dxd), np.squeeze(dxr), conv_pcc, conv_motor, \
                        np.squeeze(tau), np.squeeze(tau_r), np.squeeze(x), np.squeeze(M), np.squeeze(C), np.squeeze(A), np.squeeze(cq))
    return tau, cont

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

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start