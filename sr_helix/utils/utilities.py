#!/usr/bin/env python
# -*- coding: utf-8 -*-
from datetime import datetime
import os
import numpy as np
import scipy
import json
import time
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

# def torque_to_current(tau,l):
#     # print(f"[DEBUG] OG TAU: {tau}\n")
#     tau_clip = np.concatenate((np.zeros((1,1)), tau[1:4], np.zeros((1,1)), tau[5:8],np.zeros((1,1)), tau[9:12], np.zeros((1,1)), tau[13:]))
#     tau_cables = np.maximum(tau_clip,-30 * np.ones((16,1)))
#     # put back og joint torque values since we didn't want to clip those
#     tau_cables[0,0] = tau[0,0]
#     tau_cables[4,0] = tau[4,0]
#     tau_cables[8,0] = tau[8,0]
#     tau_cables[12,0] = tau[12,0]
#     # print(f"[DEBUG] tau cables: {tau_cables}\n")
#     arm_input = grab_arm_current(tau_cables, min_torque, max_torque)
#     # print("[DEBUG] dt time:", (time.time() - tt))
#     # print(f"[DEBUG] arm_input before mod clip: {arm_input}\n")  
#     mod_clip =  arm_input[1:4] +  arm_input[5:8] + arm_input[9:12] + arm_input[13:]
#     # print(f"mod clip: {mod_clip}\n")
#     for mod in range(len(limits)):
#         for cable in range(len(limits[0])):
#             idx = (mod * 3) + cable
#             idx2 = (mod * 4) + cable
#             # for cases where we try to extend the cables
#             # mod_clip[idx] = 0
#             # if mod != 3:
#                 # mod_clip[idx] = 0
#             if mod_clip[idx] < 0:
#                 if l[mod][cable] >= limits[mod][cable]:
#                     mod_clip[idx] = 0
#                     # arm_input[idx2] = 0
#                 # if c[mod][0] < 0:
#                 #     mod_clip[idx] = 0
#     # print(f"[DEBUG] arm_input: {arm_input}\n")                
#     # send current command to motors
#     mod_cmds = mod_clip
#     # print(f"[DEBUG] mod cmds: {mod_cmds}\n")

#     return arm_input, mod_cmds

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

    helix_controller(np.squeeze(q), np.squeeze(dq), np.squeeze(qd),  np.squeeze(zero), np.squeeze(zero), d, m, \
                        r, kb, ks, bb, bs, bm, L0, np.squeeze(Kpvec), np.squeeze(KDvec), Kpx, KDx, np.squeeze(xd), np.squeeze(dxd), np.squeeze(dxr), conv_pcc, conv_motor, \
                        np.squeeze(tau), np.squeeze(tau_r), np.squeeze(x), np.squeeze(M), np.squeeze(C), np.squeeze(A), np.squeeze(cq))
    return tau, cont
# def puppet_controller_wrapper(q,dq,qd,dqd,ddqd,d,hp,mplate,r,s,param,puppet_controller_c,Lm=Lm):
#     zero =  np.zeros((len(q),1))
#     mm = param['mm']
#     mm = param['mm']
#     hm = param['hm']
#     rm = param['rm']
#     N  = param['N']
#     kb = param['kb']
#     ks = param['ks']
#     bb = param['bb']
#     bs = param['bs']
#     e  = param['e']
#     kp = param['kp']
#     kl = param['kl']
#     kp_bot = param['kp_bot']
#     kl_bot = param['kl_bot']
#     km = param['km']
#     k_base = param['k_base']
#     kp_bot  = param['kp_bot']
#     km_bottom = param['km_bottom']
#     KD_mod = param['KD_mod']
#     KD_m = param['KD_m']
#     kc = param['kc']
#     Ka = param['Ka']
#     c_offset = param['c_offset']
#     contact = param['contact']  
#     t_test = time.time()
#     # Kp = np.diag([k_base,kp_bot,kp_bot,kl_bot,km_bottom,kp,kp,kl,km_bottom,kp,kp,kl,km,kp,kp,kl])
#     Kp = np.diag([k_base,kp,kp,kl,km_bottom,kp,kp,kl,km_bottom,kp,kp,kl,km,kp,kp,kl])
#     KD = np.diag([KD_m,KD_mod,KD_mod,KD_mod,KD_m,KD_mod,KD_mod,KD_mod,KD_m,KD_mod,KD_mod,KD_mod,KD_m,KD_mod,KD_mod,KD_mod])
#     Kpvec = np.reshape(Kp.astype(float),(256,1))
#     KDvec = np.reshape(KD.astype(float),(256,1))
#     tau = np.zeros((16,1))
#     cont = np.zeros((4,1))
#     # print(f"KD: {KD}\n")
#     # print(f"[DEBUG] make arrays: {time.time() - t_test}\n") 
#     t_test = time.time()
#     puppet_controller_c(np.squeeze(q), np.squeeze(dq), np.squeeze(qd),  np.squeeze(zero), np.squeeze(zero), d, mass_module, \
#                         mm, hm, rm, r, kb, ks, bb, bs, s, np.squeeze(Kpvec), np.squeeze(KDvec), kc, Ka, c_offset, contact, \
#                             2.0, 225.0, np.squeeze(tau), np.squeeze(cont))
#     return tau, cont

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

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start