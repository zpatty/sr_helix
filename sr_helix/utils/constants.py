#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
# KEYBOARD INPUTS
ESC_ASCII_VALUE             = 0x1b
SPACE_ASCII_VALUE           = 0x20
WKEY_ASCII_VALUE            = 0x77
HKEY_ASCII_VALUE            = 0x68
OKEY_ASCII_VALUE            = 0x6F
PKEY_ASCII_VALUE            = 0x70
SKEY_ASCII_VALUE            = 0x73
AKEY_ASCII_VALUE            = 0x61
DKEY_ASCII_VALUE            = 0x64
CKEY_ASCII_VALUE            = 0x63
BKEY_ASCII_VALUE            = 0x62      
UKEY_ASCII_VALUE            = 0x75      
NKEY_ASCII_VALUE            = 0x6E
IKEY_ASCII_VALUE            = 0x69     
QKEY_ASCII_VALUE            = 0x71 
RKEY_ASCII_VALUE            = 0x72
TKEY_ASCII_VALUE            = 0x74     
MOD1_VALUE                  = 0x31      # pressing 1 on keyboard
MOD2_VALUE                  = 0x32
MOD3_VALUE                  = 0x33

MAX_VELOCITY = 20
MAX_DELTA = 650                         # in degrees
DEG_PER_PULSE = 0.088                   # to convert degrees to motor steps

d = 0.04337                             # circumradius of hexagon plate (m)
r = 0.004                               # spool radius (m)
err = math.inf                          # start with large error
mass_module = .05                       # in grams
mj = mass_module*4
m = mj
mplate = 0.08                           # in kilograms
hp = 0.007
Lm = 0.013                              # distance between center of spool and centroid of bottom plate
L0 = 0.09                               # lengths when modules are not stressed (in m)

# min and max are in mA because Dynamixel takes in mA inputs for current control
max_torque = 250
# xm_max_torque = 1100
xm_max_torque = 20

min_torque = 5
xm_min_torque = 2
l1_0 = [L0, L0, L0]
l2_0 = [L0, L0, L0]
l3_0 = [L0, L0, L0]

limits = [l1_0, l2_0, l3_0]

