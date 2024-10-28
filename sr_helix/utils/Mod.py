#!/usr/bin/env python
# -*- coding: utf-8 -*-
from utils.dynamixel_sdk import *     
import math
DDR_POS_D_GAIN              = 80
ADDR_POS_I_GAIN             = 82
ADDR_POS_P_GAIN             = 84
ADDR_TORQUE_ENABLE          = 64
ADDR_GOAL_CURRENT           = 102
ADDR_GOAL_POSITION          = 116
ADDR_GOAL_VELOCITY          = 104
ADDR_PRESENT_POSITION       = 132
ADDR_PRESENT_VELOCITY       = 128
ADDR_PRESENT_CURRENT        = 126
ADDR_PROFILE_VELOCITY       = 112
ADDR_OPERATING_MODE         = 11
EXT_POSITION_CONTROL_MODE   = 4                 # The value of Extended Position Control Mode that can be set through the Operating Mode (11)
CURRENT_CONTROL_MODE        = 0                 # The value of Current Control Mode -- controls current(torque) regardless of speed and position 
VELOCITY_CONTROL_MODE       = 1                 # The value of Velocity Control Mode 
TORQUE_ENABLE               = 1                 # Value for enabling the torque
TORQUE_DISABLE              = 0                 # Value for disabling the torque
MAX_POSITION_VALUE          = 1048575           # Of MX with firmware 2.0 and X-Series the revolution on Extended Position Control Mode is 256 rev
MIN_POSITION_VALUE          = -1048575          # important to note that this is for EXTENDED position mode
DXL_MOVING_STATUS_THRESHOLD = 20                # Dynamixel will rotate between this value
MAX_PROFILE_VELOCITY        = 50                # ranges from 0-32,767 [0.229 rev/min]
LEN_PRESENT_POSITION        = 4                 # Data Byte Length
LEN_GOAL_CURRENT            = 2                 # Byte Length for current
LEN_GOAL_POSITION           = 4                 # Data Byte Length

def twos_comp(bit_s):
    
    # print(f"our first bit: {bit_s[0]}")
    if bit_s[0] == '1':
        inverse_s = ''.join(['1' if i == '0' else '0'
                     for i in bit_s])
        new = ('1' + inverse_s)[2:]
        dec = -int(new,2)
        return dec
    else:
        return int(bit_s,2)

def to_radians(steps):
    """
    Input: takes in dynamixel motor steps pos
    Output: radian value of dynamixel
    """
    rad_val = (steps * ((2 * math.pi)/4096))
    return rad_val

class Mod:
    """
    Class that defines a module of Dyamixels for group sync read and write
    """
    def __init__(self, packetHandler, portHandler, IDS):
        self.packetHandler = packetHandler
        self.portHandler = portHandler
        self.IDS            = IDS
        # Initialize GroupSyncRead instace for Present Position
        self.groupSyncRead = GroupSyncRead(portHandler, packetHandler, ADDR_PRESENT_POSITION, LEN_PRESENT_POSITION)
        self.groupSyncReadCurrent = GroupSyncRead(portHandler, packetHandler, ADDR_PRESENT_CURRENT, 2)
        self.groupSyncReadVel = GroupSyncRead(portHandler, packetHandler, ADDR_PRESENT_VELOCITY, 4)

        # Initialize GroupSyncWrite instance
        self.groupSyncWrite = GroupSyncWrite(portHandler, packetHandler, ADDR_GOAL_CURRENT, LEN_GOAL_CURRENT)
        for ID in self.IDS:
            # Add parameter storage for Dynamixel present position value
            dxl_addparam_result = self.groupSyncRead.addParam(ID)
            dxl_addparam_result = self.groupSyncReadVel.addParam(ID)

            if dxl_addparam_result != True:
                print("[ERROR] [ID:%03d] groupSyncRead addparam failed" %ID)
                quit()
            
            dxl_addparam_result = self.groupSyncReadCurrent.addParam(ID)
            if dxl_addparam_result != True:
                print("[ERROR] [ID:%03d] groupSyncRead addparam failed" %ID)
                quit()
        print(f"[STATUS] Mod initialized with IDS: {self.IDS}\n")
    def enable_torque(self):
        #Enable Dynamixel Torques
        for ID in self.IDS:
            dxl_comm_result, dxl_error = self.packetHandler.write1ByteTxRx(self.portHandler, ID, ADDR_TORQUE_ENABLE, TORQUE_ENABLE)
            if dxl_comm_result != COMM_SUCCESS:
                print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
            elif dxl_error != 0:
                print("%s" % self.packetHandler.getRxPacketError(dxl_error))
            else:
                print(f"[STATUS] Torque Enabled for motor {ID}\n")
    def disable_torque(self):
        for ID in self.IDS:
            dxl_comm_result, dxl_error = self.packetHandler.write1ByteTxRx(self.portHandler, ID, ADDR_TORQUE_ENABLE, TORQUE_DISABLE)
            if dxl_comm_result != COMM_SUCCESS:
                print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
            elif dxl_error != 0:
                print("%s" % self.packetHandler.getRxPacketError(dxl_error))
            else:
                print(f"[STATUS] Torque Disabled for motor {ID}\n")
    def get_position(self):
        '''
        Returns position in radians
        '''
        # Syncread present position
        # note about line below- though dxl_comm_result isn't used, you still need to call txRxPacket()
        # otherwise dynamixel won't read out its position 
        dxl_comm_result = 1
        while dxl_comm_result != COMM_SUCCESS:
            dxl_comm_result = self.groupSyncRead.txRxPacket()
        address = ADDR_PRESENT_POSITION
        pos = []
        for ID in self.IDS: 
            data = self.groupSyncRead.data_dict[ID]
            # print(f"data: {data}\n")
            d3 = bin(data[3])[2:]
            # print(f"byte: {d3}\n")
            if d3[0] == '0':
                p = DXL_MAKEDWORD(DXL_MAKEWORD(self.groupSyncRead.data_dict[ID][address - self.groupSyncRead.start_address + 0],
                                              self.groupSyncRead.data_dict[ID][address - self.groupSyncRead.start_address + 1]),
                                 DXL_MAKEWORD(self.groupSyncRead.data_dict[ID][address - self.groupSyncRead.start_address + 2],
                                              self.groupSyncRead.data_dict[ID][address - self.groupSyncRead.start_address + 3])) 
                pos.append(to_radians(p))
            else:
                # TODO: test for loop
                comb = d3
                for i in range(2,-1, -1):
                    if data[i] < 128:
                        d = (8 - len(bin(data[i])[2:])) * '0' + bin(data[i])[2:]
                        comb += d
                    else:
                        comb += bin(data[i])[2:]
                p = twos_comp(comb)
                pos.append(to_radians(p))
        return pos
    def get_velocity(self):
        '''
        Returns current velocity in radians
        '''
        # Syncread present position
        # note about line below- though dxl_comm_result isn't used, you still need to call txRxPacket()
        # otherwise dynamixel won't read out its position 
        t = time.time()
        dxl_comm_result = self.groupSyncReadVel.txRxPacket()
        # dxl_comm_result = self.groupSyncReadVel.fastSyncRead()

        # print(f"[DEBUG] dt: {time.time() - t}\n") 
        # address = ADDR_PRESENT_VELOCITY
        address = ADDR_PRESENT_VELOCITY
        # print(f"address: {address}")
        pos = []
        # print(f"dict: {self.groupSyncReadVel.data_dict}")
        for ID in self.IDS:
            # print(f"ID pos: {ID}\n")
            data = self.groupSyncReadVel.data_dict[ID]
            # print(f"data: {data}\n")
            d3 = bin(data[3])[2:]
            # print(f"byte: {d3}\n")
            if d3[0] == '0':
                p = DXL_MAKEDWORD(DXL_MAKEWORD(self.groupSyncReadVel.data_dict[ID][address - self.groupSyncReadVel.start_address + 0],
                                              self.groupSyncReadVel.data_dict[ID][address - self.groupSyncReadVel.start_address + 1]),
                                 DXL_MAKEWORD(self.groupSyncReadVel.data_dict[ID][address - self.groupSyncReadVel.start_address + 2],
                                              self.groupSyncReadVel.data_dict[ID][address - self.groupSyncReadVel.start_address + 3])) 
                pos.append(((p * 0.229)/60) * 2 * math.pi)
            else:
                # TODO: test for loop
                comb = d3
                for i in range(2,-1, -1):
                    if data[i] < 128:
                        d = (8 - len(bin(data[i])[2:])) * '0' + bin(data[i])[2:]
                        comb += d
                    else:
                        comb += bin(data[i])[2:]
                p = twos_comp(comb)
                pos.append(((p * 0.229)/60) * 2 * math.pi)
        return pos
    # TODO: try implementing speed sync
    def get_current(self):
        '''
        Returns position in radians
        '''
        # Syncread present position
        # note about line below- though dxl_comm_result isn't used, you still need to call txRxPacket()
        # otherwise dynamixel won't read out its position 
        dxl_comm_result = self.groupSyncReadCurrent.txRxPacket()
        address = ADDR_PRESENT_CURRENT
        current = []
        for ID in self.IDS:
            data = self.groupSyncReadCurrent.data_dict[ID]
            # print(f"data: {data}\n")
            c_string = bin(data[1]) + bin(data[0])[2:]
            current.append(twos_comp(c_string))
            # d3 = bin(data[3])[2:]
            # # print(f"byte: {d3}\n")
            # if d3[0] == '0':
            #     p = DXL_MAKEDWORD(DXL_MAKEWORD(self.groupSyncRead.data_dict[ID][address - self.groupSyncRead.start_address + 0],
            #                                   self.groupSyncRead.data_dict[ID][address - self.groupSyncRead.start_address + 1])) 
            #     current.append(p)
            # else:
            #     # TODO: test for loop
            #     comb = d3
            #     for i in range(2,-1, -1):
            #         if data[i] < 128:
            #             d = (8 - len(bin(data[i])[2:])) * '0' + bin(data[i])[2:]
            #             comb += d
            #         else:
            #             comb += bin(data[i])[2:]
            #     p = twos_comp(comb)
            #     current.append(p)
        return current
    def set_current_cntrl_mode(self):
        # Set operating mode to current control mode (used for torque control)
        for ID in self.IDS:
            dxl_comm_result, dxl_error = self.packetHandler.write1ByteTxRx(self.portHandler, ID, ADDR_OPERATING_MODE, CURRENT_CONTROL_MODE)
            if dxl_comm_result != COMM_SUCCESS:
                print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
            elif dxl_error != 0:
                print("%s" % self.packetHandler.getRxPacketError(dxl_error))
            else:
                print(f"[STATUS] Motor {ID} operating mode changed to current control mode.")
    def set_extended_pos_mode(self):
        # Set operating mode to extended pos mode (used to manually reset arm)
        for ID in self.IDS:
            dxl_comm_result, dxl_error = self.packetHandler.write1ByteTxRx(self.portHandler, ID, ADDR_OPERATING_MODE, EXT_POSITION_CONTROL_MODE)
            if dxl_comm_result != COMM_SUCCESS:
                print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
            elif dxl_error != 0:
                print("%s" % self.packetHandler.getRxPacketError(dxl_error))
            else:
                print(f"[STATUS] Motor {ID} operating mode changed to extended pos mode.")
    def set_velocity_mode(self):
        # Set operating mode to extended pos mode (used to manually reset arm)
        for ID in self.IDS:
            dxl_comm_result, dxl_error = self.packetHandler.write1ByteTxRx(self.portHandler, ID, ADDR_OPERATING_MODE, VELOCITY_CONTROL_MODE)
            if dxl_comm_result != COMM_SUCCESS:
                print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
            elif dxl_error != 0:
                print("%s" % self.packetHandler.getRxPacketError(dxl_error))
            else:
                print(f"[STATUS] Motor {ID} operating mode changed to velocity control mode.")
    def set_max_velocity(self, velocity):
        print(f"[STATUS] Setting max velocity... \n")
        # Set operating mode to extended pos mode (used to manually reset arm)
        for ID in self.IDS:
            # Set max velocity of dynamixel
            dxl_comm_result, dxl_error = self.packetHandler.write4ByteTxRx(self.portHandler, ID, ADDR_PROFILE_VELOCITY, velocity)
            if dxl_comm_result != COMM_SUCCESS:
                print(f"[ERROR] couln't set max velocity for motor {ID}\n")
                print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
            elif dxl_error != 0:
                print(f"[ERROR] other max velocity error for motor {ID}\n")
                print("%s" % self.packetHandler.getRxPacketError(dxl_error))
            else:
                print(f"Max velocity for motor {ID} has been set")
    def send_torque_cmd(self, torques):
        '''
        Takes in list of torques for module and sends to motors at once
        '''
        for i in range(len(torques)):    
            torque = [DXL_LOBYTE(DXL_LOWORD(torques[i])), DXL_HIBYTE(DXL_LOWORD(torques[i]))]
            # Allocate goal position value into byte array

            # Add Dynamixel goal position value to the Syncwrite parameter storage
            dxl_addparam_result = self.groupSyncWrite.addParam(self.IDS[i], torque)
            if dxl_addparam_result != True:
                print("[ERROR] [ID:%03d] groupSyncWrite addparam failed" % self.IDS[i])
                quit()
        # Syncwrite goal position
        dxl_comm_result = self.groupSyncWrite.txPacket()
        if dxl_comm_result != COMM_SUCCESS:
            print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
        # else:
            # print("[STATUS] Torque command sent")
        # Clear syncwrite parameter storage
        self.groupSyncWrite.clearParam()
    def send_pos_cmd(self, pos):
        posSyncWrite = GroupSyncWrite(self.portHandler, self.packetHandler, ADDR_GOAL_POSITION, LEN_GOAL_POSITION)
        for i in range(len(pos)):
            # Allocate goal position value into byte array
            p = [DXL_LOBYTE(DXL_LOWORD(pos[i])), DXL_HIBYTE(DXL_LOWORD(pos[i])), DXL_LOBYTE(DXL_HIWORD(pos[i])), DXL_HIBYTE(DXL_HIWORD(pos[i]))]
            # Add parameter storage for Dynamixel present position value
            dxl_addparam_result = posSyncWrite.addParam(self.IDS[i], p)
            if dxl_addparam_result != True:
                print(f"[ERROR] [ID:{self.ID[i]}] groupSyncRead addparam failed")
                quit()
        # Syncwrite goal position
        dxl_comm_result = posSyncWrite.txPacket()
        if dxl_comm_result != COMM_SUCCESS:
            print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
        else:
            print("[STATUS] Pos command sent")
        # Clear syncwrite parameter storage
        self.groupSyncWrite.clearParam()

    def send_vel_cmd(self, vel):
        velSyncWrite = GroupSyncWrite(self.portHandler, self.packetHandler, ADDR_GOAL_VELOCITY, LEN_GOAL_POSITION)
        for i in range(len(vel)):
            # Allocate goal position value into byte array
            v = [DXL_LOBYTE(DXL_LOWORD(vel[i])), DXL_HIBYTE(DXL_LOWORD(vel[i])), DXL_LOBYTE(DXL_HIWORD(vel[i])), DXL_HIBYTE(DXL_HIWORD(vel[i]))]
            # Add parameter storage for Dynamixel present position value
            dxl_addparam_result = velSyncWrite.addParam(self.IDS[i], v)
            if dxl_addparam_result != True:
                print(f"[ERROR] [ID:{self.ID[i]}] groupSyncRead addparam failed")
                quit()
        # Syncwrite goal position
        dxl_comm_result = velSyncWrite.txPacket()
        if dxl_comm_result != COMM_SUCCESS:
            print("%s" % self.packetHandler.getTxRxResult(dxl_comm_result))
        else:
            print("[STATUS] Pos command sent")
        # Clear syncwrite parameter storage
        self.groupSyncWrite.clearParam()