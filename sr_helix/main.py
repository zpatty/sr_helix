from utils.utilities import *
from utils.os_utils import *
from utils.dyn_functions import *
from utils.Mod import *

class HelicoidRobot():
    """
    Class to control the 3-module Helicoid Robot 
    Notes:
    Mod 1: M3, M6, M9
    Mod 2: M2, M5, M8
    Mod 3: M1, M4, M7
    """
    def __init__(self, robot_params):
        portHandlerMod              = PortHandler(robot_params['mod_device'])
        packetHandlerMod            = PacketHandler(robot_params['protocol_version'])
        portHandlerJoint            = PortHandler(robot_params['joint_device'])
        packetHandlerJoint          = PacketHandler(robot_params['protocol_version'])

        print("--------------------- Initializing Robot.....-----------------------------\n")
         # Open module port
        if portHandlerMod.openPort():
            print("Succeeded to open the port")
        else:
            print("Failed to open mod port")
            print("Press any key to terminate...")
            quit()

        # Set port baudrate
        if portHandlerMod.setBaudRate(robot_params['baudrate']):
            print("Succeeded to change the baudrate")
        else:
            print("Failed to change mod baudrate")
            print("Press any key to terminate...")
            quit()

        if portHandlerJoint.openPort():
            print("Succeeded to open the port")
        else:
            print("Failed to open joint port")
            print("Press any key to terminate...")
            quit()

        if portHandlerJoint.setBaudRate(robot_params['baudrate']):
            print("Succeeded to change the baudrate")
        else:
            print("Failed to change joint baudrate")
            print("Press any key to terminate...")
            quit()

        self.max_profile_velocity = robot_params['max_profile_velocity']
        self.Arm = Mod(packetHandlerMod, portHandlerMod, robot_params['Arm_IDs'])
        self.Arm.set_max_velocity(self.max_profile_velocity)
        self.Arm.set_current_cntrl_mode()
        self.Arm.enable_torque()

        self.L0 = robot_params['L0']
        self.l0 = [[self.L0, self.L0, self.L0], [self.L0, self.L0, self.L0], [self.L0, self.L0, self.L0]]
        self.l = [[0,0,0], [0,0,0], [0,0,0]]
        self.s = self.L0
        self.r = robot_params['r']
        self.d = robot_params['d']
        pos = self.Arm.get_position()
        self.th0 = [[0,0,0], [0,0,0], [0,0,0]]

        print("-------------- INIT CONFIGURATION---------------------\n")
        print(f"[INIT] th0: {self.th0}\n")
        print(f"[INIT] current angle readings: {pos}\n")
        self.Joint = Mod(packetHandlerJoint, portHandlerJoint, robot_params['Joint_IDs'])
        self.Joint.set_current_cntrl_mode()
        self.Joint.enable_torque()

        print(f"[INIT] Arm Position: {self.Arm.get_position()}\n")
        print(f"[INIT] Joint Position: {self.Joint.get_position()}\n")
        print(f"[INIT] Arm cable lengths: {self.l}")
        print(f"[INIT] Starting Q: {self.grab_q()}")
        self.update_cable_lens()
        print("-----------------------------------------------------\n")

    def send_home(self):
        """
        Set robot back to home position 
        """
        home = 9 * [0]
        self.Arm.disable_torque()
        self.Arm.set_max_velocity(self.max_profile_velocity)
        self.Arm.set_extended_pos_mode()
        self.Arm.enable_torque()
        self.Arm.send_pos_cmd(home)
        self.update_cable_lens()

    def update_cable_lens(self):
        """
        updates cable lengths of dynamixels
        """
        pos = self.Arm.get_position()
        self.l = grab_helix_cable_lens(pos, self.l, self.l0, self.th0, self.r)  

    def grab_cable_lens(self):
        pos = self.Arm.get_position()
        self.l = grab_helix_cable_lens(pos, self.l, self.l0, self.th0, self.r)  
        return self.l
    
    def grab_q(self):
        pos = self.Arm.get_position()
        mj = self.Joint.get_position()[0]
        self.l = grab_helix_cable_lens(pos, self.l, self.l0, self.th0, self.r)  
        # print(f"in grab q here is l: {self.l}")
        l1 = self.l[0]      # m1, m2, m3
        l2 = self.l[1]      # m4, m5, m6
        l3 = self.l[2]      # m7, m8, m9
        # OG
        m1 = [l1[2], l2[2], l3[2]]  # mod 1 cable lengths [3, 6, 9]
        m2 = [l1[1], l2[1], l3[1]]  # mod 2 cable lengths [2, 5, 8]
        m3 = [l1[0], l2[0], l3[0]]  # mod 3 cable lengths [1, 4, 7]

        # print(f"m1: {m1}")  
        # print(f"m2: {m2}")
        # print(f"m3: {m3}")
        q = grab_helix_q(m1, m2, m3, mj, self.s, self.d)

        return q

    def send_torque(self, arm_tau, mod_tau):
        # joint_cmds = [arm_tau[0]]
        joint_cmds = [0]
        m1 = mod_tau[:3]
        m2 = mod_tau[3:6]
        m3 = mod_tau[6:]
        # print(f"m1: {m1}")
        # print(f"m2: {m2}")
        # print(f"m3: {m3}")
        mod_curr = [m3[0], m2[0], m1[0], m3[1], m2[1], m1[1], m3[2], m2[2], m1[2]]
        # print(f"mod curr: {mod_curr}")

        self.Arm.send_torque_cmd(mod_tau)
        self.Joint.send_torque_cmd(joint_cmds)

    def set_position(self):
        pd = get_pos_cmds_from_cables()
        self.Arm.disable_torque()
        self.Arm.set_max_velocity(self.max_profile_velocity)
        self.Arm.set_extended_pos_mode()
        self.Arm.enable_torque()

        self.Joint.disable_torque()
        self.Joint.set_max_velocity(self.max_profile_velocity)
        self.Joint.set_extended_pos_mode()
        self.Joint.enable_torque()

        self.Arm.send_pos_cmd(pd) 
            
    def shutdown(self):
        self.Arm.disable_torque()
        self.Joint.disable_torque()
        print("[STATUS] Successfully shut down all motors \n")


def main():
    shell_cmd = 'sudo ' + 'utils/latency_write.sh'
    os.system(shell_cmd)
    robot_params = load_robot_config()
    cntrl_params = load_controller_config()
    print(f"control params: {cntrl_params}")
    Robot = HelicoidRobot(robot_params)
    # try:
    while True: 
        print("\nH: Home position, W: Set point, C: Current Config (or press SPACE to quit!)")
        key_input = getch()
        if key_input == chr(SPACE_ASCII_VALUE):
            print("[STATUS] Quitting program... \n")
            Robot.shutdown()
            break
        elif key_input == chr(HKEY_ASCII_VALUE):
            """
            Set Arm to position control mode and move to neutral configuration (dL=0)
            """
            Robot.send_home()
            time.sleep(3)
            print(f"[STATUS] Robot back to home position")
            l = Robot.grab_cable_lens()
            l1 = l[0]
            l2 = l[1]
            l3 = l[2]
            print(f"[STATUS] MOD 1 cable len: {[l1[2], l2[2], l3[2]]}")
            print(f"[STATUS] MOD 2 cable len: {[l1[1], l2[1], l3[1]]}")
            print(f"[STATUS] MOD 3 cable len: {[l1[0], l2[0], l3[0]]}")
            print(f"[STATUS] th: {Robot.Arm.get_position()}\n")
        elif key_input == chr(CKEY_ASCII_VALUE):
            """
            Get current robot configuration
            """
            print(f"[STATUS] current theta readings: {Robot.Arm.get_position()}")
            l = Robot.grab_cable_lens()
            l1 = l[0]
            l2 = l[1]
            l3 = l[2]
            print(f"[STATUS] MOD 1 cable len: {[l1[2], l2[2], l3[2]]}")
            print(f"[STATUS] MOD 2 cable len: {[l1[1], l2[1], l3[1]]}")
            print(f"[STATUS] MOD 3 cable len: {[l1[0], l2[0], l3[0]]}")
            q = Robot.grab_q()
            print(f"[STATUS] current q: {q}")
            dx1 = q[1]
            dy1 = q[2]
            dL1 = q[3]

            dx2 = q[4]
            dy2 = q[5]
            dL2 = q[6]

            dx3 = q[7]
            dy3 = q[8]
            dL3 = q[9]
            print(f"[STATUS] dx: {dx3}")
            print(f"[STATUS] dy: {dy3}")
            print(f"[STATUS] dL: {dL3}")
        elif key_input == chr(PKEY_ASCII_VALUE):
            """
            Position control mode --- set "cable lengths"
            """
            print("-------------- POSITION MODE---------------------\n")
            print(f"[POS] theta pre command: {Robot.Arm.get_position()}")
            print(f"[POS_MODE] pre command cable lengths: {Robot.grab_cable_lens()}")
            Robot.set_position()
            time.sleep(2)
            l = Robot.grab_cable_lens()
            print(f"L {l}")
            l1 = l[0]
            l2 = l[1]
            l3 = l[2]
            print(f"[POS_MODE] MOD 1 cable len: {[l1[2], l2[2], l3[2]]}")
            print(f"[POS_MODE] MOD 2 cable len: {[l1[1], l2[1], l3[1]]}")
            print(f"[POS_MODE] MOD 3 cable len: {[l1[0], l2[0], l3[0]]}")
            q = Robot.grab_q()
            print(f"[POS_MODE] Helix q: {q}")
            print(f"[POS MODE] current theta readings: {Robot.Arm.get_position()}")

        elif key_input == chr(WKEY_ASCII_VALUE) or (TKEY_ASCII_VALUE):
            """
            Set Arm to specific configuration
            """
            print(f"-------------- CONTROLLER MODE---------------------\n")
            Robot.Arm.disable_torque()
            Robot.Arm.set_current_cntrl_mode()
            Robot.Arm.enable_torque()
            Robot.Joint.set_current_cntrl_mode()

            helix_controller = load_helix_controller()
            qd_str = parse_setpoint()
            qd = grab_helix_qd(qd_str, Robot.d)
            # print(f"[DEBUG] qd: {qd}\n")
            xd = np.array([qd_str['xd'], qd_str['yd'], qd_str['zd']]).reshape(-1,1)
            zero = np.zeros((10,1))
            dqd = zero
            ddqd = zero
            dxd = np.zeros((3,1))
            dxr = np.zeros((3,1))
            nq = 10
            nmod = 4
            q_data = np.zeros((nq*2, 1))
            # print(f"Q DATA SIZE: {q_data.shape}")
            tau_data = np.zeros((nq,1))
            input_data = np.zeros((nq,1))
            dt_loop = np.zeros((1,1))       # hold dt data 
            timestamps = np.zeros((1,1))
            c_data = np.zeros((nmod,1))
            x_data = np.zeros((6,1))
            first_time = True
            qd_mat, dqd_mat, ddqd_mat, tvec = grab_trajectories()
            q_old = Robot.grab_q()
            t_old = time.time()
            t_0 = time.time()
            while 1:
                if kbhit():
                    # stop program 
                    print("[STATUS] Shutting off controller")
                    Robot.shutdown()
                    first_time = True
                    print(f"------------------------ outputs --------------------\n")
                    print(f"[OUTPUT] Our desired config: {qd}\n")
                    print(f"[OUTPUT] Our last recorded q: {q}\n")
                    err = q - qd
                    print(f"[OUTPUT] Our last recorded error: {err}\n")
                    print(f"[OUTPUT] dL1 error: {err[3]}\n")
                    print(f"[OUTPUT] dL2 error: {err[6]}\n")
                    print(f"[OUTPUT] dL3 error: {err[9]}\n")

                    print(f"max dt value: {np.max(dt_loop)}\n")
                    print(f"last time: {timestamps[-1]}\n")
                    config_params = json.dumps(cntrl_params, indent=14)
                    qd_params = json.dumps(qd_str, indent=14)

                    save_data(q_data, qd, tau_data, input_data, c_data, t_0, timestamps, config_params, qd_params, dt_loop, x_data=x_data)
                    break
                else:
                    # run controller
                    try:
                        q = Robot.grab_q()
                    except:
                        q = q_old
                    # print(f"[DEBUG] q: {q}")
                    if key_input == chr(TKEY_ASCII_VALUE):
                        n = np.argmax(tvec>time.time() - t_0) - 1
                        # print(f"this takes: {time.time() - tt}\n")
                        qd = np.array(qd_mat[:, n]).reshape(-1,1)
                        dqd = np.array(dqd_mat[:, n]).reshape(-1,1)
                        ddqd = np.array(ddqd_mat[:, n]).reshape(-1,1)
                    else:
                        dqd = zero
                        ddqd = zero

                    if first_time:
                        first_time = False
                        dq = np.zeros((10,1))
                        dx = np.zeros((3,1))
                        x = np.zeros((3,1))
                    else:
                        t = time.time()
                        time_elapsed = t-t_0
                        timestamps = np.append(timestamps, time_elapsed) 

                        dt = t - t_old
                        dq = diff(q, q_old, dt)
                
                    q_old = q
                    x_old = x
                    q_data=np.append(q_data, np.append(q,dq).reshape(-1,1), axis = 1) 
                    x_data=np.append(x_data, np.append(x,dx).reshape(-1,1), axis = 1) 
                    err = q - qd
                    err_dot = dq
                    # print(f"err: {err}\n")
                    tau, A, C = helix_controller_wrapper(q,dq,qd,dqd,ddqd,xd,dxd,dxr,Robot.d,Robot.r,cntrl_params,helix_controller, Lm=Lm)                
                    A = A.reshape((10,10),order = 'F')  

                    # tau = np.zeros((10,1))
                    # tau[4:] = np.zeros((6,1))
                    # tau = np.linalg.inv(A) @ tau
                    # tau[1] = 6
                    force_tau = A.dot(tau)
                    # print(f"forces tau: {force_tau}")
                    # print(f"tau size: {tau.size}")
                    # print(f"C: {C}")
                    # tau= [m0, m3, m6, m9, m2, m8, m5, m4, m7, m1]
                    a_tau = np.array([tau[0], 
                                      tau[9], tau[4], tau[1],
                                      tau[7], tau[6], tau[2],
                                      tau[8], tau[5], tau[3]])
                    # print(f"[DEBUG] tau: {tau}\n")
                    # print(f"realigned motors: {a_tau}")
                    arm_input, mod_cmds = torque_to_current(a_tau,Robot.l)
                    input_data=np.append(input_data, np.array(arm_input).reshape(-1,1), axis=1) 
                    tau_data=np.append(tau_data, tau, axis=1) 

                    # print(f"[DEBUG] mod cmds: {mod_cmds}\n")
                    # print(f"[DEBUG] arm input: {arm_input}\n")
                    # print(f"[DEBUG] joint input: {arm_input[0]}\n")

                    # mod_cmds = arm_input
                    Robot.send_torque(arm_input, mod_cmds)
                    

    # except:
    #     print("[ERROR] Shutting down robot \n")
    #     # Robot.shutdown()
    #     traceback.print_exc()

if __name__ == '__main__':
    main()