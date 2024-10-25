from utils.utilities import *
from utils.os_utils import *
from utils.dyn_functions import *
from utils.Mod import *

class HelicoidRobot():
    def __init__(self, robot_params):
        portHandlerMod              = PortHandler(robot_params['mod_device'])
        packetHandlerMod            = PacketHandler(robot_params['protocol_version'])
        portHandlerJoint            = PortHandler(robot_params['joint_device'])
        packetHandlerJoint          = PacketHandler(robot_params['protocol_version'])

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

        self.Arm = Mod(packetHandlerMod, portHandlerMod, robot_params['Arm_IDs'])
        self.Arm.set_max_velocity(robot_params['max_profile_velocity'])
        self.Arm.set_current_cntrl_mode()
        self.Arm.enable_torque()

        self.L0 = robot_params['L0']
        self.l0 = [[L0, L0, L0], [L0, L0, L0], [L0, L0, L0]]
        self.l = self.l0
        self.s = L0
        self.r = robot_params['r']
        pos = self.Arm.get_position()
        self.th0 = [pos[:3], pos[3:6], pos[6:9]]
        self.Joint = Mod(packetHandlerJoint, portHandlerJoint, robot_params['Joint_IDs'])
        self.Joint.set_current_cntrl_mode()
        self.Joint.enable_torque()

        print(f"[STATUS] Arm Position: {self.Arm.get_position()}\n")
        print(f"[STATUS] Joint Position: {self.Joint.get_position()}\n")

    def send_home(self):
        """
        Set robot back to home position 
        """
        home = 9 * [0]
        HelicoidRobot.Arm.set_extended_pos_mode()
        HelicoidRobot.Arm.send_pos_cmd(home)

    def set_position(self, pd):
        HelicoidRobot.Arm.set_extended_pos_mode()
        HelicoidRobot.Arm.send_pos_cmd(pd) 
            
    def shutdown(self):
        self.Arm.disable_torque()
        self.Joint.disable_torque()
        print("[STATUS] Successfully shut down all motors \n")


def main():
    shell_cmd = 'sudo ' + 'utils/latency_write.sh'
    os.system(shell_cmd)
    robot_params = load_robot_config()
    cntrl_params = load_controller_config()
    Robot = HelicoidRobot(robot_params)
    # try:
    while True: 
        print("\nH: Home position, W: Set point, C: Compress (or press SPACE to quit!)")
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
        elif key_input == chr(PKEY_ASCII_VALUE):
            """
            Position control mode -- set "cable lengths"
            """
        elif key_input == chr(WKEY_ASCII_VALUE):
            """
            Set Arm to specific configuration
            """
            helix_controller = load_helix_controller()
            zero = np.zeros((10,1))
            dqd = zero
            ddqd = zero
            qd_str = parse_setpoint()
            qd = grab_qd(qd_str)
            print(f"[DEBUG] qd: {qd}\n")
            xd = np.array([qd_str['xd'], qd_str['yd'], qd_str['zd']]).reshape(-1,1)
            dxd = np.zeros((3,1))
            dxr = np.zeros((3,1))
            first_time = True
            t_old = time.time()
            while 1:
                if kbhit():
                    # stop program 
                    print("[STATUS] Shutting off controller")
                    Robot.shutdown()
                    break
                else:
                    # run controller
                    pos = Robot.Arm.get_position()
                    mj = Robot.Joint.get_position()[0]
                    Robot.l = grab_helix_cable_lens(pos, Robot.l, Robot.l0, Robot.th0, r)  
                    q = grab_helicoid_q(Robot.l[0], Robot.l[1], Robot.l[2], mj, Robot.s, d)
                    print(f"[DEBUG] q: {q}")
                    print(f"[DEBUG] q shape: {q.shape}")
                    if first_time:
                        first_time = False
                        dq = np.zeros((10,1))
                        dx = np.zeros((3,1))
                        x = np.zeros((3,1))
                        dt = time.time() - t_old
                    else:
                        dq = diff(q, q_old, dt)
                    q_old = q
                    x_old = x

                    print(f"cable lens: {Robot.l}\n")   
                    tau, cont = helix_controller_wrapper(q,dq,qd,dqd,ddqd,xd,dxd,dxr,d,r,cntrl_params,helix_controller, Lm=Lm)                # Robot.Arm.send_torque_cmd(9 * [0])
                    print(f"[DEBUG] tau: {tau}\n")
    # except:
    #     print("[ERROR] Shutting down robot \n")
    #     # Robot.shutdown()
    #     traceback.print_exc()

if __name__ == '__main__':
    main()