from utils.utilities import *
from utils.os_utils import *
from utils.dyn_functions import *
from utils.Mod import *

def main():
    s = 0.19
    mj = 0.0
    L0 = 0.190
    l0 = [[L0, L0, L0], [L0, L0, L0], [L0, L0, L0]]
    l = [[0,0,0], [0,0,0], [0,0,0]]

    init_th0 =  [[0.0, 0.0, 0.0], 
                 [0.0, 0.0030679615757712823, 0.0015339807878856412], 
                 [0.0015339807878856412, 0.0, 0.0015339807878856412]]
    later_th = [0.0, 0.0, 5.663457068873787,
                0.0, 0.0030679615757712823, 5.663457068873787,
                0.0015339807878856412, 0.0, 5.661923088085902]
    r = 0.03
    newl = grab_helix_cable_lens(later_th, l, l0, init_th0, r)  
    print(newl)
if __name__ == '__main__':
    main()  