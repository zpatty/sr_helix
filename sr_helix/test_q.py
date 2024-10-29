from utils.utilities import *
from utils.os_utils import *
from utils.dyn_functions import *
from utils.Mod import *

def main():
    # l1 = self.l[0]
    # l2 = self.l[1]
    # l3 = self.l[2]
    # m1 = [l1[2], l2[2], l3[2]]
    # m2 = [l1[1], l2[1], l3[1]]
    # m3 = [l1[0], l2[0], l3[0]]
    s = 0.19
    mj = 0.0
    m1 = [0.0, 0.0, 0.0]
    m2 = [0.0, 0.0, 0.0]
    m3 = [0.0, 0.0, 0.0]

    q = grab_helix_q(m1, m2, m3, mj, s, d)
    print(q)

if __name__ == '__main__':
    main()