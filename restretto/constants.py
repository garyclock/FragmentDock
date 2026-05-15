INF_ENERGY = 1e9
LIMIT_ENERGY = 1e2
EPS = 1e-4

XS_TYPE_C_H = 0
XS_TYPE_C_P = 1
XS_TYPE_N_P = 2
XS_TYPE_N_D = 3
XS_TYPE_N_DC = 4
XS_TYPE_N_A = 5
XS_TYPE_N_DA = 6
XS_TYPE_O_P = 7
XS_TYPE_O_D = 8
XS_TYPE_O_A = 9
XS_TYPE_O_AC = 10
XS_TYPE_O_DA = 11
XS_TYPE_S_P = 12
XS_TYPE_P_P = 13
XS_TYPE_F_H = 14
XS_TYPE_CL_H = 15
XS_TYPE_BR_H = 16
XS_TYPE_I_H = 17
XS_TYPE_MET_D = 18
XS_TYPE_OTHER = 19
XS_TYPE_DUMMY = 20
XS_TYPE_SIZE = 21
XS_TYPE_H = 22

XS_STRINGS = [
    "C_H",
    "C_P",
    "N_P",
    "N_D",
    "N_DC",
    "N_A",
    "N_DA",
    "O_P",
    "O_D",
    "O_A",
    "O_AC",
    "O_DA",
    "S_P",
    "P_P",
    "F_H",
    "Cl_H",
    "Br_H",
    "I_H",
    "Met_D",
    "Other",
    "Dummy",
]

XS_VDW_RADII = [
    1.9,
    1.9,
    1.8,
    1.8,
    1.8,
    1.8,
    1.8,
    1.7,
    1.7,
    1.7,
    1.7,
    1.7,
    2.0,
    2.1,
    1.5,
    1.8,
    2.0,
    2.2,
    1.2,
    2.0,
    1.5,
]


def xs_radius(xs_type):
    if xs_type == XS_TYPE_H:
        return 1.2
    if xs_type < 0 or xs_type >= XS_TYPE_SIZE:
        raise ValueError("invalid X-Score atom type: %r" % (xs_type,))
    return XS_VDW_RADII[xs_type]


def xs_name(xs_type):
    if xs_type == XS_TYPE_H:
        return "H"
    if xs_type < 0 or xs_type >= XS_TYPE_SIZE:
        raise ValueError("invalid X-Score atom type: %r" % (xs_type,))
    return XS_STRINGS[xs_type]


def xs_is_hydrophobic(xs_type):
    return xs_type in {XS_TYPE_C_H, XS_TYPE_F_H, XS_TYPE_CL_H, XS_TYPE_BR_H, XS_TYPE_I_H}


def xs_is_acceptor(xs_type):
    return xs_type in {XS_TYPE_N_A, XS_TYPE_N_DA, XS_TYPE_O_A, XS_TYPE_O_AC, XS_TYPE_O_DA}


def xs_is_donor(xs_type):
    return xs_type in {XS_TYPE_N_D, XS_TYPE_N_DC, XS_TYPE_N_DA, XS_TYPE_O_D, XS_TYPE_O_DA, XS_TYPE_MET_D}


def xs_hbond(t1, t2):
    return (xs_is_donor(t1) and xs_is_acceptor(t2)) or (xs_is_donor(t2) and xs_is_acceptor(t1))


def xs_is_heavy(xs_type):
    return xs_type not in {XS_TYPE_H, XS_TYPE_DUMMY}


def xs_type_from_atomic_properties(atomic_num, is_acceptor=False, is_donor=False, ob_type=""):
    ob_type = (ob_type or "").upper()
    if atomic_num == 1:
        return XS_TYPE_H
    if atomic_num == 6:
        return XS_TYPE_C_H
    if atomic_num == 7:
        if is_acceptor and is_donor:
            return XS_TYPE_N_DA
        if is_acceptor:
            return XS_TYPE_N_A
        if is_donor:
            return XS_TYPE_N_D
        return XS_TYPE_N_P
    if atomic_num == 8:
        if is_acceptor and is_donor:
            return XS_TYPE_O_DA
        if is_acceptor:
            return XS_TYPE_O_A
        if is_donor:
            return XS_TYPE_O_D
        return XS_TYPE_O_P
    if atomic_num == 15:
        return XS_TYPE_P_P
    if atomic_num == 16:
        return XS_TYPE_S_P
    if atomic_num == 9:
        return XS_TYPE_F_H
    if atomic_num == 17:
        return XS_TYPE_CL_H
    if atomic_num == 35:
        return XS_TYPE_BR_H
    if atomic_num == 53:
        return XS_TYPE_I_H
    if "DUMMY" in ob_type or atomic_num == 0:
        return XS_TYPE_DUMMY
    return XS_TYPE_OTHER


def atomic_num_from_xs_type(xs_type):
    if xs_type == XS_TYPE_H:
        return 1
    if xs_type == XS_TYPE_DUMMY:
        return 1
    if xs_type in {XS_TYPE_C_H, XS_TYPE_C_P}:
        return 6
    if xs_type in {XS_TYPE_N_P, XS_TYPE_N_D, XS_TYPE_N_DC, XS_TYPE_N_A, XS_TYPE_N_DA}:
        return 7
    if xs_type in {XS_TYPE_O_P, XS_TYPE_O_D, XS_TYPE_O_A, XS_TYPE_O_AC, XS_TYPE_O_DA}:
        return 8
    if xs_type == XS_TYPE_S_P:
        return 16
    if xs_type == XS_TYPE_P_P:
        return 15
    if xs_type == XS_TYPE_F_H:
        return 9
    if xs_type == XS_TYPE_CL_H:
        return 17
    if xs_type == XS_TYPE_BR_H:
        return 35
    if xs_type == XS_TYPE_I_H:
        return 53
    return 6
