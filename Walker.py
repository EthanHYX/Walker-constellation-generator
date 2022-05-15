'''
Generate Walker Constellation
'''
# %%
import math
from astropy import units as u
from astropy import time
from astropy.time import Time
from astropy.coordinates import CartesianRepresentation
from poliastro.frames import Planes
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import numpy as np
import math


R_E = 6378
H_L = 1248  # LEO层轨道高度
EPOCH = time.Time("2022-04-12 04:00:00", scale="utc")  # UTC by default
LNUMP = 6  # LEO number of planes
LNUMS = 8  # LEO number of satellites in a plane

SLOTS = 1

create_sat = locals()
get_sat = locals()


def sat(h, ecc, inc, raan, argp, nu, ep):
    # six elements
    a = (R_E + h) * u.km
    ecc = ecc * u.one
    inc = inc * u.deg
    raan = raan * u.deg  
    argp = argp * u.deg  
    nu = nu * u.deg  
    orbit = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu, epoch=ep) 
    return orbit


def constellation(num_plane, num_sat, F, h, ecc, inc, raan_0, argp, nu_0, ep):
    allsat = [[0 for i in range(num_sat)] for i in range(num_plane)]
    for i in range(num_plane):
        for j in range(num_sat):
            raan = raan_0 + i * 360 / num_plane
            nu = nu_0 + j * 360 / num_sat + i * 360 * F / (num_sat * num_plane)
            nu = nu % 360
            if nu >= 180:
                nu -= 360
            allsat[i][j] = sat(h, ecc, inc, raan, argp, nu, ep)
    return allsat


leosats_0 = constellation(LNUMP, LNUMS, 1, H_L, 0, 50, 0, 0, 0, EPOCH)  # LEO walker constellation at time 0
time = 10 * u.min
print(leosats_0[0][0].r)  # xyz
neworbit = leosats_0[0][0].propagate(time)  # sat0_0 10-min propagation
print(neworbit.r)  # xyz
