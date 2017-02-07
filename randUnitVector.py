
def randUnitVector():
    phi = 2*math.pi*random.random()
    theta = math.acos(1-2*random.random())
    u = np.array([(math.sin(theta)*math.cos(phi)),(math.sin(theta)*math.sin(phi)),(math.cos(theta))])
    return u
import math
import random
import numpy as np
from numpy import linalg as la
def randRotation():
    R = np.zeros((3,3),float)
    u1 = randUnitVector()
    R[:,0] = np.array([u1])
    u2 = randUnitVector()
    mu = u2 - ((u2.transpose())*u1)*u1

    while la.norm(mu) > 0.1:
        u2 = randUnitVector()
        mu = u2 - ((u2.transpose()) * u1) * u1

    R[:,1] = np.array([(mu / la.norm(mu))])

    R[:,2] = np.array([(np.cross(R[:,0],R[:,1]))])
    return(R)

R = randRotation()
print(R)
