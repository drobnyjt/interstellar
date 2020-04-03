import numpy as np
import matplotlib.pyplot as plt

c = 2.99792e8
amu = 1.66053907e-27
m_H = 1.00727627*amu
m_He = 4.002602*amu
q = 1.602e-19

def ke(m, v):
    gamma = 1./np.sqrt(1. - v**2/c**2)
    breakpoint()
    return m*gamma*c**2 - m*c**2

print(ke(m_He, 0.1*c)/q/1e6)
print(ke(m_He, 0.3*c)/q/1e6)
