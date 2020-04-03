import numpy as np
import matplotlib.pyplot as plt

c = 3e8
kB = 1.381e-23
sigma = 5.67e-8
e = 1.602e-19
eV = e
lyr = 9.461e15
AU = 1.496e11
km = 1e3

tungsten = {
    'name': 'tungsten',
    'symbol': 'W',
    'n': 6.306e28,
    'A': 183.84,
    'm': 3.0527e-25,
    'Z': 74,
    'R_v': {
        'H': {
            0.2: 4e-5 * 1e10 #vac/meter-ion
        }
    }
}

copper = {
    'name': 'copper',
    'symbol': 'Cu',
    'n': 8.491e28,
    'A': 63.546,
    'm': 1.0552e-25,
    'Z': 29,
    'R_v': {
        'H': {
            0.2: 5e-5 * 1e10 #vac/meter-ion
        }
    }
}

hydrogen = {
    'name': 'hydrogen',
    'symbol': 'H',
    'n': None,
    'A': 1,
    'm': 1.674e-27,
    'Z': 1,
    'D0': {
        'W': 5e-8, #Frauenfelder
        'Cu':1.74e-6 #Magnusson
    },
    'Ea': {
        'W': 0.21*eV, #Frauenfelder
        'Cu': 0.4353*eV #Magnusson
    },
    'Eb': {
        'W': 1.04*eV, #Lots of sources
        'Cu': 0.95*eV #Korzhavyi
    }
}

def D_eff(gas, material, v_c, n_H, d, max_dpa):
    material_symbol = material['symbol']
    gas_symbol = gas['symbol']
    D0 = gas['D0'][material_symbol]
    Ea = gas['Ea'][material_symbol]
    Eb = gas['Eb'][material_symbol]
    R = material['R_v'][gas_symbol][v_c]
    n = material['n']
    m = gas['m']
    v = v_c*c

    Teq = (1.4 * n_H * m * v**3 / 2. / sigma)**(1./4.)
    print(f'Teq: {Teq}')

    D_L = D0 * np.exp(-Ea / kB / Teq)
    dpa = n_H * d * R / n

    if max_dpa:
        dpa[dpa > max_dpa] = max_dpa

    print(f'D_L: {D_L}')

    D_eff = D_L / (1. + dpa * np.exp(Eb / kB / Teq))

    return D_eff, dpa

def D_eff_dpa_T(gas, material, T, dpa):
    material_symbol = material['symbol']
    gas_symbol = gas['symbol']
    D0 = gas['D0'][material_symbol]
    Ea = gas['Ea'][material_symbol]
    Eb = gas['Eb'][material_symbol]
    n = material['n']
    m = gas['m']

    D_L = D0 * np.exp(-Ea / kB / T)

    D_eff = D_L / (1. + dpa * np.exp(Eb / kB / T))

    return D_eff

def main():
    d_max = 1
    n_d = 10000
    d = np.logspace(-3, 3, n_d)*AU
    v_c = 0.2
    n_H_cm3 = 1.0
    n_H = n_H_cm3*1e6
    materials = [tungsten, copper]

    D_eff_sig = 1.59e-23

    figure1, axis1 = plt.subplots()
    #axis2 = axis1.twinx()

    for material in materials:
        D, dpa = D_eff(hydrogen, material, v_c, n_H, d, max_dpa = 0.2)
        axis1.loglog(d/AU, D)
    axis1.loglog(d/AU, np.full(n_d, D_eff_sig), '--', color='black')

    plt.title('Hydrogen D_eff for n_H: 0.3/cm3 v: 0.2c')
    plt.legend([material['symbol'] for material in materials] + ['L_eff < 0.1 um'])
    plt.xlabel('d [AU]')
    plt.ylabel('D_eff [m2/s]')
    #plt.axis([0.0, d_max, 1e-22, 1e-12])
    #axis1.set_xlim((-0.1, d_max))
    plt.savefig('D_eff_d.png')

    dpa = np.logspace(-9, -1, n_d)
    T = 50 + 273.15

    figure2, axis2 = plt.subplots()
    for material in materials:
        D = D_eff_dpa_T(hydrogen, material, T, dpa)
        axis2.loglog(dpa, D)
    axis2.loglog(dpa, np.full(n_d, D_eff_sig), '--', color='black')

    plt.title('Hydrogen D_eff for T: 50C')
    plt.legend([material['symbol'] for material in materials] + ['L_eff < 0.1 um'])
    plt.xlabel('dpa')
    plt.ylabel('D_eff [m2/s]')
    plt.savefig('D_eff_dpa.png')
    axis2.set_xticks([1e-8, 1e-6, 1e-4, 0.005, 1e-2,  1e-1])
    plt.show()

if __name__ == '__main__':
    main()
