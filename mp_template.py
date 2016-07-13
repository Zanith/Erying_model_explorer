from __future__ import division

from mpmath import mp, exp, matrix, norm, nprint
import mpmath

__author__ = 'Kyle Vitautas Lopin'


states = [[0, 0, 0],
          ['Na', 0, 0],
          ['Ca', 0, 0],
          ['Mg', 0, 0],
          [0, 'Na', 0],
          [0, 'Ca', 0],
          [0, 'Mg', 0],
          ['Na', 'Na', 0],
          ['Ca', 'Na', 0],
          ['Mg', 'Na', 0],
          [0, 0, 'Na'],
          ['Na', 'Ca', 0],
          ['Ca', 'Ca', 0],
          ['Mg', 'Ca', 0],
          [0, 0, 'Ca'],
          ['Na', 'Mg', 0],
          ['Ca', 'Mg', 0],
          ['Mg', 'Mg', 0],
          [0, 0, 'Mg'],
          ['Na', 0, 'Na'],
          ['Ca', 0, 'Na'],
          ['Mg', 0, 'Na'],
          ['Na', 0, 'Ca'],
          ['Ca', 0, 'Ca'],
          ['Mg', 0, 'Ca'],
          ['Na', 0, 'Mg'],
          ['Ca', 0, 'Mg'],
          ['Mg', 0, 'Mg'],
          [0, 'Na', 'Na'],
          [0, 'Ca', 'Na'],
          [0, 'Mg', 'Na'],
          [0, 'Na', 'Ca'],
          [0, 'Ca', 'Ca'],
          [0, 'Mg', 'Ca'],
          [0, 'Na', 'Mg'],
          [0, 'Ca', 'Mg'],
          [0, 'Mg', 'Mg'],
          ['Na', 'Na', 'Na'],
          ['Ca', 'Na', 'Na'],
          ['Mg', 'Na', 'Na'],
          ['Na', 'Ca', 'Na'],
          ['Ca', 'Ca', 'Na'],
          ['Mg', 'Ca', 'Na'],
          ['Na', 'Mg', 'Na'],
          ['Ca', 'Mg', 'Na'],
          ['Mg', 'Mg', 'Na'],
          ['Na', 'Na', 'Ca'],
          ['Ca', 'Na', 'Ca'],
          ['Mg', 'Na', 'Ca'],
          ['Na', 'Ca', 'Ca'],
          ['Ca', 'Ca', 'Ca'],
          ['Mg', 'Ca', 'Ca'],
          ['Na', 'Mg', 'Ca'],
          ['Ca', 'Mg', 'Ca'],
          ['Mg', 'Mg', 'Ca'],
          ['Na', 'Na', 'Mg'],
          ['Ca', 'Na', 'Mg'],
          ['Mg', 'Na', 'Mg'],
          ['Na', 'Ca', 'Mg'],
          ['Ca', 'Ca', 'Mg'],
          ['Mg', 'Ca', 'Mg'],
          ['Na', 'Mg', 'Mg'],
          ['Ca', 'Mg', 'Mg'],
          ['Mg', 'Mg', 'Mg']]


def erying_rate_func(voltage, ion_concs, energy_barriers, Qs):
    mp.pretty = False
    mp.dps = 80
    # k0 = 6.1*10**12
    k0 = 1
    V = voltage
    q = 1/25.  # unit is e- / kT

    Nai = ion_concs['Nai']
    Nae = ion_concs['Nae']
    Cai = ion_concs['Cai']
    Cae = ion_concs['Cae']
    Mgi = ion_concs['Mgi']
    Mge = ion_concs['Mge']

    Q12 = Qs[0]
    Q13 = Qs[1]
    Q23 = Qs[2]

    d1 = energy_barriers['distance'][0]
    d2 = energy_barriers['distance'][1]
    d3 = energy_barriers['distance'][2]
    d4 = energy_barriers['distance'][3]
    d5 = energy_barriers['distance'][4]
    d6 = energy_barriers['distance'][5]
    d7 = energy_barriers['distance'][6]

    GNa1 = energy_barriers['GNa'][0]
    GNa2 = energy_barriers['GNa'][1]
    GNa3 = energy_barriers['GNa'][2]
    GNa4 = energy_barriers['GNa'][3]
    GNa5 = energy_barriers['GNa'][4]
    GNa6 = energy_barriers['GNa'][5]
    GNa7 = energy_barriers['GNa'][6]

    GCa1 = energy_barriers['GCa'][0]
    GCa2 = energy_barriers['GCa'][1]
    GCa3 = energy_barriers['GCa'][2]
    GCa4 = energy_barriers['GCa'][3]
    GCa5 = energy_barriers['GCa'][4]
    GCa6 = energy_barriers['GCa'][5]
    GCa7 = energy_barriers['GCa'][6]

    GMg1 = energy_barriers['GMg'][0]
    GMg2 = energy_barriers['GMg'][1]
    GMg3 = energy_barriers['GMg'][2]
    GMg4 = energy_barriers['GMg'][3]
    GMg5 = energy_barriers['GMg'][4]
    GMg6 = energy_barriers['GMg'][5]
    GMg7 = energy_barriers['GMg'][6]

    global k_0_1_Na, k_1_0_Na, k_0_1_Ca, k_1_0_Ca, k_0_1_Mg, k_1_0_Mg, k_1_2_Na, k_2_1_Na
    global k_1_2_Ca, k_2_1_Ca, k_1_2_Mg, k_2_1_Mg, k_2_3_Na, k_3_2_Na, k_2_3_Ca, k_3_2_Ca
    global k_2_3_Mg, k_3_2_Mg, k_3_4_Na, k_4_3_Na, k_3_4_Ca, k_4_3_Ca, k_3_4_Mg, k_4_3_Mg

    k_0_1_Na = mp.mpf(Nai*k0*exp(-GNa1)*exp(1*q*d1*V))
    k_1_0_Na = mp.mpf(k0*exp(GNa2-GNa1)*exp(1*q*(d2-d1)*V))
    k_0_1_Ca = mp.mpf(Cai*k0*exp(-GCa1)*exp(2*q*d1*V))
    k_1_0_Ca = mp.mpf(k0*exp(GCa2-GCa1)*exp(2*q*(d2-d1)*V))
    k_0_1_Mg = mp.mpf(Mgi*k0*exp(-GMg1)*exp(2*q*d1*V))
    k_1_0_Mg = mp.mpf(k0*exp(GMg2-GMg1)*exp(2*q*(d2-d1)*V))
    k_1_2_Na = mp.mpf(k0*exp(GNa2-GNa3)*exp(1*q*(d2-d3)*V))
    k_2_1_Na = mp.mpf(k0*exp(GNa4-GNa3)*exp(1*q*(d4-d3)*V))
    k_1_2_Ca = mp.mpf(k0*exp(GCa2-GCa3)*exp(2*q*(d2-d3)*V))
    k_2_1_Ca = mp.mpf(k0*exp(GCa4-GCa3)*exp(2*q*(d4-d3)*V))
    k_1_2_Mg = mp.mpf(k0*exp(GMg2-GMg3)*exp(2*q*(d2-d3)*V))
    k_2_1_Mg = mp.mpf(k0*exp(GMg4-GMg3)*exp(2*q*(d4-d3)*V))
    k_2_3_Na = mp.mpf(k0*exp(GNa4-GNa5)*exp(1*q*(d4-d5)*V))
    k_3_2_Na = mp.mpf(k0*exp(GNa6-GNa5)*exp(1*q*(d6-d5)*V))
    k_2_3_Ca = mp.mpf(k0*exp(GCa4-GCa5)*exp(2*q*(d4-d5)*V))
    k_3_2_Ca = mp.mpf(k0*exp(GCa6-GCa5)*exp(2*q*(d6-d5)*V))
    k_2_3_Mg = mp.mpf(k0*exp(GMg4-GMg5)*exp(2*q*(d4-d5)*V))
    k_3_2_Mg = mp.mpf(k0*exp(GMg6-GMg5)*exp(2*q*(d6-d5)*V))
    k_3_4_Na = mp.mpf(k0*exp(GNa6-GNa7)*exp(1*q*(d6-d7)*V))
    k_4_3_Na = mp.mpf(Nae*k0*exp(-GNa7)*exp(1*q*(1-d7)*V))
    k_3_4_Ca = mp.mpf(k0*exp(GCa6-GCa7)*exp(2*q*(d6-d7)*V))
    k_4_3_Ca = mp.mpf(Cae*k0*exp(-GCa7)*exp(2*q*(1-d7)*V))
    k_3_4_Mg = mp.mpf(k0*exp(GMg6-GMg7)*exp(2*q*(d6-d7)*V))
    k_4_3_Mg = mp.mpf(Mge*k0*exp(-GMg7)*exp(2*q*(1-d7)*V))

    print 'check1'

    trans_matrix = matrix([
    [-(k_0_1_Na + k_0_1_Ca + k_0_1_Mg + k_4_3_Na + k_4_3_Ca + k_4_3_Mg), k_1_0_Na, k_1_0_Ca, k_1_0_Mg, 0, 0, 0, 0, 0, 0, k_3_4_Na, 0, 0, 0, k_3_4_Ca, 0, 0, 0, k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [k_0_1_Na, -(k_1_0_Na + k_1_2_Na + Q13**1 * k_4_3_Na + Q13**2 * k_4_3_Ca + Q13**2 * k_4_3_Mg), 0, 0, k_2_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**1 * k_3_4_Na, 0, 0, Q13**2 * k_3_4_Ca, 0, 0, Q13**2 * k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [k_0_1_Ca, 0, -(k_1_0_Ca + k_1_2_Ca + Q13**2 * k_4_3_Na + Q13**4 * k_4_3_Ca + Q13**4 * k_4_3_Mg), 0, 0, k_2_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**2 * k_3_4_Na, 0, 0, Q13**4 * k_3_4_Ca, 0, 0, Q13**4 * k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [k_0_1_Mg, 0, 0, -(k_1_0_Mg + k_1_2_Mg + Q13**2 * k_4_3_Na + Q13**4 * k_4_3_Ca + Q13**4 * k_4_3_Mg), 0, 0, k_2_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**2 * k_3_4_Na, 0, 0, Q13**4 * k_3_4_Ca, 0, 0, Q13**4 * k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, k_1_2_Na, 0, 0, -(k_2_1_Na + Q12**1 * k_0_1_Na + Q12**2 * k_0_1_Ca + Q12**2 * k_0_1_Mg + k_2_3_Na + Q23**1 * k_4_3_Na + Q23**2 * k_4_3_Ca + Q23**2 * k_4_3_Mg), 0, 0, Q12**1 * k_1_0_Na, Q12**2 * k_1_0_Ca, Q12**2 * k_1_0_Mg, k_3_2_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**1 * k_3_4_Na, 0, 0, Q23**2 * k_3_4_Ca, 0, 0, Q23**2 * k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, k_1_2_Ca, 0, 0, -(k_2_1_Ca + Q12**2 * k_0_1_Na + Q12**4 * k_0_1_Ca + Q12**4 * k_0_1_Mg + k_2_3_Ca + Q23**2 * k_4_3_Na + Q23**4 * k_4_3_Ca + Q23**4 * k_4_3_Mg), 0, 0, 0, 0, 0, Q12**2 * k_1_0_Na, Q12**4 * k_1_0_Ca, Q12**4 * k_1_0_Mg, k_3_2_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_3_4_Na, 0, 0, Q23**4 * k_3_4_Ca, 0, 0, Q23**4 * k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, k_1_2_Mg, 0, 0, -(k_2_1_Mg + Q12**2 * k_0_1_Na + Q12**4 * k_0_1_Ca + Q12**4 * k_0_1_Mg + k_2_3_Mg + Q23**2 * k_4_3_Na + Q23**4 * k_4_3_Ca + Q23**4 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * k_1_0_Na, Q12**4 * k_1_0_Ca, Q12**4 * k_1_0_Mg, k_3_2_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_3_4_Na, 0, 0, Q23**4 * k_3_4_Ca, 0, 0, Q23**4 * k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, Q12**1 * k_0_1_Na, 0, 0, -(Q12**1 * k_1_0_Na + Q13**1 * k_2_3_Na + Q12**1 * Q13**1 * Q23**1 * k_4_3_Na + Q12**1 * Q13**2 * Q23**2 * k_4_3_Ca + Q12**1 * Q13**2 * Q23**2 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**1 * k_3_2_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**1 * Q23**1 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**2 * Q23**2 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**2 * Q23**2 * k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, Q12**2 * k_0_1_Ca, 0, 0, 0, -(Q12**2 * k_1_0_Ca + Q13**2 * k_2_3_Na + Q12**2 * Q13**2 * Q23**1 * k_4_3_Na + Q12**2 * Q13**4 * Q23**2 * k_4_3_Ca + Q12**2 * Q13**4 * Q23**2 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**2 * k_3_2_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**1 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_3_4_Mg, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, Q12**2 * k_0_1_Mg, 0, 0, 0, 0, -(Q12**2 * k_1_0_Mg + Q13**2 * k_2_3_Na + Q12**2 * Q13**2 * Q23**1 * k_4_3_Na + Q12**2 * Q13**4 * Q23**2 * k_4_3_Ca + Q12**2 * Q13**4 * Q23**2 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**2 * k_3_2_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**1 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_3_4_Mg, 0, 0, 0, 0, 0, 0],
    [k_4_3_Na, 0, 0, 0, k_2_3_Na, 0, 0, 0, 0, 0, -(k_3_4_Na + k_3_2_Na + k_0_1_Na + k_0_1_Ca + k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, k_1_0_Na, k_1_0_Ca, k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, Q12**2 * k_0_1_Na, 0, 0, 0, 0, 0, -(Q12**2 * k_1_0_Na + Q13**2 * k_2_3_Ca + Q12**2 * Q13**1 * Q23**2 * k_4_3_Na + Q12**2 * Q13**2 * Q23**4 * k_4_3_Ca + Q12**2 * Q13**2 * Q23**4 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**2 * k_3_2_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**1 * Q23**2 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_3_4_Mg, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, Q12**4 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, -(Q12**4 * k_1_0_Ca + Q13**4 * k_2_3_Ca + Q12**4 * Q13**2 * Q23**2 * k_4_3_Na + Q12**4 * Q13**4 * Q23**4 * k_4_3_Ca + Q12**4 * Q13**4 * Q23**4 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**4 * k_3_2_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_3_4_Mg, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, Q12**4 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * k_1_0_Mg + Q13**4 * k_2_3_Ca + Q12**4 * Q13**2 * Q23**2 * k_4_3_Na + Q12**4 * Q13**4 * Q23**4 * k_4_3_Ca + Q12**4 * Q13**4 * Q23**4 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**4 * k_3_2_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_3_4_Mg, 0, 0, 0],
    [k_4_3_Ca, 0, 0, 0, 0, k_2_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, -(k_3_4_Ca + k_3_2_Ca + k_0_1_Na + k_0_1_Ca + k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, k_1_0_Na, k_1_0_Ca, k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, Q12**2 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * k_1_0_Na + Q13**2 * k_2_3_Mg + Q12**2 * Q13**1 * Q23**2 * k_4_3_Na + Q12**2 * Q13**2 * Q23**4 * k_4_3_Ca + Q12**2 * Q13**2 * Q23**4 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**2 * k_3_2_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**1 * Q23**2 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_3_4_Mg, 0, 0],
    [0, 0, 0, 0, 0, 0, Q12**4 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * k_1_0_Ca + Q13**4 * k_2_3_Mg + Q12**4 * Q13**2 * Q23**2 * k_4_3_Na + Q12**4 * Q13**4 * Q23**4 * k_4_3_Ca + Q12**4 * Q13**4 * Q23**4 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**4 * k_3_2_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_3_4_Mg, 0],
    [0, 0, 0, 0, 0, 0, Q12**4 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * k_1_0_Mg + Q13**4 * k_2_3_Mg + Q12**4 * Q13**2 * Q23**2 * k_4_3_Na + Q12**4 * Q13**4 * Q23**4 * k_4_3_Ca + Q12**4 * Q13**4 * Q23**4 * k_4_3_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**4 * k_3_2_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_3_4_Na, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_3_4_Ca, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_3_4_Mg],
    [k_4_3_Mg, 0, 0, 0, 0, 0, k_2_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(k_3_4_Mg + k_3_2_Mg + k_0_1_Na + k_0_1_Ca + k_0_1_Mg), 0, 0, 0, 0, 0, 0, k_1_0_Na, k_1_0_Ca, k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, Q13**1 * k_4_3_Na, 0, 0, 0, 0, 0, Q13**1 * k_2_3_Na, 0, 0, k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, -(Q13**1 * k_3_4_Na + Q13**1 * k_3_2_Na + k_1_0_Na + Q23**1 * k_1_2_Na), 0, 0, 0, 0, 0, 0, 0, 0, Q23**1 * k_2_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, Q13**2 * k_4_3_Na, 0, 0, 0, 0, 0, Q13**2 * k_2_3_Na, 0, k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q13**2 * k_3_4_Na + Q13**2 * k_3_2_Na + k_1_0_Ca + Q23**2 * k_1_2_Ca), 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_2_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, Q13**2 * k_4_3_Na, 0, 0, 0, 0, 0, Q13**2 * k_2_3_Na, k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q13**2 * k_3_4_Na + Q13**2 * k_3_2_Na + k_1_0_Mg + Q23**2 * k_1_2_Mg), 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_2_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, Q13**2 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**2 * k_2_3_Ca, 0, 0, k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, -(Q13**2 * k_3_4_Ca + Q13**2 * k_3_2_Ca + k_1_0_Na + Q23**2 * k_1_2_Na), 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_2_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, Q13**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**4 * k_2_3_Ca, 0, k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, -(Q13**4 * k_3_4_Ca + Q13**4 * k_3_2_Ca + k_1_0_Ca + Q23**4 * k_1_2_Ca), 0, 0, 0, 0, 0, 0, 0, 0, Q23**4 * k_2_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, Q13**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**4 * k_2_3_Ca, k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q13**4 * k_3_4_Ca + Q13**4 * k_3_2_Ca + k_1_0_Mg + Q23**4 * k_1_2_Mg), 0, 0, 0, 0, 0, 0, 0, 0, Q23**4 * k_2_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, Q13**2 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**2 * k_2_3_Mg, 0, 0, k_0_1_Na, 0, 0, 0, 0, 0, 0, -(Q13**2 * k_3_4_Mg + Q13**2 * k_3_2_Mg + k_1_0_Na + Q23**2 * k_1_2_Na), 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_2_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, Q13**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**4 * k_2_3_Mg, 0, k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, -(Q13**4 * k_3_4_Mg + Q13**4 * k_3_2_Mg + k_1_0_Ca + Q23**4 * k_1_2_Ca), 0, 0, 0, 0, 0, 0, 0, 0, Q23**4 * k_2_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, Q13**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q13**4 * k_2_3_Mg, k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, -(Q13**4 * k_3_4_Mg + Q13**4 * k_3_2_Mg + k_1_0_Mg + Q23**4 * k_1_2_Mg), 0, 0, 0, 0, 0, 0, 0, 0, Q23**4 * k_2_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, Q23**1 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**1 * k_1_2_Na, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**1 * k_3_4_Na + Q23**1 * k_2_1_Na + Q12**1 * Q13**1 * Q23**1 * k_0_1_Na + Q12**2 * Q13**2 * Q23**1 * k_0_1_Ca + Q12**2 * Q13**2 * Q23**1 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**1 * Q23**1 * k_1_0_Na, Q12**2 * Q13**2 * Q23**1 * k_1_0_Ca, Q12**2 * Q13**2 * Q23**1 * k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, Q23**2 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_1_2_Ca, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**2 * k_3_4_Na + Q23**2 * k_2_1_Ca + Q12**2 * Q13**1 * Q23**2 * k_0_1_Na + Q12**4 * Q13**2 * Q23**2 * k_0_1_Ca + Q12**4 * Q13**2 * Q23**2 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**1 * Q23**2 * k_1_0_Na, Q12**4 * Q13**2 * Q23**2 * k_1_0_Ca, Q12**4 * Q13**2 * Q23**2 * k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, Q23**2 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_1_2_Mg, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**2 * k_3_4_Na + Q23**2 * k_2_1_Mg + Q12**2 * Q13**1 * Q23**2 * k_0_1_Na + Q12**4 * Q13**2 * Q23**2 * k_0_1_Ca + Q12**4 * Q13**2 * Q23**2 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**1 * Q23**2 * k_1_0_Na, Q12**4 * Q13**2 * Q23**2 * k_1_0_Ca, Q12**4 * Q13**2 * Q23**2 * k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, Q23**2 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_1_2_Na, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**2 * k_3_4_Ca + Q23**2 * k_2_1_Na + Q12**1 * Q13**2 * Q23**2 * k_0_1_Na + Q12**2 * Q13**4 * Q23**2 * k_0_1_Ca + Q12**2 * Q13**4 * Q23**2 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**2 * Q23**2 * k_1_0_Na, Q12**2 * Q13**4 * Q23**2 * k_1_0_Ca, Q12**2 * Q13**4 * Q23**2 * k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, Q23**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**4 * k_1_2_Ca, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**4 * k_3_4_Ca + Q23**4 * k_2_1_Ca + Q12**2 * Q13**2 * Q23**4 * k_0_1_Na + Q12**4 * Q13**4 * Q23**4 * k_0_1_Ca + Q12**4 * Q13**4 * Q23**4 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_1_0_Na, Q12**4 * Q13**4 * Q23**4 * k_1_0_Ca, Q12**4 * Q13**4 * Q23**4 * k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, Q23**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**4 * k_1_2_Mg, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**4 * k_3_4_Ca + Q23**4 * k_2_1_Mg + Q12**2 * Q13**2 * Q23**4 * k_0_1_Na + Q12**4 * Q13**4 * Q23**4 * k_0_1_Ca + Q12**4 * Q13**4 * Q23**4 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_1_0_Na, Q12**4 * Q13**4 * Q23**4 * k_1_0_Ca, Q12**4 * Q13**4 * Q23**4 * k_1_0_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, Q23**2 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**2 * k_1_2_Na, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**2 * k_3_4_Mg + Q23**2 * k_2_1_Na + Q12**1 * Q13**2 * Q23**2 * k_0_1_Na + Q12**2 * Q13**4 * Q23**2 * k_0_1_Ca + Q12**2 * Q13**4 * Q23**2 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**2 * Q23**2 * k_1_0_Na, Q12**2 * Q13**4 * Q23**2 * k_1_0_Ca, Q12**2 * Q13**4 * Q23**2 * k_1_0_Mg, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, Q23**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**4 * k_1_2_Ca, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**4 * k_3_4_Mg + Q23**4 * k_2_1_Ca + Q12**2 * Q13**2 * Q23**4 * k_0_1_Na + Q12**4 * Q13**4 * Q23**4 * k_0_1_Ca + Q12**4 * Q13**4 * Q23**4 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_1_0_Na, Q12**4 * Q13**4 * Q23**4 * k_1_0_Ca, Q12**4 * Q13**4 * Q23**4 * k_1_0_Mg, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, Q23**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q23**4 * k_1_2_Mg, 0, 0, 0, 0, 0, 0, 0, 0, -(Q23**4 * k_3_4_Mg + Q23**4 * k_2_1_Mg + Q12**2 * Q13**2 * Q23**4 * k_0_1_Na + Q12**4 * Q13**4 * Q23**4 * k_0_1_Ca + Q12**4 * Q13**4 * Q23**4 * k_0_1_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_1_0_Na, Q12**4 * Q13**4 * Q23**4 * k_1_0_Ca, Q12**4 * Q13**4 * Q23**4 * k_1_0_Mg],
    [0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**1 * Q23**1 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**1 * Q23**1 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**1 * Q13**1 * Q23**1 * k_3_4_Na + Q12**1 * Q13**1 * Q23**1 * k_1_0_Na), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**1 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**1 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**2 * Q23**1 * k_3_4_Na + Q12**2 * Q13**2 * Q23**1 * k_1_0_Ca), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**1 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**1 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**2 * Q23**1 * k_3_4_Na + Q12**2 * Q13**2 * Q23**1 * k_1_0_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**1 * Q23**2 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**1 * Q23**2 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**1 * Q23**2 * k_3_4_Na + Q12**2 * Q13**1 * Q23**2 * k_1_0_Na), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**2 * Q23**2 * k_3_4_Na + Q12**4 * Q13**2 * Q23**2 * k_1_0_Ca), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**2 * Q23**2 * k_3_4_Na + Q12**4 * Q13**2 * Q23**2 * k_1_0_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**1 * Q23**2 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**1 * Q23**2 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**1 * Q23**2 * k_3_4_Na + Q12**2 * Q13**1 * Q23**2 * k_1_0_Na), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**2 * Q23**2 * k_3_4_Na + Q12**4 * Q13**2 * Q23**2 * k_1_0_Ca), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_4_3_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**2 * Q23**2 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**2 * Q23**2 * k_3_4_Na + Q12**4 * Q13**2 * Q23**2 * k_1_0_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**2 * Q23**2 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**2 * Q23**2 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**1 * Q13**2 * Q23**2 * k_3_4_Ca + Q12**1 * Q13**2 * Q23**2 * k_1_0_Na), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**4 * Q23**2 * k_3_4_Ca + Q12**2 * Q13**4 * Q23**2 * k_1_0_Ca), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**4 * Q23**2 * k_3_4_Ca + Q12**2 * Q13**4 * Q23**2 * k_1_0_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**2 * Q23**4 * k_3_4_Ca + Q12**2 * Q13**2 * Q23**4 * k_1_0_Na), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**4 * Q23**4 * k_3_4_Ca + Q12**4 * Q13**4 * Q23**4 * k_1_0_Ca), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**4 * Q23**4 * k_3_4_Ca + Q12**4 * Q13**4 * Q23**4 * k_1_0_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**2 * Q23**4 * k_3_4_Ca + Q12**2 * Q13**2 * Q23**4 * k_1_0_Na), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**4 * Q23**4 * k_3_4_Ca + Q12**4 * Q13**4 * Q23**4 * k_1_0_Ca), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_4_3_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**4 * Q23**4 * k_3_4_Ca + Q12**4 * Q13**4 * Q23**4 * k_1_0_Mg), 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**2 * Q23**2 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**1 * Q13**2 * Q23**2 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**1 * Q13**2 * Q23**2 * k_3_4_Mg + Q12**1 * Q13**2 * Q23**2 * k_1_0_Na), 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**4 * Q23**2 * k_3_4_Mg + Q12**2 * Q13**4 * Q23**2 * k_1_0_Ca), 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**4 * Q23**2 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**4 * Q23**2 * k_3_4_Mg + Q12**2 * Q13**4 * Q23**2 * k_1_0_Mg), 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**2 * Q23**4 * k_3_4_Mg + Q12**2 * Q13**2 * Q23**4 * k_1_0_Na), 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**4 * Q23**4 * k_3_4_Mg + Q12**4 * Q13**4 * Q23**4 * k_1_0_Ca), 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**4 * Q23**4 * k_3_4_Mg + Q12**4 * Q13**4 * Q23**4 * k_1_0_Mg), 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**2 * Q13**2 * Q23**4 * k_0_1_Na, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**2 * Q13**2 * Q23**4 * k_3_4_Mg + Q12**2 * Q13**2 * Q23**4 * k_1_0_Na), 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_0_1_Ca, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**4 * Q23**4 * k_3_4_Mg + Q12**4 * Q13**4 * Q23**4 * k_1_0_Ca), 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_4_3_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Q12**4 * Q13**4 * Q23**4 * k_0_1_Mg, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(Q12**4 * Q13**4 * Q23**4 * k_3_4_Mg + Q12**4 * Q13**4 * Q23**4 * k_1_0_Mg)]
    ])
    # print 'check starting steady state'
    # null_ss = steady_state_null_space(trans_matrix)
    # print 'check2'
    # steady_state = null_ss

    return trans_matrix


def steady_state_null_space(_matrix):
    null_space_ss = null(_matrix)
    return null_space_ss / norm(null_space_ss, 1)


def null(_matrix):
    # u, s, vh = np.linalg.svd(_matrix)
    # null_mask = (s <= eps)
    # null_space = np.compress(null_mask, vh, axis=0)
    # return np.transpose(null_space)
    U, S, V = mp.svd_r(_matrix)
    two_smallest(S)
    return V[len(V)-1, :]


def steady_state_lowest_eig(matrix):
    values, vectors = mp.eig(matrix)
    print 'values: '
    nprint(values, 8)
    real_values = []
    for num in values:
        real_values.append(num.real)
    two_smallest(real_values)
    _index = real_values.index(max(real_values))
    ss = vectors[:, _index]
    return ss


def two_smallest(num_list):
    smallest_num = mp.mpf(1000.0)
    second_smallest_num = mp.mpf(999.0)
    for num in num_list:
        num = mpmath.fabs(num)
        if num < smallest_num:
            second_smallest_num = smallest_num
            smallest_num = num
        if smallest_num < num < second_smallest_num:
            second_smallest_num = num
    print 'two smallest'
    nprint(smallest_num, 10)
    nprint(second_smallest_num, 10)


def split_vector(_complex_vector):
    # _real_vector = matrix()
    # _imag_vector = []
    # for num in _complex_vector:
    #     # _real_vector.append(num.real)
    #     _imag_vector.append(num.imag)
    # return _real_vector, _imag_vector
    return _complex_vector.apply(mp.re), _complex_vector.apply(mp.im)


def vector_divide(vector, scalar):
    new_vector = []
    for num in vector:
        new_vector.append(num / scalar)
    return new_vector


def get_plus_minus_elements(vector):
    plus_sums = mp.mpf(0.0)
    minus_sums = mp.mpf(0.0)
    for num in vector:
        if num > 0:
            plus_sums += num
        else:
            minus_sums -= num
    print "sums"
    print plus_sums
    print minus_sums


def norm_vector(_vector):
    return vector_divide(_vector, norm(_vector, p=1))


def print_mp_vector(_vector, size):
    print 'check'
    for num in _vector:
        nprint(num, size)


if __name__ == "__main__":
    ion_conc = {'Nai': 0.12, 'Nae': 0.12,
                'Cai': 0.001, 'Cae': 0.002,
                'Mgi': .001, 'Mge': 0}
    energy_barriers_values = {'distance': [0.1, 0.2, 0.4, 0.5, 0.6, 0.75, 0.9],
                              'GNa': [2, -12, 2, -12, 2, -12, 8],
                              'GCa': [9, -12, 9, -12, 9, -12, 6],
                              'GMg': [9, -2, 9, -2, 9, -2, 9]}
    # energy_barriers_values = {'distance': [0.1, 0.2, 0.4, 0.5, 0.6, 0.75, 0.9],
    #                           'GNa': [2, -2, 2, -2, 2, -2, 2],
    #                           'GCa': [2, -2, 2, -2, 2, -2, 2],
    #                           'GMg': [2, -2, 2, -2, 2, -2, 2]}
    Q_values = [1.0, 1.0, 1.0]
    matrix = erying_rate_func(0, ion_conc, energy_barriers_values, Q_values)
    ss_null = steady_state_null_space(matrix)
    print 'ss_null'
    print_mp_vector(ss_null, 10)


    ss_eig = steady_state_lowest_eig(matrix)
    print 'ss_eig'
    nprint(ss_eig, 10)

    ss2, ss_im = split_vector(ss_eig)

    ss2 = ss2 / mpmath.norm(ss2, 1)

    print "ss2 eig"
    # print_mp_vector(ss2)
    # nprint(ss2[0], 10)
    # print "ss_im"
    print_mp_vector(ss2, 10)
    print "diffs ss_eig"
    get_plus_minus_elements(ss2)
    print "diffs ss_null"
    get_plus_minus_elements(ss_null)
    print 'nul test'
    nprint(ss_null)
    null_test = matrix * ss_null.transpose()
    print "null_test:"
    nprint(mpmath.norm(null_test, 1), 10)
    nprint(null_test, 10)
    print len(ss2)
    print len(matrix), len(matrix[0, :])
    eig_test = matrix * ss2
    print "eig test"
    nprint(mpmath.norm(eig_test, 1), 10)

    nprint(eig_test, 10)
    print 'end'

    # normed_ss = norm_vector(ss2)
    # print "normed"
    # nprint(normed_ss, 10)


