# -*- coding: utf-8 -*-
import numpy as np
import csv
import matplotlib.pyplot as plt

# constant setting
# temperature of Cpu K
Tcpu = 150 + 273
# temperature of cooler K
Tc = 25 + 273
# temperature of air K
Ta = 25 + 273
# hean trans coefficient between Cpu and Cooler
alpha0 = 50000
# hean trans coefficient between Cooler and air
alpha1 = 0.01
#dimension m
UNIT = 0.05
# interval number
INT_NUM = 8
# interval
h = UNIT / INT_NUM
# boundary setting
S0=1
S11=2
S12=3
S13=4
S14=5
IN=-1
OUT=0
# Cooler model reading
position_vet = []
with open('PPC.csv') as model:
    reader = csv.reader(model)
    for row in reader:
        position_vet.append(row)

y_vet_count=len(position_vet)
x_vet_count = len(position_vet[0])
y_count = (y_vet_count-1) * INT_NUM + 1
x_count = (x_vet_count-1) * INT_NUM + 1
T0_p = np.zeros([y_count,x_count])
T0 = np.zeros([y_count,x_count])
for y in range(1,y_vet_count):
    for x in range(1,x_vet_count):
        vet_1 = position_vet[y][x]
        vet_2 = position_vet[y-1][x]
        vet_3 = position_vet[y-1][x-1]
        vet_4 = position_vet[y][x-1]
        if vet_1 == 'i' or vet_4 == 'i' or (vet_1 == 's12' and vet_3 =='s13') or (vet_1 == 's12' and vet_3 =='i') or (vet_1 == 's14' and vet_3 =='s13') or (vet_1 == 's14' and vet_3 =='i') or (vet_1 == 's14' and vet_3 =='s12'):
            T0_p[(y - 1) * INT_NUM:y * INT_NUM, (x - 1) * INT_NUM:x * INT_NUM] = IN

        if (vet_1 == 's12' and vet_4 == 's12' and vet_2 == 's14' and vet_3 == 's13'):
            T0_p[y * INT_NUM,(x-1) * INT_NUM: x * INT_NUM] = S12
            T0_p[(y-1) * INT_NUM:y * INT_NUM, x * INT_NUM] = S14
            T0_p[(y-1) * INT_NUM:y * INT_NUM, (x-1) * INT_NUM] = S13
        if (vet_1 == 's13' and vet_4 == 's14' and vet_2 == 's12' and vet_3 == 's12'):
            y_s12 = y-1
            T0_p[(y-1) * INT_NUM:y * INT_NUM, x * INT_NUM] = S13
            T0_p[(y-1) * INT_NUM:y * INT_NUM, (x-1) * INT_NUM] = S14
            T0_p[(y - 1) * INT_NUM, (x - 1) * INT_NUM: x * INT_NUM] = S12
        if (vet_1 == 's14' and vet_2 == 's14' and vet_3 == 's13' and vet_4 == 's13'):
            T0_p[(y-1) * INT_NUM:y * INT_NUM, x * INT_NUM] = S14
            T0_p[(y-1) * INT_NUM:y * INT_NUM, (x-1) * INT_NUM] = S13
        if (vet_1=='s14' and vet_2 =='s14') or (vet_1=='s14' and vet_2 =='s12') or (vet_1=='s14' and vet_2 =='s11'):
            T0_p[(y-1) * INT_NUM:y * INT_NUM, x * INT_NUM] = S14
        if (vet_1=='s12' and vet_2 =='s14'):
            T0_p[(y-1) * INT_NUM:y * INT_NUM+1, x * INT_NUM] = S14
        if (vet_3=='s12' and vet_4 =='s13'):
            T0_p[(y-1) * INT_NUM:y * INT_NUM, (x-1) * INT_NUM] = S13
        if y == 1:
            vet_p = position_vet[0][x - 1]
            vet_n = position_vet[0][x]
            if vet_p == 's11' and vet_n == 's11':
                T0_p[0,(x-1)*INT_NUM:x*INT_NUM] = S11
            if (vet_p == 's11' and vet_n == 's0'):
                T0_p[0,(x-1)*INT_NUM:x*INT_NUM] = S0
                x_c_s = x-1
            if (vet_p == 's0' and vet_n == 's0'):
                T0_p[0,(x-1)*INT_NUM:x*INT_NUM] = S0
            if (vet_p == 's0' and vet_n == 's11'):
                T0_p[0,(x-1)*INT_NUM:x*INT_NUM] = S0
                x_c_e = x
        if x == 1:
            vet_p = position_vet[y-1][0]
            vet_n = position_vet[y][0]
            if vet_p == 's11' and vet_n == 's13' or (vet_p == 's13' and vet_n == 's13') or (vet_p == 's13' and vet_n == 's12') :
                T0_p[(y - 1) * INT_NUM:y * INT_NUM, 0 ] = S13
print(T0_p)
# initialize
for y in range(0,y_count):
    for x in range(0, x_count):
        if T0_p[y,x] == OUT:
            T0[y,x] = Ta
        else:
            if T0_p[y,x] == S0:
                T0[y,x] = Tcpu
            else:
                T0[y,x] = Tc
print(T0)


def s0_cal(alpha,h,t0,t1,t2,t3):
    return 1/(1-h*alpha*.5)*(-h*alpha*.5*t0+0.25*(2.*t1+t2+t3))
def s1_cal(alpha,h,t0,t1,t2,t3):
    return 1/(1+h*alpha*.5)*(h*alpha*.5*t0+0.25*(2.*t1+t2+t3))

dt = 1
T = np.arange(0,300,dt)
tmp=[]
err = 0.05
e = 1
first = True
while e > err:
    for y in range(0,y_s12*INT_NUM+1):
        for x in range(x_c_s*INT_NUM,x_c_e*INT_NUM+1):
            if y == 0:
                if x == x_c_e*INT_NUM:
                    tt = Ta
                else:
                    tt = T0[0,x+1]
                T0[y,x] = s0_cal(alpha0,h,Tcpu,T0[1,x],tt,T0[0,x-1])
            else:
                if x == x_c_e * INT_NUM:
                    tt = Ta
                else:
                    tt = T0[y, x + 1]
                if T0_p[y,x] != S12:
                    T0[y,x] = 0.25*(T0[y-1,x]+T0[y,x-1]+tt+T0[y+1,x])
                else:
                    T0[y, x] = s1_cal(alpha1, h, Ta, T0[y+1, x], tt, T0[y, x - 1])
    for x in range(x_c_e*INT_NUM,(x_vet_count-1)*INT_NUM+1):
        for y in range(0, y_s12 * INT_NUM+1):
            if x == (x_vet_count-1)*INT_NUM:
                tt = Ta
            else:
                tt = T0[y, x + 1]
            if T0_p[y, x] == S12:
                T0[y, x] = s1_cal(alpha1, h, Ta, T0[y + 1, x], tt, T0[y, x - 1])
                continue
            if T0_p[y, x] == S11:
                T0[y, x] = s0_cal(alpha1, h, Ta, T0[y + 1, x], tt, T0[y, x - 1])
                continue
            if T0_p[y, x] == S14:
                T0[y, x] = s1_cal(alpha1, h, Ta, tt, T0[y + 1, x], T0[y - 1, x])
                continue
            T0[y, x] = 0.25 * (T0[y - 1, x] + T0[y, x - 1] + tt + T0[y + 1, x])
    for x in range(x_c_s*INT_NUM,0,-1):
        for y in range(0, y_s12 * INT_NUM+1):
            if x == 0:
                tt = Ta
            else:
                tt = T0[y, x - 1]
            if T0_p[y, x] == S12:
                T0[y, x] = s1_cal(alpha1, h, Ta, T0[y + 1, x], T0[y, x + 1], tt)
                continue
            if T0_p[y, x] == S11:
                T0[y, x] = s0_cal(alpha1, h, Ta, T0[y + 1, x], T0[y, x + 1], tt)
                continue
            if T0_p[y, x] == S13:
                T0[y, x] = s0_cal(alpha1, h, Ta, tt, T0[y + 1, x], T0[y - 1, x])
                continue
            T0[y, x] = 0.25 * (T0[y - 1, x] + tt + T0[y, x + 1] + T0[y + 1, x])
    for y in range(y_s12*INT_NUM,(y_vet_count-1)*INT_NUM+1):
        for i in range(0,x_vet_count,2):
            for x in range(i*INT_NUM,(i+1)*INT_NUM+1):
                if y == (y_vet_count-1) * INT_NUM:
                    tty = Ta
                else:
                    tty = T0[y+1, x]
                if x == (x_vet_count-1) * INT_NUM:
                    ttx = Ta
                else:
                    ttx = T0[y, x+1]
                if T0_p[y, x] == S12:
                    T0[y, x] = s1_cal(alpha1, h, Ta, tty, ttx, T0[y, x - 1])
                    continue
                if T0_p[y, x] == S11:
                    T0[y, x] = s0_cal(alpha1, h, Ta, tty, ttx, T0[y, x - 1])
                    continue
                if T0_p[y, x] == S13:
                    T0[y, x] = s0_cal(alpha1, h, Ta, T0[y, x - 1], tty, T0[y - 1, x])
                    continue
                if T0_p[y, x] == S14:
                    T0[y, x] = s1_cal(alpha1, h, Ta, ttx, tty, T0[y - 1, x])
                    continue
                T0[y, x] = 0.25 * (T0[y - 1, x] + T0[y, x - 1] + T0[y, x + 1] + tty)
    if first == False:
        e=np.linalg.norm(T0-tmp)
        print(e)
    else:
        first = False
    tmp = T0.copy()
fig = plt.figure()
pcm = plt.pcolormesh(tmp)
plt.colorbar()
pcm.set_array(tmp.ravel())
plt.draw()
plt.show()