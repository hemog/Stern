import time
import logging
import math
from scipy.optimize import fsolve, fmin
from pyXSteam.XSteam import XSteam
import numpy as np

P6 = 0.5
T7_values = np.arange(448, 873, 50)

def calculate_eta_th(T7_value):
    steam_table = XSteam(XSteam.UNIT_SYSTEM_BARE) # m/kg/sec/K/MPa/W
      
    state_point_1 = dict()
    state_point_1["P"] = 0.0075 ## MPa
    state_point_1["h"] = steam_table.hL_p(state_point_1["P"])
    state_point_1["s"] = steam_table.sL_p(state_point_1["P"])
    state_point_1["T"] = steam_table.tsat_p(state_point_1["P"])

    state_point_2 = dict()
    state_point_2["P"] = P6
    state_point_2["T"] = state_point_1["T"]
    state_point_2["h"] = steam_table.hL_t(state_point_2["T"])
    state_point_2["s"] = steam_table.sL_t(state_point_2["T"])

    state_point_3 = dict()
    state_point_3["P"] = state_point_2["P"]
    state_point_3["h"] = steam_table.hL_p(state_point_3["P"])
    state_point_3["s"] = steam_table.sL_p(state_point_3["P"])

    state_point_4 = dict()
    state_point_4["P"] = 9.5
    state_point_4["s"] = state_point_3["s"]
    state_point_4["T"] = steam_table.t_ps(state_point_4["P"],state_point_4["s"])
    state_point_4["h"] = steam_table.hL_t(state_point_4["T"])

    state_point_5 = dict()
    state_point_5["P"] = state_point_4["P"]
    state_point_5["T"] = 873
    state_point_5["h"] = steam_table.h_pt(state_point_5["P"],state_point_5["T"])
    state_point_5["s"] = steam_table.s_pt(state_point_5["P"],state_point_5["T"])

    eta_hpt = 1
    state_point_6 = dict()
    state_point_6["P"] = state_point_3["P"]  ## MPa
    state_point_6["s_s"] = state_point_5["s"]
    state_point_6["h_s"] = steam_table.h_ps(state_point_6["P"], state_point_6["s_s"])
    state_point_6["h"] = state_point_5["h"] - eta_hpt * (
        state_point_5["h"] - state_point_6["h_s"]
    )
    state_point_6["x"] = steam_table.x_ph(state_point_6["P"], state_point_6["h"])
    state_point_6["T"]= steam_table.t_ps(state_point_6["P"],state_point_5["s"])

    state_point_7 = dict()
    state_point_7["P"] = state_point_6["P"]
    state_point_7["T"] = T7_value
    state_point_7["h"] = steam_table.h_pt(state_point_7["P"], state_point_7["T"])
    state_point_7["s"] = steam_table.s_pt(state_point_7["P"], state_point_7["T"])

    state_point_8 = dict()
    state_point_8["P"] = state_point_1["P"]
    eta_lpt = 1
    state_point_8["s"] = state_point_7["s"]
    state_point_8["h"] = steam_table.h_ps(state_point_8["P"], state_point_8["s"])
    state_point_8["x"] = steam_table.x_ps(state_point_8["P"], state_point_8["s"])

    

    m = ((state_point_3["h"]-state_point_2["h"]))/((state_point_6["h"])-state_point_2["h"])

    w_p1 = (state_point_2["h"]-state_point_1["h"])*(1-m)
    w_p2 = (state_point_4["h"]-state_point_3["h"])
    w_p = w_p1 + w_p2

    w_hpt = (state_point_5["h"]-state_point_6["h"])
    w_lpt = (state_point_7["h"]-state_point_8["h"])*(1-m)
    w_t = w_hpt + w_lpt

    w_net = (w_hpt + w_lpt) - (w_p1 + w_p2)

    q_H = state_point_5["h"]-state_point_4["h"] + (state_point_7["h"]-state_point_6["h"])*(1-m)
    q_L = (state_point_8["h"]-state_point_1["h"])*(1-m)
    q_net = q_H - q_L


    eta_th = 1- (q_L/q_H)
    
    
    return eta_th, m, w_net, q_net, w_p, w_t, q_L, q_H
    



print("Fixed Value of P6 is 0.5 MPa")
for T7_value in T7_values:
    eta_th, m, w_net, w_p, w_t, q_L, q_H, q_net = calculate_eta_th(T7_value)
        
    print(f"For T7 = {T7_value} K, eta_th = {eta_th}, m = {m}, w_net = {w_net}, Qnet = {q_net}, w_p = {w_p}, w_t = {w_t}, q_L = {q_L}, q_H = {q_H}")
