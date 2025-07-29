import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import pandas as pd

# ---- Fixed constants & geometry ----
g = 9.81                      # Gravity (m/s^2)
A_t = 0.32 * 0.26             # Tank area (m²)
rho = 997.0                   # Water density (kg/m³)
mu = 1.000e-3                 # Water dynamic viscosity (Pa·s)
h_offset = 0.1                # Total starting head (m)
sin_theta = 1/150             # Tube slope

# Tube diameters: match other group
D_main = 0.00794              # Main tube inner diameter (m)
D_out = 0.0111125             # T-joint outlet diameter (m)

# Minor-loss coefficients
kSharpEdge = 0.5
k_t_joint = 1.0

# Friction factor correlation (piecewise)
def f_factor(Re):
    if Re < 2300:
        return 64.0 / Re
    elif Re > 4000:
        return 0.3164 * Re**-0.25
    else:
        f_lam = 64.0 / Re
        f_turb = 0.3164 * Re**-0.25
        w = (Re - 2300) / 1700.0
        return (1 - w) * f_lam + w * f_turb

# Create the ODE RHS for a given K_total, L, and pipe geometry
def make_dhdt(K_total, L, diameter, area_pipe, outlet_factor=1.0):
    def dhdt(t, h_arr):
        h = h_arr[0]
        # Converge on friction factor f_guess
        f_guess = 0.02
        for _ in range(3):
            H = h + L * sin_theta + (h_offset - L * sin_theta)
            V = math.sqrt(2 * g * H / (1 + K_total + f_guess * L / diameter))
            Re = rho * diameter * V / mu
            f_new = f_factor(Re)
            if abs(f_new - f_guess) < 1e-5:
                break
            f_guess = f_new
        # dh/dt positive: surface drop increases from 0 → 0.08
        return [(area_pipe / A_t) * V * outlet_factor]
    return dhdt

# Compute drain time for a given L and whether the T-joint is present
def drain_time(L, K_total, t_joint=False):
    if not t_joint:
        # Straight tube
        diameter = D_main
        area_pipe = math.pi * diameter**2 / 4
        outlet_factor = 1.0
        K_tot = K_total + kSharpEdge
    else:
        # T-joint: two outlets
        diameter = D_out
        area_pipe = 2 * math.pi * (diameter**2) / 4
        # Factor to scale V → effective dh/dt
        outlet_factor = (math.pi * D_main**2 / 4) / area_pipe
        K_tot = K_total + kSharpEdge + k_t_joint

    rhs = make_dhdt(K_tot, L, diameter, area_pipe, outlet_factor)
    sol = solve_ivp(
        rhs,
        [0, 5000],
        [0.0],
        events=lambda t, y: y[0] - 0.08,
        max_step=0.25
    )
    return float(sol.t_events[0][0])

# ---- Calibration on the 0.30 m straight-tube case ----
L_calib = 0.30
t_exp_calib = 3 * 60 + 34  # 214 s

def objective(K_total):
    return drain_time(L_calib, K_total, t_joint=False) - t_exp_calib

# Find the best K_total
K_best = brentq(objective, 0.0, 5.0)

# ---- Predictions for both configurations ----
lengths = [0.20, 0.30, 0.40, 0.60]
t_exp = [3*60 + 19, 3*60 + 34, 4*60 + 26, 4*60 + 48]

results = []
for L, t_e in zip(lengths, t_exp):
    t_straight = drain_time(L, K_best, t_joint=False)
    t_tjoint   = drain_time(L, K_best, t_joint=True)
    err_pct = 100 * (t_straight - t_e) / t_e
    results.append({
        "L (m)": L,
        "t_exp (s)": t_e,
        "t_pred_straight (s)": round(t_straight, 1),
        "Error_straight (%)": round(err_pct, 2),
        "t_pred_tjoint (s)": round(t_tjoint, 1)
    })

# Display results
df = pd.DataFrame(results)
print(f"Calibrated K_total: {K_best:.4f}\n")
print(df.to_string(index=False))
