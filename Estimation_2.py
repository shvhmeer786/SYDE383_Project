import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import pandas as pd

# ---- Fixed constants & geometry ----
g = 9.81                      # Gravity (m/s^2)
A_t = 0.32 * 0.26             # Tank area (m²)
rho = 997.0                   # Water density (kg/m³)
mu = 1.000e-3                 # Water dynamic viscosity (Pa·s)
h_offset = 0.1                # Total starting head (m)
sin_theta = 1/150             # Tube slope

# Tube diameters
D_main = 0.00794              # Main tube diameter (m)
D_out  = 0.0111125            # T-joint outlet diameter (m)

# Minor-loss coefficients (base)
kSharpEdge = 0.5
k_t_joint  = 1.0

# Friction factor correlation (piecewise)
def f_factor(Re):
    if Re < 2300:
        return 64.0 / Re
    elif Re > 4000:
        return 0.3164 * Re**-0.25
    else:
        f_lam  = 64.0 / Re
        f_turb = 0.3164 * Re**-0.25
        w      = (Re - 2300) / 1700.0
        return (1 - w) * f_lam + w * f_turb

# Build the ODE RHS
def make_dhdt(K_total, alpha, L, diameter, area_pipe, outlet_factor=1.0):
    def dhdt(t, h_arr):
        h = h_arr[0]
        f_guess = 0.02
        for _ in range(3):
            H = h + L * sin_theta + (h_offset - L * sin_theta)
            V = math.sqrt(2 * g * H / (1 + K_total + alpha * f_guess * L / diameter))
            Re = rho * diameter * V / mu
            f_new = f_factor(Re)
            if abs(f_new - f_guess) < 1e-5:
                break
            f_guess = f_new
        return [(area_pipe / A_t) * V * outlet_factor]
    return dhdt

# Compute drain time for straight or T-joint
def drain_time(L, K_total, alpha, t_joint=False):
    if not t_joint:
        diameter    = D_main
        area_pipe   = math.pi * diameter**2 / 4
        outlet_factor = 1.0
        K_tot = K_total + kSharpEdge
    else:
        diameter    = D_out
        area_pipe   = 2 * math.pi * (diameter**2) / 4
        outlet_factor = (math.pi * D_main**2 / 4) / area_pipe
        K_tot = K_total + kSharpEdge + k_t_joint

    rhs = make_dhdt(K_tot, alpha, L, diameter, area_pipe, outlet_factor)
    sol = solve_ivp(
        rhs,
        [0, 5000],
        [0.0],
        events=lambda t, y: y[0] - 0.08,
        max_step=0.25
    )
    return float(sol.t_events[0][0])

# Experimental data for straight tubes
lengths = np.array([0.20, 0.30, 0.40, 0.60])
t_exp   = np.array([3*60+19, 3*60+34, 4*60+26, 4*60+48])

# Residuals for least-squares calibration
def residuals(params):
    K_total, alpha = params
    t_model = [drain_time(L, K_total, alpha, t_joint=False) for L in lengths]
    return np.array(t_model) - t_exp

# Initial guess and bounds for [K_total, alpha]
x0     = [0.34, 1.0]
bounds = ([0, 0.5], [5, 2])

# Perform calibration
res = least_squares(residuals, x0, bounds=bounds)
K_calib, alpha_calib = res.x

# Generate predictions post-calibration
results = []
for L, t_e in zip(lengths, t_exp):
    t_straight = drain_time(L, K_calib, alpha_calib, t_joint=False)
    t_tjoint   = drain_time(L, K_calib, alpha_calib, t_joint=True)
    err_pct = 100 * (t_straight - t_e) / t_e
    results.append({
        "L (m)": L,
        "t_exp (s)": t_e,
        "t_pred_straight (s)": round(t_straight, 1),
        "Error_straight (%)": round(err_pct, 2),
        "t_pred_tjoint (s)": round(t_tjoint, 1)
    })

# Display calibration and prediction results
print(f"Calibrated parameters:\n  K_total = {K_calib:.4f}\n  alpha    = {alpha_calib:.4f}\n")
df = pd.DataFrame(results)
print(df.to_string(index=False))
