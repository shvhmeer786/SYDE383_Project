# Tank Drainage Analysis Project

## Overview
This project implements a numerical model to predict tank drainage times through various tube configurations. The model accounts for fluid mechanics principles including friction losses, minor losses, and flow regime transitions.

## Description
The simulation models water draining from a rectangular tank through tubes of different lengths and configurations:
- **Straight tube configuration**: Single outlet tube
- **T-joint configuration**: Two outlet tubes connected via a T-junction

## Key Features
- **Friction Factor Modeling**: Implements piecewise correlation for laminar, transitional, and turbulent flow regimes
- **Minor Loss Coefficients**: Accounts for sharp-edge entrance losses and T-joint losses
- **Numerical Integration**: Uses scipy's `solve_ivp` for ODE solving
- **Model Calibration**: Calibrates against experimental data using root-finding algorithms
- **Predictive Analysis**: Generates predictions for various tube lengths and configurations

## Technical Details

### Physical Parameters
- Tank dimensions: 0.32 m × 0.26 m
- Water density: 997 kg/m³
- Dynamic viscosity: 1.000×10⁻³ Pa·s
- Tube slope: sin(θ) = 1/150
- Main tube diameter: 7.94 mm
- T-joint outlet diameter: 11.11 mm

### Loss Coefficients
- Sharp edge entrance: K = 0.5
- T-joint: K = 1.0

### Flow Regime Correlations
- **Laminar (Re < 2300)**: f = 64/Re
- **Turbulent (Re > 4000)**: f = 0.3164/Re^0.25
- **Transitional (2300 ≤ Re ≤ 4000)**: Linear interpolation

## Usage

### Prerequisites
```bash
pip install numpy scipy pandas
```

### Running the Analysis
```bash
python Estimation.py
```

### Output
The script provides:
1. Calibrated total loss coefficient (K_total)
2. Comparison table showing:
   - Tube lengths tested
   - Experimental drainage times
   - Predicted times for straight tube configuration
   - Prediction errors
   - Predicted times for T-joint configuration

## Methodology

### Mathematical Model
The drainage process is governed by the differential equation:
```
dh/dt = (A_pipe/A_tank) × V × outlet_factor
```

Where velocity V is determined by energy balance:
```
V = √(2gH / (1 + K_total + fL/D))
```

### Calibration Process
1. Uses experimental data from 0.30 m straight tube (t = 214 s)
2. Employs Brent's root-finding method to determine optimal K_total
3. Minimizes difference between predicted and experimental drainage time

### Validation
Model predictions are validated against experimental data for tube lengths:
- 0.20 m
- 0.30 m (calibration case)
- 0.40 m  
- 0.60 m

## Results Format
```
Calibrated K_total: [value]

   L (m)  t_exp (s)  t_pred_straight (s)  Error_straight (%)  t_pred_tjoint (s)
    0.20        199                 [xx]              [x.xx]             [xxx]
    0.30        214                 [xx]              [x.xx]             [xxx]
    0.40        266                 [xx]              [x.xx]             [xxx]
    0.60        288                 [xx]              [x.xx]             [xxx]
```

## Project Structure
```
SYDE383_Project/
├── Estimation.py    # Main analysis script
└── README.md       # Project documentation
```

## Author
SYDE 383 Project - Fluid Mechanics Analysis

## License
This project is for educational purposes as part of SYDE 383 coursework. 