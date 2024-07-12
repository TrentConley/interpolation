import pandas as pd
import numpy as np
from scipy.interpolate import LinearNDInterpolator

# Uncomment the code below if you have one master CSV with all the data that you need
csv_file = "testing_coeffs_2.csv"

# Read CSV into a pandas DataFrame
df = pd.read_csv(csv_file)

# Extract columns into a NumPy array
data_array = df[["AoA_deg_Used", "AoS_deg", "delta_El_L", "delta_Rud_R", "CX_body"]].values

# Extract the input parameters and coefficients
inputs = data_array[:, :4]  # First 4 columns: AoA, AoS, Elevon_Left, Rudder_Right
coefficients = data_array[:, 4]  # Last column: CX_body

def interpolated_coeffs(AOA_input, AOS_input, Elevon_Left_input, Rudder_right_and_left_input):
    # This returns the interpolated coefficients for any AoA, AoS, Elevon left, and Rudder right/left as long as it is within bounds (does not extrapolate and will return nan)
    input_points = [AOA_input, AOS_input, Elevon_Left_input, Rudder_right_and_left_input]
    interpolator_func = LinearNDInterpolator(inputs, coefficients)
    return interpolator_func(input_points)[0]

query_points = [
    ([0, 0, 0, 0], -0.0185284),
    ([0.001, 0, 0, 0], -0.018527),
    ([30.001, 0, 0, 0], 0.0425226),
    ([2, 15, 0, 0], -0.0145791),
    ([2, 13.43, 13, 12], -0.0234039),
    ([5.67, 8.91, 14.32, 3.45], -0.0155981),
    ([10.23, 7.89, 2.34, 19.56], -0.00652823),
    ([4.56, 12.34, 6.78, 9.01], -0.0132203),
    ([11.11, 13.13, 14.14, 15.15], -0.00353603),
    ([16.16, 17.17, 18.18, 19.19], -0.00494464),
    ([1.23, 2.34, 3.45, 4.56], -0.0194466),
    ([7.89, 8.9, 9.01, 10.12], -0.00528553),
    ([13.14, 14.15, 15.16, 16.17], -0.00132123),
    ([18.19, 19.2, 0.21, 1.22], 0.00182494)
]

for (q, actual) in query_points:
    coeff = interpolated_coeffs(q[0], q[1], q[2], q[3])
    if np.isnan(coeff):
        print(f"Query: {q} | Coeff: NaN | Expected: {actual} | Error: NaN%")
    else:
        err = np.abs((coeff - actual) / actual) * 100
        print(f"Query: {q} | Coeff: {coeff:.8f} | Expected: {actual} | Error: {err:.2f}%")
