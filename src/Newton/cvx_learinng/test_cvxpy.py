import numpy as np
import matplotlib.pyplot as plt

# Given data
B_T = np.array([
    [2., 2., 0., 1., 1., 0.],
    [1., 0., 1., 0., 1., 0.]
])
c = np.array([0.16001367, 0.05333789])

# Solve B^T @ x = c symbolically/numerically
def feasible_region(phi_3, phi_4, phi_5):
    # Solve for phi_0, phi_1, phi_2 using the derived symbolic expressions
    phi_0 = (c[1] - phi_4)  # From the second equation
    phi_2 = phi_0 - (c[1] - phi_4)  # Substitute into the first equation
    phi_1 = (c[0] - 2*phi_0 - phi_3 - phi_4) / 2
    return phi_0, phi_1, phi_2

# Generate a grid of values for phi_3, phi_4, phi_5
phi_3_vals = np.linspace(0, 1, 10)
phi_4_vals = np.linspace(0, 1, 10)
phi_5_vals = np.linspace(0, 1, 10)

# Compute points in the feasible region
phi_0_list, phi_1_list, phi_2_list = [], [], []
for phi_3 in phi_3_vals:
    for phi_4 in phi_4_vals:
        for phi_5 in phi_5_vals:
            phi_0, phi_1, phi_2 = feasible_region(phi_3, phi_4, phi_5)
            phi_0_list.append(phi_0)
            phi_1_list.append(phi_1)
            phi_2_list.append(phi_2)

# Convert to arrays for plotting
phi_0_list = np.array(phi_0_list)
phi_1_list = np.array(phi_1_list)
phi_2_list = np.array(phi_2_list)

# Plot the feasible region
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(phi_0_list, phi_1_list, phi_2_list, s=1, alpha=0.5)
ax.set_xlabel('phi_0')
ax.set_ylabel('phi_1')
ax.set_zlabel('phi_2')
plt.title("Feasible Region for B.T @ x = c")
plt.show()