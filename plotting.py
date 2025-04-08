
import numpy as np
from matplotlib import rcParams
import corner

rcParams["font.size"] = 13
rcParams["text.usetex"] = True

# Load MCMC samples (each row is one sample, each column is a parameter)
samples1 = np.loadtxt("MCMCchainsF.dat")
samples2 = np.loadtxt("MCMCchainsW.dat")

# Find the index of the sample with the maximum likelihood
best_fit_index1 = np.argmax(samples1[:, -1])
best_fit_index2 = np.argmax(samples2[:, -1])

best_fit_point1 = samples1[best_fit_index1]
best_fit_point2 = samples2[best_fit_index2]

# Print the best fit point
print("Best fit point for FDM:")
print(best_fit_point1)

# Print the best fit point
print("Best fit point for WDM:")
print(best_fit_point2)

# samples without likelihoods
samples1 = samples1[:, :-1]
samples2 = samples2[:, :-1]

# Cut off large m values
samples1 = samples1[samples1[:, 7] <= 1.0]
samples2 = samples2[samples2[:, 7] <= 1.0]

# Scale Mc by 10^11
samples1[:, 0] /= 1e11
samples2[:, 0] /= 1e11

# Scale epsilon by fB
samples1[:, 1] /= 0.1565
samples2[:, 1] /= 0.1565

# Compute the 95% lower bound on m22
param1 = samples1[:, 7]
param2 = samples2[:, 7]
lower1, upper1 = np.percentile(param1, [5.0, 95.0])
lower2, upper2 = np.percentile(param2, [5.0, 95.0])

# Print the constraints
print("Lower bound for FDM:")
print(10**lower1)

# Print the constraints
print("Lower bound for WDM:")
print(10**lower2)

# Define parameter labels (modify based on the number of parameters)
num_params = samples1.shape[1]
labels = [r'$M_c/10^{11}M_\odot$', r'$\epsilon$', r'$\alpha$', r'$\beta$', r'$\gamma$', r'$z_{\rm br}$', r'$\sigma_{\rm UV}$', r'$\log_{10} m_{j}$']

figure = corner.corner(
    samples1,
    range = [(2.9,4.7), (0.34,0.42), (0.5,1.1), (0.2,0.6), (0.4,3.0), (9.6,11), (0.02,0.24), (-0.2,1.0)],
    labels=labels,
    title_fmt=".2f",
    levels = [0.68, 0.95],
    quantiles = [],
    plot_contours = True,
    fill_contours = True,
    color = 'C0',
    plot_datapoints = False,
    hist_kwargs = {'density': True},
    smooth = False
)

corner.corner(
    samples2,
    range = [(2.9,4.7), (0.34,0.42), (0.5,1.1), (0.2,0.6), (0.4,3.0), (9.6,11), (0.02,0.24), (-0.2,1.0)],
    labels=labels,
    title_fmt=".2f",
    levels = [0.68, 0.95],
    quantiles = [],
    plot_contours = True,
    fill_contours = True,
    color = 'C1',
    plot_datapoints = False,
    hist_kwargs = {'density': True},
    smooth = False,
    fig=figure  # Overlay on the same figure
)

# Save the plot as a PDF
output_filename = "corner.pdf"
figure.savefig(output_filename, format="pdf", bbox_inches="tight")
