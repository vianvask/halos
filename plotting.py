
import numpy as np
from matplotlib import rcParams
import corner

rcParams["font.size"] = 18
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

# scale Mc by 10^11
samples1[:, 1] /= 1e11
samples2[:, 1] /= 1e11

# scale epsilon by fB
samples1[:, 2] /= 0.1565
samples2[:, 2] /= 0.1565

# reorder columns
samples1 = samples1[:, [0,1,2,3,4,5,7,6,8,9,10,11]]
samples2 = samples2[:, [0,1,2,3,4,5,7,6,8,9,10,11]]

# Compute the 95% lower bound on m22
lower1, upper1 = np.percentile(samples1[:, 11], [5, 95])
lower2, upper2 = np.percentile(samples2[:, 11], [5, 95])

# Print the constraints
print("Lower bound for FDM:")
print(10**lower1)

# Print the constraints
print("Lower bound for WDM:")
print(10**lower2)

# define parameter labels
labels = [r'$\log_{10}(M_t/M_\odot)$', r'$M_c/10^{11}M_\odot$', r'$\epsilon$', r'$\alpha$', r'$\beta$', r'$\gamma$', r'$f_\kappa$', r'$z_\kappa$', r'$z_{\rm fb}$', r'$z_*$', r'$\sigma_{\rm UV}$', r'$\log_{10} m$']

figure = corner.corner(
    samples1,
    range = [(6,10), (3,5), (0.36,0.44), (0.64,1.1), (0.26,0.6), (0.05,0.6), (0.05,0.7), (9.8,11.4), (7,25), (20,40), (0.05,0.2), (-0.2,1.2)],
    labels = labels,
    levels = [0.68, 0.95],
    plot_contours = True,
    fill_contours = True,
    color = 'C2',
    plot_datapoints = False,
    plot_density=False,
    smooth = False
)

corner.corner(
    samples2,
    range = [(6,10), (3,5), (0.37,0.44), (0.7,1.1), (0.26,0.6), (0.05,0.6), (0.05,0.7), (9.8,11.4), (7,25), (20,40), (0.05,0.2), (-0.2,1.2)],
    labels = labels,
    levels = [0.68, 0.95],
    plot_contours = True,
    fill_contours = True,
    color = 'C1',
    plot_datapoints = False,
    plot_density = False,
    smooth = False,
    fig=figure  # Overlay on the same figure
)

axes = np.array(figure.axes).reshape((12, 12))

for i in range(12):  # loop over parameters
    q161, q501, q841 = np.percentile(samples1[:, i], [16, 50, 84])
    errm1 = q501 - q161
    errp1 = q841 - q501
    q162, q502, q842 = np.percentile(samples2[:, i], [16, 50, 84])
    errm2 = q502 - q162
    errp2 = q842 - q502
    axes[i, i].set_title(rf"$\begin{{array}}{{c}} {q501:.2f}^{{+{errp1:.2f}}}_{{-{errm1:.2f}}} \\ {q502:.2f}^{{+{errp2:.2f}}}_{{-{errm2:.2f}}} \end{{array}}$")

# change the title of m22
upper1, lower1 = np.percentile(samples1[:, 11], [95, 5])
upper2, lower2 = np.percentile(samples2[:, 11], [95, 5])
axes[11, 11].set_title(rf"$\begin{{array}}{{c}} > {lower1:.2f} \,(95\%) \\ > {lower2:.2f} \,(95\%) \end{{array}}$")

# change the titles of gamma
upper1, lower1 = np.percentile(samples1[:, 5], [95, 5])
upper2, lower2 = np.percentile(samples2[:, 5], [95, 5])
axes[5, 5].set_title(rf"$\begin{{array}}{{c}} < {upper1:.2f} \,(95\%) \\ < {upper2:.2f} \,(95\%) \end{{array}}$")

# change the title of sigmaUV
upper1, lower1 = np.percentile(samples1[:, 10], [95, 5])
upper2, lower2 = np.percentile(samples2[:, 10], [95, 5])
axes[10, 10].set_title(rf"$\begin{{array}}{{c}} < {upper1:.2f} \,(95\%) \\ < {upper2:.2f} \,(95\%) \end{{array}}$")

# change the titles of M_t
upper1, lower1 = np.percentile(samples1[:, 0], [95, 5])
upper2, lower2 = np.percentile(samples2[:, 0], [95, 5])
axes[0, 0].set_title(rf"$\begin{{array}}{{c}} < {upper1:.2f} \,(95\%) \\ < {upper2:.2f} \,(95\%) \end{{array}}$")

from matplotlib.patches import Patch
legend_elements = [
    Patch(color = 'C2', label = rf"Fuzzy DM"),
    Patch(color = 'C1', label = rf"Warm DM")
]

# Add legend above the upper-right part of the plot
figure.legend(
    handles = legend_elements,
    loc = "upper left",  # Position the legend in the upper left of the plot area
    bbox_to_anchor = (0.5, 0.8),  # Position the legend outside the plot, in the empty space
    fontsize = 30
)

# save the plot as a PDF
output_filename = "corner.pdf"
figure.savefig(output_filename, format="pdf", bbox_inches="tight")
