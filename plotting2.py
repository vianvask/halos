
import numpy as np
from matplotlib import rcParams
import corner

rcParams["font.size"] = 13
rcParams["text.usetex"] = True

# Load MCMC samples (each row is one sample, each column is a parameter)
samples = np.loadtxt("MCMCchains.dat")

# Find the index of the sample with the maximum likelihood
best_fit_index = np.argmax(samples[:, -1])

# samples without likelihoods
samples = samples[:, :-1]
best_fit_point = samples[best_fit_index]

# Print the best fit point
print("Best fit point:")
print(best_fit_point)

# Cut off large m3 values
samples = samples[samples[:, 6] <= 1.2]

# Scale Mc by 10^11
samples[:, 0] /= 1e11

# Scale epsilon by fB
samples[:, 1] /= 0.1565

# Compute the 95% lower bound on m3
param = samples[:, 6]
lower, upper = np.percentile(param, [5.0, 95.0])

# Define parameter labels (modify based on the number of parameters)
num_params = samples.shape[1]
labels = [r'$M_c/10^{11}M_\odot$', r'$\epsilon$', r'$\alpha$', r'$\beta$', r'$\gamma$', r'$z_{\rm br}$', r'$\log_{10} m_{3}$']

# Create a corner plot with shaded confidence levels in blue
fig = corner.corner(
    samples,
    range = [(4,7.2), (0.2,0.25), (0.7,1.2), (0.1,0.6), (0,3), (9.6,11), (0,1.2)],
    labels = labels,
    title_fmt=".3f",
    show_titles = True,
    quantiles = [],
    levels = [0.68, 0.95],
    plot_contours = True,
    fill_contours = True,
    color = 'C1',
    plot_datapoints = False,
    smooth = False,
    smooth1d = False
)

# Get the list of axes from the figure
axes = np.array(fig.axes).reshape((7, 7))  # Reshape based on number of parameters

# Change the title of m22
axes[6, 6].set_title(rf"$\log_{{10}} m_{{3}} > {lower:.3f} \,(95\%)$")

# Save the plot as a PDF
output_filename = "WDM_corner.pdf"
fig.savefig(output_filename, format="pdf", bbox_inches="tight")
