
import numpy as np
from matplotlib import rcParams
import corner

rcParams["font.size"] = 13
rcParams["text.usetex"] = True

# Load MCMC samples (each row is one sample, each column is a parameter)
samples = np.loadtxt("MCMCchains_200.dat")[:, :-1]

# Scale the first column by 10^11
samples[:, 0] /= 1e11
samples[:, 1] /= 0.1565

# Define parameter labels (modify based on the number of parameters)
num_params = samples.shape[1]
labels = [r'$M_c/10^{11}M_\odot$', r'$\epsilon$', r'$\alpha$', r'$\beta$', r'$\gamma$', r'$z_{\rm br}$', r'$m/10^{-22}{\rm eV}$']

# Create a corner plot with shaded confidence levels in blue
fig = corner.corner(
    samples,
    labels = labels,
    show_titles = True,
    quantiles = [0.16, 0.5, 0.84],
    levels = [0.68, 0.95],
    plot_contours = True,
    fill_contours = True,
    color='C0',
    plot_datapoints = False,
    smooth = None,
    smooth1d = None
)

# Save the plot as a PDF
output_filename = "mcmc_corner_plot.pdf"
fig.savefig(output_filename, format="pdf", bbox_inches="tight")
