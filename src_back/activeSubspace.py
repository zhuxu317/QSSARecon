#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import numpy as np
import matplotlib
# headless
os.environ['MPLBACKEND'] = 'Agg'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import active_subspaces as ac

UQ_name   = "UQ_OH_analysis_noise_10"
# UQ_name   = "UQ_analysis"

# input_species = ["T", "X_H2", "X_O2", "X_H2O", "X_N2", "X_NH3"]
input_species = ["T", "X_H2", "X_O2", "X_H2O", "X_N2", "X_NH3", "X_OH"]
case_name = "N_30"
# directories using .format for Python 2
truth_path = "SIM_results/NH3_NP_CF_reduce/{0}.csv".format(case_name)
inp_dir    = "figs/NH3_NP_CF/{0}/{1}/inputs".format(UQ_name, case_name)
err_dir    = "figs/NH3_NP_CF/{0}/{1}/error_outputs".format(UQ_name, case_name)
plot_dir   = "figs/NH3_NP_CF/{0}/{1}/plots".format(UQ_name, case_name)
nin = len(input_species)


for d in (plot_dir,):
    if not os.path.isdir(d):
        os.makedirs(d)

def make_opts(save_dir, figtype=".png", fontsize=14):
    return {
        "savefigs": True, 
        "figtype": figtype,
        "myfont": {"size": fontsize},
        "save_dir": save_dir,
    }

# —————————————————————————————————————————————————————
# 0) Find the grid value at which TRUE HRR is maximal
# —————————————————————————————————————————————————————
# read header to find column indices
with open(truth_path) as f:
    hdr = f.readline().strip().split(',')
grid_i = hdr.index('grid')
hrr_i = hdr.index('HRR')

truth = np.genfromtxt(truth_path, delimiter=',', skip_header=1)
# pick the grid where HRR is maximal
idx_max = np.nanargmax(truth[:, hrr_i])
grid_max = truth[idx_max, grid_i]
print "True HRR peaks at grid =", grid_max

# —————————————————————————————————————————————————————
# 1) Gather one X,f per iteration at the grid_max (using isclose)
# —————————————————————————————————————————————————————
inp_files = sorted(glob.glob(os.path.join(inp_dir,  "uncertainty_inputs_iter*.csv")))
err_files = sorted(glob.glob(os.path.join(err_dir,  "error_uncertainty_outputs_iter*.csv")))

if len(inp_files) != len(err_files):
    raise RuntimeError("Need same count of input+error files")

X_list = []
f_list = []

for inp_path, err_path in zip(inp_files, err_files):
    # --- load input file header + data ---
    with open(inp_path) as f:
        hdr_inp = f.readline().strip().split(',')
    data_inp = np.genfromtxt(inp_path, delimiter=',', skip_header=1)
    gi_inp   = hdr_inp.index('grid')

    # mask on approximate equality
    mask_inp = np.isclose(data_inp[:, gi_inp], grid_max, atol=1e-6)
    if mask_inp.sum() != 1:
        raise RuntimeError(
            "In {} found {} matches for grid≈{:.7f}".format(inp_path, mask_inp.sum(), grid_max)
        )
    row_inp = data_inp[mask_inp, :][0]

    # extract X dynamically from input_species
    cols = [ hdr_inp.index(sp) for sp in input_species ]
    Xi   = row_inp[cols].reshape((1, nin))
    
    # --- load error file header + data ---
    with open(err_path) as f:
        hdr_err = f.readline().strip().split(',')
    data_err = np.genfromtxt(err_path, delimiter=',', skip_header=1)
    gi_err   = hdr_err.index('grid')

    mask_err = np.isclose(data_err[:, gi_err], grid_max, atol=1e-6)
    if mask_err.sum() != 1:
        raise RuntimeError(
            "In {} found {} matches for grid≈{:.7f}".format(err_path, mask_err.sum(), grid_max)
        )
    row_err = data_err[mask_err, :][0]

    # extract f = HRR error
    qi = hdr_err.index('HRR')
    fi = row_err[qi].reshape((1,1))

    X_list.append(Xi)
    f_list.append(fi)

# After gathering all Xi, fi: build X, f
# —————————————————————————————————————————————————————
X = np.vstack(X_list)   # M×6
f = np.vstack(f_list)   # M×1

# drop invalid f (NaN/Inf)
mask = np.isfinite(f).flatten()
if mask.sum() < len(f):
    print "Dropping %d invalid rows" % (len(f) - mask.sum())
X = X[mask, :]
f = f[mask, :]

# —————————————————————————————————————————————————————
# 2) OLS subspace
# —————————————————————————————————————————————————————
ss = ac.subspaces.Subspaces()
ss.compute(X=X, f=f, nboot=1000, sstype='OLS')
opts = make_opts(plot_dir)
print("eigenvecs", ss.eigenvecs[0].reshape(nin, 1))
ac.utils.plotters.eigenvectors(ss.eigenvecs[0].reshape(nin, 1), opts=opts)
ac.utils.plotters.eigenvalues(ss.eigenvals, ss.e_br, opts=opts)
ac.utils.plotters.subspace_errors(ss.sub_br, opts=opts)
ac.utils.plotters.sufficient_summary(X.dot(ss.W1), f.flatten(), opts=opts)
print "Finished OLS subspace plots."

# # —————————————————————————————————————————————————————
# # 3) QPHD subspace
# # —————————————————————————————————————————————————————
# ss = ac.subspaces.Subspaces()
# ss.compute(X=X, f=f, nboot=800, sstype='QPHD')
# opts2 = make_opts(plot_dir)
# opts2['suffix'] = '_qp'
# ac.utils.plotters.eigenvalues(ss.eigenvals, ss.e_br, opts=opts2)
# ac.utils.plotters.subspace_errors(ss.sub_br, opts=opts2)
# ac.utils.plotters.sufficient_summary(X.dot(ss.W1), f.flatten(), opts=opts2)
# print "Finished QPHD subspace plots."


##save data to replot!
# 1) Eigenvalues
eigvals = ss.eigenvals               # shape (n,)
np.savetxt(os.path.join(plot_dir, "eigenvalues.csv"),
           eigvals,
           delimiter=",",
           header="eigenvalue",
           comments="")

# 2) Bootstrap error bars on eigenvalues
e_br = ss.e_br                       # shape (n,2) if you kept lower/upper
np.savetxt(os.path.join(plot_dir, "eigenvalue_errorbars.csv"),
           e_br,
           delimiter=",",
           header="lower,upper",
           comments="")

# 3) Leading eigenvector(s), pure Python/NumPy CSV
W1 = ss.eigenvecs[0].reshape((nin,))   # flatten to length-nin
out_path = os.path.join(plot_dir, "eigenvector_PC1.csv")
with open(out_path, "w") as fh:
    # header line: species,PC1
    fh.write("species,PC1\n")
    for sp, val in zip(input_species, W1):
        fh.write("{},{}\n".format(sp, val))

# 4) Subspace bootstrap errors (if you want them)
sub_err = ss.sub_br                 # shape (n_boot, subspace_dim) or similar
np.savetxt(os.path.join(plot_dir, "subspace_errors.csv"),
           sub_err,
           delimiter=",",
           header="errboot1,errboot2,...",
           comments="")

# 5) Sufficient summary data (projected coords & f), NumPy + CSV
Xproj = X.dot(ss.W1)                  # shape (M, 1) or (M, d)
out_path2 = os.path.join(plot_dir, "sufficient_summary.csv")
with open(out_path2, "w") as fh:
    # header: y_pc1[,y_pc2,...],f_true
    hdr = ["y_pc{}".format(i+1) for i in range(Xproj.shape[1])]
    hdr.append("f_true")
    fh.write(",".join(hdr) + "\n")
    for row_proj, val_f in zip(Xproj, f.flatten()):
        # row_proj could be length >1
        line = ",".join(str(x) for x in row_proj) + "," + str(val_f) + "\n"
        fh.write(line)
        

# 4) Fit & plot a quadratic response surface (Py2), *1D only*
RS = ac.utils.response_surfaces.PolynomialApproximation(2)
# pick only the first PC direction:
w1 = ss.W1[:, 0].reshape(-1, 1)    # (nin,1)
y  = X.dot(w1)                     # (M,1)
RS.train(y, f)                     # now you get 3 poly weights

# make a 1D grid
y_line = np.linspace(y.min(), y.max(), 300)
f_line = RS.predict(y_line.reshape(-1,1))[0]


# # —————————————————————————————————————————————————————
# # 4) Fit & plot a quadratic response surface (Py2)
# # —————————————————————————————————————————————————————
# RS = ac.utils.response_surfaces.PolynomialApproximation(2)
# y  = X.dot(ss.W1)               # (M,1)
# # ensure f is two-dimensional:
# RS.train(y, f)                  # or f.reshape(-1,1)

# project onto first eigenvector
w1 = ss.W1[:, 0].reshape(-1,1)    # (nin,1)
y1 = X.dot(w1)                    # (M,1)
RS = ac.utils.response_surfaces.PolynomialApproximation(2)
RS.train(y1, f)                   # now truly a 1-D fit



# # build smooth curve
# y_line = np.linspace(y.min(), y.max(), 300)
# f_line = RS.predict(y_line.reshape(-1,1))[0]

# plt.figure(figsize=(6,4))
# plt.scatter(y, f.flatten(), alpha=0.6, label="data")
# plt.plot(y_line, f_line, 'r-', linewidth=2, label="2nd-order fit")
# plt.xlabel("Active variable (PC1)")
# plt.ylabel("HRR error")
# plt.title("Response Surface on QPHD Active Variable")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()

# out_fig = os.path.join(plot_dir, "sufficient_summary_qphd_rs.png")
# plt.savefig(out_fig, dpi=300)
# plt.close()

# # 4.5 save fit curve to CSV
# out_csv = os.path.join(plot_dir, "sufficient_summary_qphd_rs.csv")
# import numpy as np
# # stack the two arrays into shape (300,2)
# data = np.column_stack((y_line, f_line))
# # write with header
# np.savetxt(
#     out_csv,
#     data,
#     delimiter=",",
#     header="y_pc1,f_fit",
#     comments=""
# )

# print "Wrote response surface plot to %s" % out_fig

## —————————————————————————————————————————————————————
# 5) Predict‐vs‐True data save, plot & R² (Py2 style)
# —————————————————————————————————————————————————————
# 5.1 Re‐project onto the first PC direction
w1    = ss.W1[:, 0].reshape(-1, 1)    # (nin,1)
y1    = X.dot(w1)                     # (M,1)

# 5.2 Predict (after RS.train on y1,f)
# f_true = f.flatten()                  # (M,)
# f_pred = RS.predict(y1)[0].flatten()  # (M,)
f_true = f.flatten()                  # (M,)
f_pred  = RS.predict(y1)[0].flatten() # now also (M,)

# 5.3 save data to CSV for later replotting
out_csv2 = os.path.join(plot_dir, "predict_vs_true_qphd.csv")
data    = np.column_stack((f_true, f_pred))
np.savetxt(
    out_csv2,
    data,
    delimiter=",",
    header="f_true,f_pred",
    comments=""
)
print "Wrote predict-vs-true data to %s" % out_csv2

# 5.4 compute R² via corrcoef
r  = np.corrcoef(f_true, f_pred)[0,1]
R2 = r**2

# 5.5 make the plot
plt.figure(figsize=(6,4))
plt.scatter(f_true, f_pred,
            alpha=0.6, edgecolor='k', label="points")
mn, mx = min(f_true.min(), f_pred.min()), max(f_true.max(), f_pred.max())
plt.plot([mn, mx], [mn, mx], 'r--', linewidth=1, label="1:1 line")

plt.xlabel("True HRR error")
plt.ylabel("Predicted HRR error")
plt.title("Predict vs True ($R^2$ = %.2f)" % R2)
plt.legend(loc='best', frameon=False)
plt.grid(True)
plt.tight_layout()

out_fig2 = os.path.join(plot_dir, "predict_vs_true_qphd.png")
plt.savefig(out_fig2, dpi=300)
plt.close()

print "Wrote predict-vs-true plot to %s" % out_fig2
print "R² = %.4f" % R2