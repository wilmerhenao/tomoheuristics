# Set definitions
param numvoxels integer > 0;
param numloops := 5;

set PROJECTIONS = {0..numloops * 178 - 1};
set PROJECTIONSM1 = {0..numloops * 178 - 2};
set LEAVES = {0..79};
set VOXELS;
set KNJPARAMETERS within {n in LEAVES, k in PROJECTIONS, j in VOXELS};

# Parameters
param D {KNJPARAMETERS} >= 0;
param Mbar = numloops * 10;
param U = 10;
param thethreshold {VOXELS} >= 0;
param quadHelperOver {VOXELS} >= 0;
param quadHelperUnder {VOXELS} >= 0;

# Variables
var betas {n in LEAVES, k in PROJECTIONS} binary;
var mu{n in LEAVES, k in PROJECTIONSM1} >= 0;
var xi{n in LEAVES, k in PROJECTIONS} >= 0;
var y{k in PROJECTIONS} >= 0;
var z {j in VOXELS} >= 0;
var z_plus {j in VOXELS} >= 0;
var z_minus {j in VOXELS} >= 0;

# Objective
minimize Total_Impact: sum {j in VOXELS} (quadHelperUnder[j] * z_minus[j] * z_minus[j] + quadHelperOver[j] * z_plus[j] * z_plus[j]);

# Constraints
subject to doses_to_j {j in VOXELS}: z[j] = sum{thisloop in 0..(numloops-1)} sum{ (n, k, j) in KNJPARAMETERS}( D[n, k, j] * xi[n, k + 178 * thisloop]);# + sum{ (n, k, j) in KNJPARAMETERS}( D[n, k, j] * xi[n, k + 178]);
positive_only {j in VOXELS}: z_plus[j] - z_minus[j] = z[j] - thethreshold[j];
Mlimits {n in LEAVES}: sum{k in PROJECTIONSM1} mu[n,k] <= Mbar;
abs_greater {n in LEAVES, k in PROJECTIONSM1}: mu[n,k] >= betas[n, k+1] - betas[n,k];
abs_smaller {n in LEAVES, k in PROJECTIONSM1}: mu[n,k] >= -(betas[n, k+1] - betas[n,k]);
xi_U {n in LEAVES, k in PROJECTIONS}: xi[n,k] <= betas[n,k] * U;
xi_l_y {n in LEAVES, k in PROJECTIONS}: xi[n,k] <= y[k];
xi_g_y {n in LEAVES, k in PROJECTIONS}: xi[n,k] >= y[k] - (1 - betas[n,k]) * U;
###########
options solver gurobi;
#display Total_Impact;
