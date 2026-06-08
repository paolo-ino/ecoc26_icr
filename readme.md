# A Simplified Model for Linear Mode Coupling in Multimode Fibers
This repo contains the code for the submission "A Simplified Model for Linear Mode Coupling in Multimode Fibers" to the 2026 European Conference on Optical Communications. 

Two multisectional models with the following transfer matrix for the $i$-th section are considered:

$$ \mathbf{K}_i \mathbf{R} $$

where $\mathbf{R}$ is a block diagonal matrix with random unitary blocks uniformly distributed over the $\mathrm{U}(N)$ group, whose blocks have the sizes of the mode groups, and

$$ \mathbf{K}_i = e^{j \mathbf{B}_0 + \Delta \beta \, g(\mathrm{XT}, G) \, (\mathbf{P}_i-\mathbf{P}_i^\dagger)}
$$ 

where $\mathbf{B}_0$ is the diagonal matrix of propagation constants $\beta_0^{(n)}$ ($n$ is the mode index) for a graded-index multimode fiber (GIMMF) which are assumed to be constant within a mode group, $ \Delta \beta $ is the mismatch between the propagation constants of two consecutive groups and can be computed analytically or fixed (e.g., we set it to $ 10^4 $), $ \mathbf{P}_i $ is a matrix whose elements are independent complex Gaussian random variables $ P_{ab} \sim \mathcal{CN}(0, h_{A-B}^2) $. 

In the paper, two choices for $h_{A-B}$ are considered:

1) $ h_{A-B} = 1 $, i.e., same variance for all coupling coefficients. We call this the "uniform" choice. For $\mathrm{XT} < -4\,\mathrm{dB} $, an approximate expression for $g$ is $$ g(\mathrm{XT}, G) \approx \sqrt{\frac{\mathrm{XT}_t - 0.91 + 0.31G}{2M}}$$ where $\mathrm{XT}_t$ is the target crosstalk; 

2) $ h_{A-B} = e^{-2.25|A-B|} $. We call this the non-uniform choice. We don't provide an approximate expression for $g$ in this case, but rather a numerical one in the dataset "__________"

## Code
The main script is "concatenation.m", which allows to generate the transfer matrix of a multisectional channel based on one of the two aforementioned models. The code is for the frequency-independent, the inclusion of frequency-dependence effects is immediate and explained in the paper.

Toggle "variance_style" between "uniform" and "exp_perMode" to choose one of the two models.

The other script in the repo is "tuningParameter.m", which allows to obtain the numerically expression for $g(\mathrm{XT}, G)$. The dataset "__________" already contains $g(\mathrm{XT}, G)$  for the non-uniform case for $\mathrm{XT} \in [-80, 0] \, \mathrm{dB}$ and $G = 2, 3, \dots, 12$.