# A Simplified Model for Linear Mode Coupling in Multimode Fibers
This repo contains the code for the submission "A Simplified Model for Linear Mode Coupling in Multimode Fibers" to the 2026 European Conference on Optical Communications.

Two multisectional models with the following transfer matrix for the $i$-th section are considered:

$$ \mathbf{K}_i \mathbf{R} $$

where $\mathbf{R}$ is a block diagonal matrix with random unitary blocks uniformly distributed over the $\mathrm{U}(N)$ group, whose blocks have the sizes of the mode groups, and

$$ \mathbf{K}_i = e^{j \mathbf{B}_0 + \Delta \beta \, g(\mathrm{XT}, G) \, (\mathbf{P}_i-\mathbf{P}_i^\dagger)}
$$ 

where 

$$ \mathbf{B}_0 $$ 

is the diagonal matrix containing the product between the section length $L$ and the propagation constants $\beta_0^{(n)}$ ($n$ is the mode index) for a graded-index multimode fiber (GIMMF) which are assumed to be constant within a mode group and to differ by $\Delta \beta$ from the neighbouring groups (the term $\Delta \beta \, L$ can also be regarded as an optimization parameter or kept fixed to a certain value, e.g., $10^4$), and 

$$ \mathbf{P}_i $$

is a matrix whose elements are independent complex Gaussian random variables $P_{ab} \sim \mathcal{CN}(0, h_{A-B}^2)$. 

In the paper, two choices for $h_{A-B}$ are considered:

1) $h_{A-B} = 1$, i.e., same variance for all coupling coefficients. We call this the "uniform" choice. For $\mathrm{XT} \le -5\,\mathrm{dB}$, an approximate expression for $g$ is $$ g(\mathrm{XT}, G) \approx \sqrt{\frac{\mathrm{XT}_t - 0.91 + 0.31G}{2M}}$$ where $\mathrm{XT}_t$ is the target crosstalk; 

2) $h_{A-B} = e^{-2.25|A-B|}$. We call this the non-uniform choice. We don't provide an approximate expression for $g$ in this case, but rather a numerical one (see the first dataset listed in the Datasets section below) for $\mathrm{XT} \in [-80, 0] \, \mathrm{dB}$ and $G = 2, 3, \dots, 12$.

The models and the code are for a frequency-independent channel, the inclusion of frequency-dependence effects is immediate and explained in the paper.

## Code
The repo is built around an abstract parent class `MultiSectionalModel` (see [localLibrary/](localLibrary/)), which provides a common framework for simulating mode coupling in few-mode or multimode optical fibres using a transfer-matrix approach. The subclass `expmModel` implements the linear coupling model described above.

Other multisectional models can be implemented as subclasses of `MultiSectionalModel`, inheriting all fundamental functionalities such as Monte-Carlo simulations and computation of crosstalk metrics.

### Example scripts

- [example1_single_section.m](example1_single_section.m): Compute the transfer matrix and crosstalk statistics for a single section of the two models.
- [example2_concatenation.m](example2_concatenation.m): Generate the transfer matrix of a multisectional channel.
- [example3_tuningParameter.m](example3_tuningParameter.m): Obtain the numerical relation $g(\mathrm{XT}, G)$ between the tuning parameter and the target crosstalk.

### Datasets (in `databases/`)

- `single_segment_dB_vs_tuning_parameter_expm_exp_per_mode_a=-2.25_dBeta0=1e4_with_pol_with_igsc.mat`:
  $g(\mathrm{XT}, G)$ for the non-uniform case, $\mathrm{XT} \in [-80, 0]\,\mathrm{dB}$, $G = 2,\dots,12$.
- `single_segment_dB_vs_tuning_parameter_expm_unif_dBeta0=1e4_with_pol_with_igsc.mat`:
  $g(\mathrm{XT}, G)$ for the uniform case, with inter-group coupling.
- `single_segment_dB_vs_tuning_parameter_expm_unif_dBeta0=1e4_with_pol_no_igsc.mat`:
  $g(\mathrm{XT}, G)$ for the uniform case, without inter-group coupling.


## Requirements
Tested on MATLAB R2025b.

## Citation
A BibTeX entry will be added if the paper is accepted. For the moment, the following can be used:

P. Carniello, F. M. Ferreira, F. A. Barbosa, M. J. Li, N. Hanik, "A Simplified Model for Linear Mode Coupling in Multimode Fibers", submitted to 2026 European Conference on Optical Communications (ECOC), 2026
