![GitHub Actions](https://github.com/rikenbit/hmmTensor/actions/workflows/build_test_push.yml/badge.svg)

# hmmTensor
R package for Hidden Markov Model (HMM) by Matrix/Tensor Decomposition

## Installation

~~~~
git clone https://github.com/rikenbit/hmmTensor/
R CMD INSTALL hmmTensor
~~~~
or type the code below in the R console window
~~~~
library(devtools)
devtools::install_github("rikenbit/hmmTensor")
~~~~

## Functions

- **Standard HMM**
  - `Forward`: Forward algorithm with scaling
  - `Backward`: Backward algorithm
  - `Viterbi`: Viterbi decoding (log-space)
  - `BaumWelch`: Baum-Welch EM algorithm (multi-sequence support)
- **Tensor-based HMM**
  - `Seq2Prob`: Convert observation sequences to co-occurrence matrix/tensor (order 2/3, lag, smoothing)
  - `HMM`: HMM parameter estimation via matrix/tensor decomposition (solvers: symNMF, SVD, CP, TT)
- **Data Generation**
  - `toyModel`: Synthetic HMM data (simple, weather, leftright)

## References

- **HMM**: Rabiner, L. R., "A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition", Proceedings of the IEEE, 77(2), 257-286, 1989
- **NMF for HMM**: Kuang, D. et al., "Symmetric Nonnegative Matrix Factorization for Graph Clustering", SDM, 2012
- **Spectral/SVD**: Hsu, D. et al., "A Spectral Algorithm for Learning Hidden Markov Models", Journal of Computer and System Sciences, 78(5), 1460-1480, 2012
- **Tensor Decomposition**: Anandkumar, A. et al., "Tensor Decompositions for Learning Latent Variable Models", JMLR, 15, 2773-2832, 2014
- **Beta-divergence (MU rule)**: Fevotte, C. and Idier, J., "Algorithms for Nonnegative Matrix Factorization with the Beta-Divergence", Neural Computation, 23(9), 2011
- **Beta-divergence (convergence)**: Nakano, M. et al., "Convergence-guaranteed Multiplicative Algorithms for Nonnegative Matrix Factorization with Beta-divergence", IEEE MLSP, 283-288, 2010

## Contributing

If you have suggestions for how `hmmTensor` could be improved, or want to report a bug, open an issue! We'd love all and any contributions.

For more, check out the [Contributing Guide](CONTRIBUTING.md).

## Authors
- Koki Tsuyuzaki
