# Note On Matrix Functions

## Overview

This repository accompanies my note on matrix functions. It is home for implementations focusing on matrix functions—a critical component of numerical linear algebra with broad applications from solving ordinary and partial differential equations to control theory and machine learning. We extend the notion of scalar functions to act on matrices and tackle both theoretical and computational challenges involved in such extension. This work provides rigorous definitions, algorithmic approaches, and empirical validations for computing matrix functions.

## Table of Contents

1. [Introduction](#introduction)
2. [Features](#features)
3. [Getting Started](#getting-started)
4. [Usage](#usage)
5. [Documentation](#documentation)
6. [Numerical Experiments](#numerical-experiments)
7. [Contributing](#contributing)
8. [Citations](#citations)

## Introduction

This repository is structured around a comprehensive document that aims to:

1. Establish a formal mathematical foundation for matrix functions.
2. Explore the underlying spectral decomposition for evaluating differentiable matrix functions.
3. Address computational challenges, with a focus on the role of Krylov subspaces.
4. Validate theoretical claims through numerical experiments.

For more details, refer to the [document](./paper/matrix-functions.pdf).

## Features

1. **Rigorous Definitions**: Different approaches such as natural definition, spectrum-based definition, and interpolation-based definition are provided.
2. **Computational Efficiency**: Detailed examination of matrix-vector products and their efficient evaluation using Krylov subspaces.
3. **Algorithmic Stance**: Dense and sparse matrix algorithms, as well as two versions of Arnoldi Method.
4. **Numerical Experiments**: Code to reproduce the empirical validation of algorithmic choices and computational efficiency.

## Getting Started

To clone this repository, you can run:

```bash
git clone https://github.com/nathan-rousselot/matrix-functions-krylov.git
```

## Usage

After cloning the repository, you can run the numerical experiments by executing the demonstrations file in [/src](./src/). You can try the functions on your own matrices, and can explore various kind of structure by navigating in Matrix Market.

## Documentation

For detailed mathematical proofs, algorithmic explanations, and empirical results, refer to [document](./paper/matrix-functions.pdf).

## Numerical Experiments

The `numerical_experiments/` directory contains Python notebooks and scripts to run experiments which validate the theoretical results. These are critical for examining the practical efficacy of the computational strategies discussed.


## Citations

If you use this research or codebase in your work, please cite this repository.

```bibtex
@article{rousselot2023matrix,
  title={Note on Matrix Functions},
  author={Nathan Rousselot},
  year={2023}
}
```

---

This repository provides a focused exploration of matrix functions, offering both theoretical explanations and practical implementations. It aims to serve students, researchers, and professionals interested in numerical linear algebra. The code provided in this repository is for educational and research purposes only. It comes with no warranty or guarantee of any kind, either expressed or implied. Users are advised to use the code at their own risk and discretion. I am not responsible for any consequences arising from the use of this software.
