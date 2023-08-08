# Studies of discrete ambiguities in Photoproduction

[![GitHub page](https://img.shields.io/badge/GitHub-README.md-yellowgreen)](https://github.com/mmikhasenko/PhotoAmbiguities.jl)
[![arXiv article](https://img.shields.io/badge/article-%20hep--ph%3A2212.11767-brightgreen)](https://inspirehep.net/literature/2673390)


## Overview
The "PhotoAmbiguities.jl" repository is a Julia-based codebase dedicated to the in-depth analysis of ambiguities in photoproduction, as discussed in the associated paper titled "Ambiguities in Partial Wave Analysis of Two Spinless Meson Photoproduction." While the paper offers an academic treatment of the problem, demonstrating that one solution is always superior to others, challenges arise due to the Barrlett zeros present in the unpolarized case. These zeros can lead to local minima that might, due to statistical fluctuations, exhibit a higher likelihood for finite statistics samples. This repository delves into this effect, identifying local minima on high-statistic samples and subsequently running pseudoexperiments to observe the fluctuation of the likelihood difference in small-statistic samples.

## Features
- **Notebooks**: Jupyter notebooks for interactive analysis and visualization.
- **Plots**: Visual representations of data and results, including the "dminima plot" and "two minimas."
- **Scripts**: Utility scripts for various tasks, including Minuit adjustments.
- **Source Files**: Core Julia code for the main functionalities of the project.

## Features

### Julia Module

- [`PhotoAmbiguities.jl/src`](https://github.com/mmikhasenko/PhotoAmbiguities.jl/tree/master/src) - Core Julia code for the main functionalities of the project.
- [`PhotoAmbiguities.jl/scripts`](https://github.com/mmikhasenko/PhotoAmbiguities.jl/tree/master/scripts) - Utility scripts for various tasks, including Minuit adjustments.
- [`PhotoAmbiguities.jl/plots`](https://github.com/mmikhasenko/PhotoAmbiguities.jl/tree/master/plots) - Visual representations of data and results, including the "dminima plot" and "two minimas."

### Pluto Notebooks (Julia)

- [`SD_waves.jl`](notebooks/SD_waves.html) [[code]](https://github.com/mmikhasenko/PhotoAmbiguities.jl/tree/master/notebooks/SD_waves.jl) Exploration of the simplest model with S and D waves. The polarization value and statistics is tuned in the notebook.
- [`SPD_waves.jl`](notebooks/SPD_waves.html) [[code]](https://github.com/mmikhasenko/PhotoAmbiguities.jl/tree/master/notebooks/SPD_waves.jl) Pseudoexperiments on the wave set with higher waveset.
- [`prototype.jl`](notebooks/prototype.html) [[code]](https://github.com/mmikhasenko/PhotoAmbiguities.jl/tree/master/notebooks/prototype.jl) first implementation of the code in the code in the module. It is kept for educational reason.


## Installation
To use this repository, ensure you have Julia installed. Clone the repository and navigate to its directory:

```bash
git clone https://github.com/mmikhasenko/PhotoAmbiguities.jl.git
cd PhotoAmbiguities.jl
```

## Usage

### Using Pluto

Start [Pluto]() server, open notebook
Open notebook with Pluto

### Julia in VSCode, or terminal

1. Activate the Julia environment:
```julia
using Pkg
Pkg.activate(".")
```

2. Run the desired notebooks or scripts by `include("PATH_TO_FILE")`.

## Contributing
Contributions are welcome! Please fork the repository and create a pull request with your changes.

## Acknowledgments
Special thanks to the colleagues from [JPAC team](https://www.jpac-physics.org/), and the [GlueX experiment](http://www.gluex.org/).
