# A homoclinic route to chaos in omnivore communities
This repository contains code for the paper titled "A homoclinic route to chaos in omnivore communities" authored by Yiyuan Niu, Ju Kang, Wei Tao, Xin Wang.

All codes are written in Matlab.

## Authors:
- Yiyuan Niu<sup>1&dagger;</sup>,
- Ju Kang<sup>2&dagger;&ddagger;</sup> (Corresponding Author),
- Wei Tao<sup>1</sup>,
- Xin Wang<sup>1</sup>* (Corresponding Author)

## Authors' Affiliations:
1. School of Physics, Sun Yat-sen University, Guangzhou 510275. China
2. School of Ecology, Sun Yat-sen University, Shenzhen 518107, China

<sup>&ddagger;</sup> For Correspondence: kangj29@mail.sysu.edu.cn

*For Correspondence: wangxin36@mail.sysu.edu.cn

## Content
```
.
├── code/        MATLAB source code for simulations and analysis
├── data/        Precomputed datasets used to generate manuscript figures
├── LICENSE
└── README.md
```

---

## Code Directory (`code/`)

This folder contains all MATLAB scripts and functions required to reproduce the simulations and figures.

### Main Figure Scripts

These are the entry points used to generate each figure:

| Script                   | Description                     |
| ------------------------ | ------------------------------- |
| `main_eq1_Figure2a.m`    | Simulation for Fig. 2a          |
| `main_eq1_Figure2bd.m`   | Simulation Fig. 2b and 2d       |
| `main_eq1_Figure2c.m`    | Simulation Fig. 2c              |
| `main_eq1_Figure3abcd.m` | Simulation Fig. 3a–d            |
| `main_eq1_Figure3ef.m`   | Simulation Fig. 3e–f            |
| `main_eq1_Figure4abcd.m` | Simulation Fig. 4a–d            |

Run these scripts directly to reproduce the corresponding results.

---

### Core Model Implementation Funtions

| File                                  | Purpose                                                |
| ------------------------------------- | ------------------------------------------------------ |
| `simulate_equation_1.m`               | Numerical integration of Equation (1)                  |
| `compute_rel_abundance.m`             | Post-processing: time-averaged abundance & composition |
| `compute_jac_eig_eq1.m`               | Jacobian and eigenvalue analysis                       |
| `unpack_args.m`                       | Parameter handling utility                             |
| `lyapunov_benettin.m`                 | Lyapunov exponent (Benettin method)                    |
| `lyapunov_benettin_jacobian_method.m` | Jacobian-based Lyapunov calculation                    |
| `findpeaks.m`                         | Peak detection for oscillatory dynamics                |
| `plot_time_series.m`                  | Visualization helper                                   |

---

## Data Directory (`data/`)

Contains simulation outputs used to construct the manuscript figures.
Large datasets are stored in compressed form to reduce repository size.

| File                         | Description                                     |
| ---------------------------- | ----------------------------------------------- |
| `figure_2a.csv`              | Data for Fig. 2a                                |
| `figure_2b.7z`               | Compressed dataset for Fig. 2b                  |
| `figure_2c_trajectory_*.csv` | Trajectories for Fig. 2c                        |
| `figure_2d.7z`               | Compressed dataset for Fig. 2d                  |
| `figure_3ab.7z`              | Data for Fig. 3a–b                              |
| `figure_3cd.7z`              | Data for Fig. 3c–d                              |
| `figure_3e.csv`              | Data for Fig. 3e                                |
| `figure_3f_*.csv`            | Species trajectories for Fig. 3f                |
| `figure_4abc.7z.001–006`     | **Split archive** containing Fig. 4a–c datasets |
| `figure_4e.csv`              | Summary table used in Fig. 4e                   |

---
The Fig. 4 dataset is divided into multiple parts using **7-Zip split compression** to comply with GitHub’s 100 MB file size restriction.

All parts must be present to reconstruct the original file.

Extract these data using 7-Zip, and make sure the following files are together:

```
figure_4abc.7z.001
figure_4abc.7z.002
figure_4abc.7z.003
figure_4abc.7z.004
figure_4abc.7z.005
figure_4abc.7z.006
```
---

## License
This project is licensed under the MIT License 
- see the [LICENSE](./LICENSE) file for details.