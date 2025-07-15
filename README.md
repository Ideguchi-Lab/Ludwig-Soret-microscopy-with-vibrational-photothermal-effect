# Ludwig-Soret-microscopy-with-vibrational-photothermal-effect

This repository includes source codes used in the paper;

Toda, Keiichiro and Takuro Ideguchi. 2025. “Ludwig-Soret microscopy with vibrational photothermal effect” arXiv [Physics.Bio-Ph]. arXiv. https://arxiv.org/pdf/2502.04578.


## Overview

We can access the results in Fig. 3, 4 and 5 in the main manuscript and related results.
The instructions are included in each directory(see `how_to.txt`).

## Software Requirements

Parentheses indicates standard installation time. This will depend on the user environment and experience.

### Python(30 min)

We recommend one to use virtual environment to run this codebase.

You can install necessary libraries via

```sh
(YOUR_VIRTUAL_ENV) pip install -r requirements.txt
```

Separately, one should install `cupy` following the [official guide](https://docs.cupy.dev/en/stable/install.html). The auther used `cupy==13.0`. One might have to take care of the CUDA driver if already installed.

### Visual Studio Code(10 min)

One needs `Python` extension to run Python code in Jupyter display. Please ensure `# %%` block is activated in the python code.

### Matlab(3 hour)

Please follow the [official guide](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html).

See the `Test Environment` section below for necessary tools.

## How to prepare data

### Data for Fig.2

The experimental data is included in `water1.zip`, so please unzip this.

```sh
unzip water1.zip
```

One has to edit a `path` variable in the `fig2/thermal_diffusivity_experiment_water.py#L22` to the directory where the data is extracted.

### Data for Fig.3

The experimental data is included in `sample1.zip`, `bg.zip`,`av.zip` so please unzip this.

One has to edit `path` variables in the `dual_pump_v2_small.m`, `impulsive_OFF_small.m`, `original_refractive_index_small.m`, `rapid_temperature_rise_by_continuous_heating_small.m`, and `ref_v2_small.m` to the directory where the data is extracted.

## Test Environment

### Python code

- Python 3.11.0, (numpy 1.26.4, cupy 13.0.0)
- Ubuntu 22.04.4 LTS
- CPU: 13th Gen Intel(R) Core(TM) i9-13900K
- GPU: NVIDIA Corporation GV100 [TITAN V]（CUDA after 11.8）

### Matlab code

- Matlab 2019a
- Windows 10 Pro for Workstations (64 bit), version 22H2
- Intel(R) Xeon(R) CPU E3-1240 v5
- MATLAB                                                version 9.6           (R2019a)
- Simulink                                              version 9.3           (R2019a)
- Control System Toolbox                                version 10.6          (R2019a)
- Curve Fitting Toolbox                                 version 3.5.9         (R2019a)
- DSP System Toolbox                                    version 9.8           (R2019a)
- Data Acquisition Toolbox                              version 4.0           (R2019a)
- Image Acquisition Toolbox                             version 6.0           (R2019a)
- Image Processing Toolbox                              version 10.4          (R2019a)
- Instrument Control Toolbox                            version 4.0           (R2019a)
- MATLAB Compiler                                       version 7.0.1         (R2019a)
- MATLAB Compiler SDK                                   version 6.6.1         (R2019a)
- Mapping Toolbox                                       version 4.8           (R2019a)
- Optimization Toolbox                                  version 8.3           (R2019a)
- Signal Processing Toolbox                             version 8.2           (R2019a)
- Simulink Control Design                               version 5.3           (R2019a)
- Statistics and Machine Learning Toolbox               version 11.5          (R2019a)
- Symbolic Math Toolbox                                 version 8.3           (R2019a)
