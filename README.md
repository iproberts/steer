# About

This repository contains code related to the following paper on STEER, a beam selection methodology for full-duplex millimeter wave (mmWave) wireless communication systems.

[1] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "STEER: Beam Selection for Full-Duplex Millimeter Wave Communication Systems," _IEEE Trans. Commun._, Oct. 2022, [PDF](https://ianproberts.com/pdf/pub/steer.pdf).

Using the code in this repo, which is based on the work presented in [1], users can run STEER in order to conduct research on full-duplex mmWave communication systems.

The work of [1] was inpsired by phenomena observed when analyzing our nearly 6.5 million measurements of self-interference in the following papers.

[2] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Beamformed Self-Interference Measurements at 28 GHz: Spatial Insights and Angular Spread," _IEEE Trans. Wireless Commun._, Nov. 2022, [PDF](https://ianproberts.com/pdf/pub/bfsi.pdf), [GitHub](https://ianproberts.com/bfsi).

[3] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Spatial and Statistical Modeling of Multi-Panel Millimeter Wave Self-Interference," Submitted to _IEEE J. Sel. Areas Commun._, 2023, [PDF](https://ianproberts.com/pdf/pub/simodel.pdf), [GitHub](https://ianproberts.com/simodel).

These measurements of self-interference were taken at 28 GHz in an anechoic chamber using two colocated 256-element phased arrays mounted on separate sides of an equilateral triangular platform. Please see [2] and [3] for details for a summary of these measurements.

If you use this code or our paper in your work, please cite [1] with the following BibTeX.

```
@ARTICLE{roberts_steer_2022,
    author={I. P. Roberts and A. Chopra and T. Novlan and S. Vishwanath and J. G. Andrews},
    title={\textsc{Steer}: Beam Selection for Full-Duplex Millimeter Wave Communication Systems},
    journal={IEEE Trans. Commun.},
    year=2022,
    month={Oct.},
    volume={70},
    number={10},
    pages={6902--6917},
}
```

Related work can be found at https://ianproberts.com.

# What is Self-Interference? 

When a transceiver (a wireless device) attempts to transmit and receive at the same time using the same frequency spectrum, some of its transmitted signal will leak into its receiver, corrupting reception of a desired signal.
This undesired leakage (or coupling) is called "self-interference". 
Transmit and receiving at the same time using the same frequency spectrum is called "full-duplex" operation.

This work is particularly focused on self-interference in full-duplex systems operating at mmWave frequencies (roughly 30 GHz to 100 GHz).
In mmWave systems, dense antenna arrays containing dozens or even hundreds of individual antennas are used to overcome high path loss at these high carrier frequencies. 
They do this by forming very directional beams, focusing energy in a particular direction to increase received signal power. 

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222208767-c359fd9d-0fe0-4814-a56b-46d1f3fd306d.svg"/>
</p>

This work is interested in how much self-interference is coupled in full-duplex mmWave systems when using particular transmit and receive beams. Some transmit and receive beams will couple higher self-interference than others; this depends on the steering direction of the beams and on the (unknown) underlying self-interference channel. 

# What is STEER?

Modern mmWave communication systems rely on beam alignment to deliver sufficient beamforming gain to close the link between devices. 
In [1], we present the first beam selection methodology for full-duplex mmWave systems, called STEER, that delivers high beamforming gain while significantly reducing the self-interference coupled between the transmit and receive beams. 

STEER does not necessitate changes to conventional beam alignment methodologies nor additional over-the-air feedback, making it compatible with existing cellular standards. 
Instead, STEER uses conventional beam alignment to identify the general directions beams should be steered, and then it makes use of a minimal number of self-interference measurements to jointly select transmit and receive beams that deliver high gain in these directions while coupling low self-interference. 

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222922193-67646503-455b-45bf-a573-348fc06e3703.svg"/>
</p>

In [1], STEER was implemented on an actual industry-grade 28 GHz phased array platform, which showed that full-duplex operation with beams selected by STEER can notably outperform both half-duplex and full-duplex operation with beams chosen via conventional beam selection. 
For instance, STEER can reliably reduce self-interference by more than 20 dB and improve SINR by more than 10 dB, compared to conventional beam selection. 
Our experimental results highlight that beam alignment can be used not only to deliver high beamforming gain in full-duplex mmWave systems but also to mitigate self-interference to levels near or below the noise floor, rendering additional self-interference cancellation unnecessary with STEER. 

# How STEER Works

STEER leverages our observations from a recent measurement campaign of mmWave self-interference [2], [3], which showed that small shifts in
the steering directions of the transmit and receive beams (on the order of one degree) can lead to noteworthy reductions in self-interference. 
STEER makes use of self-interference measurements across small spatial neighborhoods to jointly select transmit and receive beams at the full-duplex device that offer reduced self-interference while delivering high beamforming gain on the uplink and downlink.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222922156-2ad6a1c4-55b7-4ee6-9b21-78f32d1b044a.svg"/>
</p>

Put simply, STEER slightly shifts the transmit and receive beams at a full-duplex mmWave transceiver to reduce self-interference---with the goal being that these slight shifts do not prohibitively degrade downlink or uplink quality (i.e., SNR).
By maintaining high SNR and reducing self-interference, high SINRs can be achieved with STEER.

# Model Summary

Our statistical model of self-interference is based on two characteristics observed in our measurements:
1. On a large scale (at a high level), there is a connection between the steering directions of the transmit and receive beams and the degree of self-interference incurred. Broadly speaking, some transmit and receive directions tend to incur high self-interference while others tend to incur low self-interference.
2. On a small scale (within small spatial neighborhoods), the system incurs seemingly random amounts of self-interference. Slightly shifting the transmit and receive steering directions can dramatically alter the degree of self-interference coupled.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/221431446-ae3a8393-2c4b-41e8-a66c-fcdf258f63e4.svg"/>
</p>

We leverage these large-scale and small-scale characteristics to construct a stochastic model of self-interference that both statistically and spatially aligns with our measurements. 
A block diagram summarizing our model is shown above.
For particular transmit and receive beams, a mean parameter is estimated, which dictates the location of the distribution from which self-interference is drawn.
The variance of this distribution is dictated by the mean parameter and other model parameters.
This approach allows our model to capture the large-scale spatial trends in self-interference along with the small-scale variability observed over small spatial neighborhoods.
With appropriate parameterization, our model has the potential to be extended to other systems and environments beyond our own. 

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/221431451-9f7bec04-4659-4b9c-95d5-97deeb3f2345.svg"/>
</p>

To construct our model, we uncovered a coarse geometric approximation of the self-interference channel from within our measurements, which suggests that the dominant coupling between the transmit and receive arrays manifests as clusters of rays in a far-field manner (as illustrated above), rather than in a idealized near-field, spherical-wave fashion.
This is a novel finding that can steer future work aiming to model self-interference MIMO channels in full-duplex mmWave systems.

# Contents

This repo contains the following MATLAB code:
 - a main script `main.m` illustrating example usage
 - an `array` object that can be used to create and interface with arbitrary antenna arrays

# Example Usage

Suppose a full-duplex mmWave base station employs a codebook of transmit beams that it uses to serve downlink and a codebook of receive beams that it uses to serve uplink. 

The degree of self-interference incurred by the base station will depend on the transmit beam and receive beam that it uses to serve downlink and uplink users. Each transmit-receive beam pair will couple a unique amount of self-interference. 

We can draw a statistical realization of this coupling for each beam pair across the codebooks using the script `main.m`, which implements our model in MATLAB.

### Set System and Model Parameters

The first step is to set the desired system and model parameters. By default, the parameters provided in [1, Table II] can be used. 

```
% System parameters
EIRP_dBm = 60; % transmit array EIRP
P_noise_dBm = -68; % receive array integrated noise power (includes amplification)

% Model parameters: location and scale parameters for mu
G_dB = -129.00; % location parameter
xi = 0.502; % scale parameter

% Model parameters: cluster centers and angular spreads for H_bar
aod_az_list = [-174 126 -118 126]; % AoD azimuth
aod_el_list = [0 0 0 0]; % AoD elevation
aoa_az_list = [-122 -122 -122 118]; % AoA azimuth
aoa_el_list = [0 0 0 0]; % AoA elevation
spread_az = 4; % angular spread in azimuth
spread_el = 3; % angular spread in elevation
spread_az_res = 1; % angular spread resolution in azimuth
spread_el_res = 1; % angular spread resolution in elevation

% Model parameters: estimator and variance parameters for sigma^2
alpha = -0.733; % slope
beta = 42.53; % bias
nu_squared = 126.091; % variance
```

### Define Transmit and Receive Codebooks

The transmit and receive codebooks used at the full-duplex base station can be set by defining the steering directions of the codebooks' beams. The steering direction of each beam contains a component in azimuth and elevation (see [1, Fig. 2]). 

In the example below, the transmit and receive codebooks are identical. Each codebook has beams that span in azimuth from -56 deg. to 56 deg. in 8 deg. steps and in elevation from -8 deg. to 8 deg. in 8 deg. steps. This amounts to a total of 45 beams in each codebook, meaning there are 2025 transmit-receive beam pairs.

```
% transmit directions (e.g., transmit codebook)
tx_dir_az = flip([-56:8:56]); % flips are not necessary, just for plotting convenience
tx_dir_el = flip([-8:8:8]);

% receive directions (e.g., receive codebook)
rx_dir_az = tx_dir_az; % assume same for simplicity
rx_dir_el = tx_dir_el;
```

The azimuth-elevation of each transmit and receive steering direction can then be populated as follows.

```
% each transmit az-el pair
tx_dir_az_el_deg = [repelem(tx_dir_az(:),length(tx_dir_el)),repmat(tx_dir_el(:),length(tx_dir_az),1)];
num_tx = length(tx_dir_az_el_deg(:,1));

% each receive az-el pair
rx_dir_az_el_deg = [repelem(rx_dir_az(:),length(rx_dir_el)),repmat(rx_dir_el(:),length(rx_dir_az),1)];
num_rx = length(rx_dir_az_el_deg(:,1));
```


# Questions and Feedback

Feel free to reach out to the corresponding author of [1] with any questions or feedback.

# Acknowledgments

This work has been supported by the National Science Foundation Graduate Research Fellowship Program (Grant No. DGE-1610403). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
