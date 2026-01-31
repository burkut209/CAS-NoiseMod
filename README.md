# CAS-NoiseMod (Channel-Adaptive Secure Noise Modulation)

This repository contains a **reproducible MATLAB implementation** of **CAS-NoiseMod** (Channel-Adaptive Secure Noise Modulation), a physical-layer security scheme designed for low-complexity **IoT networks**.

The main idea:
- Information bits are carried by **variance switching** of Gaussian waveforms (NoiseMod concept).
- Before transmission, the waveform is **pre-scaled by the legitimate channel magnitude** (channel-adaptive pre-scaling).
- The legitimate receiver benefits from channel knowledge/reciprocity, while an eavesdropper suffers from **mismatched detection** and approaches **near-random guessing BER** in many regimes.

This repo is designed to reproduce key paper-style results (Figure 4 parts a–c) using a single unified MATLAB script.

---

## Contents

- `CAS_NoiseMod.m`  
  A **single unified MATLAB script** that:
  1) Simulates **CAS-NoiseMod** BER for Bob (legitimate receiver)  
  2) Simulates **CAS-NoiseMod** BER for Eve (eavesdropper) with **mismatched threshold**  
  3) Computes **Bob theoretical BEP** using a **moment-matched Gamma approximation** + averaging over Rayleigh fading  
  4) Optionally simulates the **baseline conventional NoiseMod** (no pre-scaling)  
  5) Generates and saves figures consistent with the paper’s Figure 4:
     - Fig. 4(a): Bob BER vs δ for multiple α (simulation + theory)
     - Fig. 4(b): 3D surface of Bob vs Eve BER over (α, δ)
     - Fig. 4(c): 2D comparison of CAS-NoiseMod vs baseline NoiseMod (Bob & Eve)

- `CAS_NOISEMOD (10).pdf`  
  The paper: **"Channel-Adaptive Secure Noise Modulation for IoT Networks"** (IEEE Communications Letters submission/manuscript).

---

## Requirements

### Software
- MATLAB (any reasonably recent version should work)

### Toolboxes
The script uses:
- `gamcdf(...)` → typically requires **Statistics and Machine Learning Toolbox**
- `integral(...)` → base MATLAB (usually available)

If you do not have the Stats toolbox, the script will fail at the theoretical BEP part (Gamma CDF). You can still run the Monte-Carlo simulations.

---

## Quick Start (How to Run)

1. Clone the repository:
   ```bash
   git clone <your-repo-url>
   cd <your-repo-folder>
