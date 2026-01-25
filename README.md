# CAS-NoiseMod: Channel-Adaptive Secure Noise Modulation for IoT Networks

This repository contains the MATLAB simulation code for the paper:

**"Channel-Adaptive Secure Noise Modulation for IoT Networks"**  
*Serdar Prencuva, Yusuf Islam Tek, and Ertugrul Basar*  
IEEE Wireless Communications Letters (Manuscript ID: WCL2025-2845)

---

## Overview

CAS-NoiseMod is a physical layer security (PLS) scheme that leverages channel-adaptive pre-scaling to achieve secure wireless communication without cryptographic encryption. The scheme encodes information bits through variance switching of Gaussian waveforms and scales the transmitted signals by the magnitude of the legitimate channel, creating an asymmetric advantage for the legitimate receiver over eavesdroppers.

---

## Repository Contents

| File | Description |
|------|-------------|
| `Complete.m` | Unified MATLAB script that generates all simulation figures |

---

## Requirements

- **MATLAB R2020a** or later (tested on R2023b)
- **Statistics and Machine Learning Toolbox** (for `gamcdf` function)

---

## How to Run

1. Clone or download this repository
2. Open MATLAB and navigate to the repository folder
3. Run the script:

```matlab
Complete
```

The script will execute automatically and generate the following outputs:

---

## Outputs

### Figure 4: Bob's BER Performance
- **Files:** `Fig4_Bob_BER_multi_alpha.pdf`, `Fig4_Bob_BER_multi_alpha.eps`
- **Description:** BER performance of the legitimate receiver (Bob) across different variance ratio (α) values with N = 120, showing both simulation results and theoretical predictions.

### Figure 5: CAS-NoiseMod vs Conventional NoiseMod Comparison
- **Files:** `Fig5_2D_CAS_vs_Baseline.pdf`, `Fig5_2D_CAS_vs_Baseline.eps`, `Fig5_2D_CAS_vs_Baseline.png`
- **Description:** BER comparison between CAS-NoiseMod and conventional NoiseMod (without channel-adaptive pre-scaling) for both Bob and Eve at α = 3.5. This figure demonstrates that the security advantage (Eve's BER near 0.5) is specifically induced by the channel-adaptive mechanism.

---

## Simulation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `N` | 120 | Samples per bit (block length) |
| `P` | 1 | Average transmit power |
| `nbits_F4` | 10⁵ | Number of bits for Figure 4 simulations |
| `nbits_F5` | 10⁶ | Number of bits for Figure 5 simulations |
| `alpha_list_F4` | [1, 1.5, 3, 10] | Variance ratios for Figure 4 |
| `alpha_F5` | 3.5 | Fixed variance ratio for Figure 5 |
| `delta_dB` | -10 to 40 dB | SNR-like parameter range |

---

## Reproducibility

A fixed random seed (`rng(3, 'twister')`) is used to ensure reproducibility of all simulation results.

---

## Key Equations Reference

The code implements the following key equations from the paper:

- **Eq. (2):** Variance values for binary signaling: σ²_L = 2P/(1+α), σ²_H = α·σ²_L
- **Eq. (3):** CAS-NoiseMod pre-scaling: s[n] = ρ·x[n], where ρ = |h|
- **Eq. (5):** Received signal: y[n] = h·s[n] + w[n]
- **Eq. (6):** Received variance (dual-channel effect): σ̃² = ρ⁴·σ²_b + σ²_w
- **Eq. (7):** Decision statistic: sample variance estimate
- **Eq. (8)/(17):** Optimal threshold computation
- **Eq. (21):** BEP averaging over Rayleigh fading

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{prencuva2025casnoismod,
  author  = {Prencuva, Serdar and Tek, Yusuf Islam and Basar, Ertugrul},
  title   = {Channel-Adaptive Secure Noise Modulation for IoT Networks},
  journal = {IEEE Wireless Communications Letters},
  year    = {2025},
  note    = {Manuscript ID: WCL2025-2845}
}
```

---

## Contact

For questions or issues, please contact:

- **Serdar Prencuva** - sprencuva19@ku.edu.tr  
- Department of Electrical and Electronics Engineering, Koç University, Istanbul, Türkiye

---

## Acknowledgments

This work is supported by the Scientific and Technological Research Council of Türkiye (TÜBİTAK) through the 1515 Frontier R&D Laboratories Support Program for the Türk Telekom 6G R&D Lab (Project No. 5249902) and under Grant No. 124E146.

---

## License

This code is provided for academic and research purposes. Please cite the paper if you use this code in your work.
