# Channel-Adaptive Secure Noise Modulation for IoT Networks (CAS-NoiseMod)

## Article

**Title**  
**Channel-Adaptive Secure Noise Modulation for IoT Networks**

**Authors**  
- Serdar Prencuva  
- Yusuf Islam Tek  
- Ertugrul Basar  

**Published in**  
**IEEE Communications Letters**

**What the article proposes (simple but detailed)**  
The article introduces **Channel-Adaptive Secure Noise Modulation (CAS-NoiseMod)**, a physical-layer security method designed for **low-complexity IoT devices**.

The scheme combines two ideas:

1) **Noise Modulation (NoiseMod): information is carried by variance switching**  
Each transmitted bit selects one of two Gaussian waveform variances:
- Bit **0** → low variance  \( \sigma_L^2 \)  
- Bit **1** → high variance \( \sigma_H^2 \)

A key design knob is:
- \( \alpha = \sigma_H^2 / \sigma_L^2 \)  (how strongly the two bit-classes are separated)

2) **Channel-adaptive pre-scaling: the waveform is scaled using the legitimate channel magnitude before transmission**  
Let the legitimate channel be \( h \) (Rayleigh fading), and let \( \rho = |h| \).  
CAS-NoiseMod pre-scales the NoiseMod waveform:
- \( s[n] = \rho \, x[n] \)

This creates a **double-channel effect** at the legitimate receiver: after propagation through the same channel, the received statistics depend strongly on the legitimate channel realization. The legitimate receiver benefits because the pre-scaling is aligned with the legitimate link, while a **passive eavesdropper** (with an independent channel) observes signals whose statistics do not match the correct decision rule. As a result, the eavesdropper’s detection performance can degrade toward **near-random guessing** (BER close to 0.5) across a wide parameter range.

The article also provides a **theoretical BEP (bit error probability)** analysis for the legitimate receiver under Rayleigh fading using a **moment-matched Gamma approximation** and channel averaging.

---

## Code

This repository provides a unified MATLAB implementation in:

- **`CAS_NoiseMod.m`**

The script reproduces the core numerical results for CAS-NoiseMod, including:
- **Monte-Carlo BER simulations**
- **Closed-form-style theoretical BEP evaluation (via Gamma approximation + fading averaging)**
- **Comparison against conventional NoiseMod (no pre-scaling)**

### What the code is doing (step-by-step)

#### 1) Global parameters and definitions
The script sets:
- `N` = samples per bit (default `120`)
- `P` = average transmit power (default `1`)
- `alpha` values (variance ratio sweep)
- `delta` values in dB (noise robustness sweep)

It converts the paper-style variance definitions:
- \( \sigma_L^2 = \frac{2P}{1+\alpha} \)
- \( \sigma_H^2 = \alpha \sigma_L^2 \)

It also uses the paper’s “delta” parameter by setting receiver noise variance as:
- \( \sigma_w^2 = \sigma_L^2 / \delta \)

So:
- larger \( \delta \) → smaller noise variance → lower BER

---

#### 2) CAS-NoiseMod simulation (Bob, and Eve with mismatched threshold)
For each bit:
1) Draw channels (Rayleigh fading):
- \( h \sim \mathcal{CN}(0,1) \) for Bob  
- \( h_e \sim \mathcal{CN}(0,1) \) for Eve (independent)

2) Generate a **real Gaussian** waveform for the chosen bit:
- bit 0: \( x[n] \sim \mathcal{N}(0,\sigma_L^2) \)
- bit 1: \( x[n] \sim \mathcal{N}(0,\sigma_H^2) \)

3) Apply **channel-adaptive pre-scaling**:
- \( \rho = |h| \)
- \( s[n] = \rho \, x[n] \)

4) Add complex AWGN at each receiver:
- \( w_B[n] \sim \mathcal{CN}(0,\sigma_w^2) \)
- \( w_E[n] \sim \mathcal{CN}(0,\sigma_w^2) \)

5) Form received signals:
- Bob: \( y_B[n] = h\,s[n] + w_B[n] \)
- Eve: \( y_E[n] = h_e\,s[n] + w_E[n] \)

6) Compute the detection statistic (energy / sample second moment):
- \( \hat{S} = \frac{1}{N}\sum_{n=1}^{N} |y[n]|^2 \)

7) Decide the bit by comparing \( \hat{S} \) to a threshold.

**Bob (matched threshold under CAS):**  
Under CAS, Bob’s conditional received variances become:
- \( C_0 = \rho^4 \sigma_L^2 + \sigma_w^2 \)
- \( C_1 = \rho^4 \sigma_H^2 + \sigma_w^2 \)

The script uses the standard log-likelihood-based variance threshold:
- \( \gamma = \frac{C_1 C_0}{C_1 - C_0}\,\ln\!\left(\frac{C_1}{C_0}\right) \)
- Decision: \( \hat{b} = \mathbb{1}[\hat{S} > \gamma] \)

**Eve (mismatched threshold under CAS):**  
Eve does *not* have the correct pre-scaling alignment with the legitimate channel magnitude \( \rho = |h| \). The script models this by letting Eve build her threshold using **her own** channel magnitude \( \rho_e = |h_e| \), i.e., she uses:
- \( C_{0,e}^{(assumed)} = \rho_e^4 \sigma_L^2 + \sigma_w^2 \)
- \( C_{1,e}^{(assumed)} = \rho_e^4 \sigma_H^2 + \sigma_w^2 \)
- \( \gamma_e^{(mis)} = \frac{C_{1,e}^{(assumed)} C_{0,e}^{(assumed)}}{C_{1,e}^{(assumed)} - C_{0,e}^{(assumed)}}\ln\!\left(\frac{C_{1,e}^{(assumed)}}{C_{0,e}^{(assumed)}}\right) \)

Then Eve decides:
- \( \hat{b}_e = \mathbb{1}[\hat{S}_E > \gamma_e^{(mis)}] \)

This mismatch is what creates the security gap in the simulations.

---

#### 3) Theoretical BEP for Bob (moment-matched Gamma + fading averaging)
The script includes a function:

- `BEP_Bob_CAS_Theory(delta_dB_vec, N, alpha, P)`

It evaluates Bob’s BEP by averaging over:
- \( r = |h|^2 \sim \text{Exp}(1) \)

Under CAS, Bob’s conditional variances (as functions of \( r \)) are written in the code as:
- \( C_0(r) = r^2 \sigma_L^2 + \sigma_w^2 \)
- \( C_1(r) = r^2 \sigma_H^2 + \sigma_w^2 \)

The distribution of the statistic \( \hat{S} \) is approximated with a **Gamma distribution** whose parameters are **moment-matched** using an “impropriety factor” term \( \zeta \). Then the error probability is computed using Gamma CDF evaluations at the threshold \( \gamma(r) \), and finally averaged over the exponential fading density \( e^{-r} \) via numerical integration.

This provides the “theory” curves for Bob that the script plots alongside the Monte-Carlo results.

---

#### 4) Conventional NoiseMod baseline (no pre-scaling)
For comparison, the script also simulates the conventional (non-adaptive) NoiseMod baseline:
- transmit: \( s[n] = x[n] \)
- receive: \( y[n] = h\,x[n] + w[n] \)

In this baseline, both Bob and Eve can use **matched thresholds** based on their own channels:
- \( C_0 = |h|^2\sigma_L^2 + \sigma_w^2 \)
- \( C_1 = |h|^2\sigma_H^2 + \sigma_w^2 \)
- same threshold form \( \gamma = \frac{C_1 C_0}{C_1 - C_0}\ln(C_1/C_0) \)

This baseline is used to show that **without** channel-adaptive pre-scaling, the security advantage largely disappears.

---

### What results the script produces (conceptually)
The script generates three paper-aligned result sets:

1) **Bob BER vs \( \delta \)** for several \( \alpha \) values (simulation + theory)  
2) **3D surfaces** of Bob and Eve BER across \( (\alpha,\delta) \) (shows the security gap region)  
3) **2D comparison** at a fixed \( \alpha \): CAS-NoiseMod vs conventional NoiseMod (Bob and Eve)
