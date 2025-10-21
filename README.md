# Non-Binary Polar Decoder (Frequency-Domain Processing)

This MATLAB project implements a **non-binary polar decoder** over GF(2^p) using **frequency-domain processing** techniques.  
The decoder alternates **Check Node (CN)** and **Variable Node (VN)** operations between frequency and probability domains:  
- CN updates are computed efficiently in the **frequency domain** using the **Fast Walsh–Hadamard Transform (FWHT)**.  
- VN updates use **precomputed Hadamard-domain tables** to accelerate message propagation without repeatedly applying FWHT.  

This hybrid approach reduces decoding complexity while preserving accuracy for large Galois Field sizes.

---

## 🧩 Project Overview

The main script **`NB_polar_CCSK.m`** simulates the complete transmission and decoding chain for CCSK-modulated non-binary polar codes:
1. **Encoding** of random information symbols using `encode.m`.
2. **CCSK modulation** over GF(q) (where q = 2^p).
3. **AWGN channel** simulation at given SNR values.
4. **Soft-input decoding** using the proposed **frequency-domain non-binary polar decoder**.
5. **Frame Error Rate (FER)** computation over multiple Monte Carlo runs.

---

## ⚙️ Decoder Description (`polar_nb_dec.m`)

The function `polar_nb_dec` performs decoding over GF(2^p) using probabilistic or LLR-based input messages.

### Key Features
- **Frequency-domain CN processing:**  
  Check Nodes perform convolutions in GF(2^p) via the FWHT.
- **Table-based VN processing:**  
  Variable Nodes reuse precomputed Hadamard-domain lookup tables to avoid redundant FWHT operations.
- **Support for both LLR and probability inputs.**
- **Frozen-symbol assignment:**  
  Follows the least-to-most reliable channel order as specified by the reliability sequence.

### Inputs
| Parameter | Description |
|------------|-------------|
| `prb1` | [q × N] matrix of symbol probabilities or LLRs |
| `Hadamard` | [q × q] normalized Hadamard matrix (FWHT lookup) |
| `reliab_seq` | [1 × N] channel reliability sequence (least → most reliable) |
| `frozen_symbols` | Sequence of known symbols assigned to the **least reliable channels** on the encoder input side. These symbols must be known by the decoder to be able to decode. |
| `is_LLR` | `true` if `prb1` contains LLRs, otherwise `false` for probabilities |

### Output
| Parameter | Description |
|------------|-------------|
| `decw` | [1 × N] decoded word over GF(2^p) |

---

## 🧮 Data Description (`data.rar`)

The file **`data.rar`** contains **reliability sequences** and related data for multiple configurations.

**Contents:**
- Code lengths **N = 2, 4, 8, ..., 512**
- Field sizes **GF(8), GF(16), GF(32), GF(64), GF(128), GF(256)**
- A **diverse set of SNR values**, specific to each configuration.  

Each text file inside contains:
1. Reliability sequence based on channel **entropies** (least → most reliable)
2. Reliability sequence based on channel **error probabilities**
3. **Channel entropies**
4. **Channel error probabilities**
5. Number of Monte Carlo samples (`N_obs`)

⚠️ **Note:**  
Reliability sequences are **sorted from least to most reliable channels**.  
That means the **first element corresponds to the least reliable** sub-channel.

---

## 📂 Repository Structure


---

## 📜 License & Attribution

This project is released under the **CeCILL-B license** (a French open-source license compatible with the GNU GPL).  
See the full text here: [CeCILL-B License](https://cecill.info/licences/Licence_CeCILL-B_V1-en.html).

**Author:** Abdallah Abdallah  
**Affiliation:** Lab-STICC, UMR 6285, Université Bretagne Sud  
**Funding:** ANR Project MIOT (Grant ANR-24-CE93-0017)  
**Website:** [https://project.inria.fr/miot/](https://project.inria.fr/miot/)

