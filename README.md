# 🧠 Open Field Test (OFT) Analysis – TrAQ-Compatible MATLAB Pipeline

A **robust, publication-ready MATLAB script** for analyzing **Open Field Test (OFT)** experiments using **TrAQ tracking outputs**.
This pipeline provides **locomotion metrics**, **central vs peripheral zone analysis**, and **high-quality visualizations**, while respecting **manually defined peripheral width** and using a **speed calculation identical to the validated EPM pipeline**.

---

## 📌 Key Features

* ✅ **TrAQ-compatible input handling**
* ✅ **Manual peripheral width (cm) respected exactly**
* ✅ **Central vs peripheral zone analysis**
* ✅ **Locomotion metrics (distance, speed, mobility)**
* ✅ **Speed calculation identical to EPM analysis**

  * Hard cap: **60 cm/s**
  * Smoothing: **3-point moving mean**
* ✅ **Exact time conservation**
* ✅ **Robust handling of coordinate formats**
* ✅ **Paper-ready visualizations**
* ✅ **Extensive diagnostics & safeguards**

---

## 📁 Repository Structure

```text
.
├── TrAQ_OpenField_Analysis.m
├── README.md
├── example_data/
│   └── Out_Stim2.mat
└── figures/
    ├── *_locomotion_overview.png
    ├── *_spatial_zone_analysis.png
    └── *_summary.csv
```

---

## 🎥 Supported Input

This script expects a **TrAQ output file** (e.g. `Out_Stim2.mat`) containing a `track` structure.

### Supported position fields (auto-detected)

* `track.Centroid` **(preferred)**
* `track.Head`
* `track.Body_Centroid`
* `track.Tail`
* `track.Position`

The script automatically:

* Detects coordinate format (`Nx2`, `2xN`, or `Nx1x2`)
* Removes invalid (NaN) frames
* Truncates safely if frame counts mismatch

---

## 🧭 Arena & Zone Definition

### Arena

* Square open field
* Real size defined by:

  ```matlab
  arena_cm = 70;   % arena side length (cm)
  ```

### Zones

* **Peripheral zone**: user-defined border width (cm)
* **Central zone**: remaining inner area

```matlab
peripheral_width_cm = 10;  % manual, explicit, respected
```

⚠️ If the chosen peripheral width makes the central zone invalid, the script:

* Warns the user
* Falls back to a minimal valid central zone

---

## 🏃 Locomotion Metrics

Computed per frame and summarized across the session:

* Total distance traveled (cm)
* Instantaneous speed (cm/s)
* Mean speed (overall & moving only)
* Time moving vs immobile
* Movement percentage

### Speed calculation (validated)

Identical to the EPM pipeline:

* Frame-to-frame displacement
* Pixel → cm conversion
* Hard cap at **60 cm/s**
* **3-point moving mean smoothing**

This ensures **cross-test comparability** between OFT and EPM.

---

## 📊 Zone-Specific Metrics

For **central** and **peripheral** zones separately:

* Time spent (s, %)
* Distance traveled (cm, %)
* Mean speed (cm/s)

These metrics are commonly used as:

* Anxiety-related measures (center avoidance)
* Exploratory behavior indicators

---

## 🚀 How to Run

1. Place your TrAQ output file (`Out_Stim2.mat`) in the working directory
2. Open MATLAB
3. Adjust parameters at the top of the script if needed:

   ```matlab
   arena_cm = 70;
   peripheral_width_cm = 10;
   ```
4. Run:

   ```matlab
   TrAQ_OpenField_Analysis
   ```

---

## 📤 Outputs

### Console

* Detailed diagnostics
* Locomotion summary
* Zone-specific statistics

### Files (auto-saved)

* **Summary CSV**
* **Per-frame detailed CSV**
* **High-quality figures**

  * Locomotion overview
  * Speed distributions
  * Trajectory with zone boundaries
  * Zone occupancy over time
* Export formats:

  * `.fig`
  * `.png / .bmp`
  * **Vector `.pdf` for publications**

---

## 📈 Example Metrics

| Metric                   | Example |
| ------------------------ | ------- |
| TotalTime_s              | 600     |
| TotalDistance_cm         | 3200    |
| TimeCentral_s            | 85      |
| TimePeripheral_s         | 515     |
| MeanSpeedCentral_cm_s    | 3.1     |
| MeanSpeedPeripheral_cm_s | 5.4     |

---

## 🧪 Validation & Fixes Applied

This version includes multiple **explicit fixes and safeguards**:

* ✔ Peripheral width strictly enforced (no silent overrides)
* ✔ Central zone validity checks
* ✔ Speed calculation synchronized with EPM analysis
* ✔ Robust handling of frame/time mismatches
* ✔ Identical movement thresholds across tests
* ✔ No double-counting or time leakage

---

## 📝 Recommended Methods Text (for papers)

> “Open Field Test analysis was performed using TrAQ tracking outputs. Locomotor speed was calculated from frame-to-frame displacement, converted to cm/s, capped at 60 cm/s, and smoothed using a 3-frame moving average. Central and peripheral zones were defined using a manually specified peripheral border width. Time spent, distance traveled, and mean speed were computed separately for each zone.”

---

## 🧑‍🔬 Intended Use

* Behavioral neuroscience
* Anxiety-related behavioral phenotyping
* Pharmacological studies
* Cross-paradigm comparison with EPM results
* TrAQ-based pipelines requiring transparency and reproducibility

---

## 📜 License

MIT License (recommended)
Free for academic and commercial use with attribution.

---

## 📬 Contact

For questions, extensions, or batch-processing versions:

* Open a GitHub issue
* Or contact the author directly

---

### ✅ Final Note

This script is intentionally **conservative and explicit**.
If results appear lower or stricter than quick scripts, this reflects **correct handling of zones, speed, and time**, not under-counting.

---

If you want, I can also:

* Add a **comparison table vs EPM**
* Create a **batch-analysis version**
* Add **group-level statistics**
* Write a **Methods supplement** combining OFT + EPM

Just tell me.
