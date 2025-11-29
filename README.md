# Robust TrAQ Open Field Analysis

This project provides a professional-grade MATLAB script for the analysis of open field test data. It is designed to be robust, configurable, and easy to use, offering a comprehensive analysis of locomotion and zone-based metrics.

## Description

The main script, `professional_TrAQ_OFT.m`, takes tracking data from a MATLAB `.mat` file and performs a detailed analysis, including:

-   **Speed and Distance Calculation:** Computes the speed and total distance traveled by the subject.
-   **Movement Analysis:** Differentiates between periods of movement and immobility based on a configurable speed threshold.
-   **Zone Analysis:** Divides the arena into central and peripheral zones and calculates the time spent and distance traveled in each.
-   **Visualization:** Generates detailed plots for speed analysis, zone analysis, and the subject's trajectory.
-   **Data Export:** Saves the analysis results to CSV files for further processing.

This script is a refactored and professionalized version of the original `TrAQ_OFT_Compelete.m` script.

## Usage

To run the analysis, simply call the `professional_TrAQ_OFT` function from the MATLAB command window. You can run it with default parameters or provide your own.

### Default Usage

```matlab
professional_TrAQ_OFT
```

### Custom Usage

You can customize the analysis by providing key-value pairs for the parameters:

```matlab
professional_TrAQ_OFT('dataFile', 'my_tracking_data.mat', 'arena_cm', 50, 'peripheral_width_cm', 5);
```

### Configurable Parameters

| Parameter             | Description                                                                 | Default Value     |
| --------------------- | --------------------------------------------------------------------------- | ----------------- |
| `dataFile`            | The name of the `.mat` file containing the tracking data.                   | `'Out_Stim2.mat'` |
| `arena_cm`            | The real-world dimension of the arena in centimeters.                       | `70`              |
| `peripheral_width_cm` | The width of the peripheral border of the arena in centimeters.             | `10`              |
| `movement_threshold`  | The speed threshold in cm/s to define movement.                             | `2.0`             |
| `smoothing_window`    | The size of the moving average window for speed smoothing.                  | `5`               |
| `histogram_bins`      | The number of bins for the histograms in the generated plots.               | `40`              |

## Output

The script generates the following outputs:

1.  **Console Summary:** A summary of the analysis results is printed to the MATLAB command window.
2.  **Figures:** Two PNG files containing the visualizations are saved to the current directory:
    -   `TrAQ_Analysis_Periph<width>cm_<timestamp>_speed.png`
    -   `TrAQ_Analysis_Periph<width>cm_<timestamp>_zones.png`
3.  **Data Files:** Two CSV files are saved to the current directory:
    -   `TrAQ_Analysis_Periph<width>cm_<timestamp>_detailed_data.csv`
    -   `TrAQ_Analysis_Periph<width>cm_<timestamp>_summary.csv`
