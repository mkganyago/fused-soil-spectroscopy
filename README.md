# Fused soil spectroscopy and feature selection

Data and R scripts accompanying the manuscript *Evaluating the Performance of the Fused VNIR-SWIR-MIR Spectra and Feature Selection for Predicting Soil Properties Using Machine Learning Algorithms* by Nomaxabiso Ngogodo, Mahlatse Kganyago, Goodman Jezile, and Kgabo Humphrey Thamaga.

The study uses archived South African soil samples from the Agricultural Research Council's Soil, Climate and Water laboratory. It evaluates partial least squares regression (PLSR) and radial-basis-function support vector machine (SVM) regression on fused VNIR–SWIR–MIR spectra. The manuscript compares full-spectrum and feature-selected models under several spectral preprocessing treatments.

## Files

| Path | Contents |
| --- | --- |
| `Spectra/combined_data1.csv` | Fused spectral measurements. The `SNM` column identifies the sample and the numeric column names identify spectral bands. |
| `Wet_Chem/` | Laboratory reference measurements in eight CSV files, one each for pH, sand, silt, clay, copper (Cu), manganese (Mn), selenium (Se), and zinc (Zn). Use `SNM` to match observations to the spectral table. |
| `Scipts/PLS.R` | PLS regression workflow. |
| `Scipts/SVM.R` | SVM regression workflow. |
| `Scipts/Feature SELECTION.R` | Competitive adaptive reweighted sampling (CARS) and PLSR model evaluation. |

## Running the scripts

Use R and install the packages listed at the top of the script you intend to run. The package declarations across the three scripts include `pls`, `caret`, `dplyr`, `ggplot2`, `corrplot`, `prospectr`, `readr`, `moments`, `kernlab`, `tidyverse`, `bestNormalize`, `DescTools`, and `gridExtra`.

Review the input paths and expected columns before running. The feature-selection script expects `combined_data1.csv` and property CSVs in its working directory. The PLS and SVM scripts reference inputs absent from this folder, including separate silt, MIR, and sample-availability files in `PLS.R` and a local `lead_full_duplicates.csv` file in `SVM.R`. Those workflows require their expected inputs or code adaptation before they can be run from this deposit.
