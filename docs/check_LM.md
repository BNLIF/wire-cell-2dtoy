# Light Matching Functions Documentation

## Overview
This document explains three related functions in the ToyFiducial class that perform light matching checks: `check_LM()`, `check_LM_cuts()`, and `check_LM_bdt()`. These functions evaluate whether a cluster matches with flash data using different approaches.

## Common Features
All three functions share some common elements:

1. **Input Parameters**:
   - `bundle`: A FlashTPCBundle object containing cluster and flash information
   - `cluster_length`: A reference parameter to store the calculated cluster length

2. **Initial Checks**:
   - Calculate cluster length using geometry parameters
   - Check if predicted PE > 25 and cluster length > 10 cm
   - Return code 1 for low energy events that fail these criteria

3. **Return Values**:
   - 0: Pass (good match)
   - 1: Low energy event
   - 2: Light mismatch

## check_LM()
The original light matching function implementing basic geometric and light matching criteria.

### Key Features:
1. **Position Flags**:
   - `flag_anode`: Whether cluster is close to PMTs
   - `flag_boundary`: Whether cluster is at x boundary

2. **Flash Type Check**:
   - Special handling for flash type 2 (beam discriminator flash)

3. **Matching Criteria**:
   For non-boundary cases:
   ```
   log10(total_pred_pe/total_flash_pe) > -0.55 AND
   ks_dis < 0.25 AND
   ks_dis - 0.15/1.4 * log10(total_pred_pe/total_flash_pe) < 0.32
   ```

   For boundary cases:
   ```
   ks_dis < 0.6 AND
   (log10(total_pred_pe/total_flash_pe) > -1.4 AND flag_anode OR 
    log10(total_pred_pe/total_flash_pe) > -0.55) AND
   log10(total_flash_pe) + 1.8*log10(total_pred_pe/total_flash_pe) > 1.25 AND
   max_meas_pe/total_flash_pe + 0.17/0.5*log10(total_pred_pe/total_flash_pe) > -0.05
   ```

## check_LM_cuts()
A simplified version with modified acceptance criteria.

### Key Features:
1. **Boundary Cases**:
   - At anode:
     ```
     log(total_pred_pe/total_meas_pe) > -1.8 AND
     ks_dis < 0.8
     ```
   - At cathode:
     ```
     log(total_pred_pe/total_meas_pe) > -1.8 AND
     ks_dis < 0.45
     ```

2. **Non-boundary Cases**:
   ```
   log(total_pred_pe/total_meas_pe) > -1.4 AND
   ks_dis < 0.25
   ```

## check_LM_bdt()
The most sophisticated version using Boosted Decision Trees (BDT) for classification.

### Key Features:
1. **BDT Input Variables**:
   - Total predicted PE
   - Total measured PE
   - Maximum measured PE
   - KS distance
   - Chi-square
   - Number of degrees of freedom
   - Cluster length
   - Flags for anode and boundary position

2. **Classification Process**:
   - Creates LMBDT object with input parameters
   - Uses pre-trained BDT model parameters from input files
   - Applies score threshold based on background efficiency (0.005)

### BDT Implementation Details:
- Loads BDT parameters from external files:
  - `lmHistParams.txt`: Histogram parameters
  - `lmSigEff.txt`: Signal efficiency curve
  - `lmBgdEff.txt`: Background efficiency curve
- Creates TH1D histograms for signal and background efficiency
- Uses BDT score with 0.5% background acceptance rate as threshold

## Performance Considerations
- `check_LM()`: Most flexible but complex criteria
- `check_LM_cuts()`: Simplified, more stringent cuts
- `check_LM_bdt()`: Most sophisticated, potentially best discrimination but depends on training quality

## Usage Notes
1. All functions require proper initialization of the FlashTPCBundle with:
   - Valid flash data
   - Cluster information
   - Predicted PE values

2. BDT version requires additional input files and proper setup of the LMBDT class

3. Consider computational cost:
   - `check_LM()`: Medium
   - `check_LM_cuts()`: Low
   - `check_LM_bdt()`: Highest