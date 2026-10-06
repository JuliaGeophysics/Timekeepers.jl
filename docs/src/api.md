# API Reference

This page gives all the exported functions and types, in groups by their
purpose.

```@meta
CurrentModule = Timekeepers
```

## Contents

```@contents
Pages = ["api.md"]
Depth = 2
```

## Data model

```@docs
TimekeeperChannel
TimekeeperRun
components
default_components
sampling_rate
start_time
end_time
duration_seconds
```

## TimeArray interop

```@docs
to_timearray
from_timearray
```

## Format-agnostic I/O

```@docs
read_timekeeper
write_timekeeper
default_data_dir
```

## LEMI-424

```@docs
LEMI424_COLUMNS
read_lemi424
load_lemi424
write_lemi424
```

## GEOMAG

```@docs
read_geomag
load_geomag
write_geomag
```

## Metronix ATS

```@docs
METRONIX_CHANNEL_MAP
read_metronix
load_metronix
write_metronix
is_metronix_site
metronix_site_rates
load_metronix_site
metronix_site_runs
write_metronix_site
write_metronix_site_masked
split_metronix_site
metronix_site_is_split
```

## Masking and cleaning

```@docs
TimekeeperMask
mask_interval!
unmask_interval!
clear_mask!
masked_samples
combine_masks
cleaned_timearray
good_segments
sample_weights
write_cleaned
write_mask
read_mask
```

## Interactive explorer

```@docs
TKApp
run_tkapp
```

## Survey dashboard

```@docs
SurveyRun
SurveySite
Survey
scan_survey
site_rates
site_components
has_magnetic
survey_rates
recording_intervals
recording_seconds
overlap_intervals
overlap_seconds
overlap_matrix
site_distance
site_references
common_window
reference_plan
write_reference_plan
read_reference_plan
TKDash
run_tkdash
```

## Calibration

```@docs
SensorCalibration
read_calibration
find_calibration
sensor_response
```

## Transfer functions

```@docs
TransferFunction
estimate_tf
default_references
apparent_resistivity
impedance_phase
rotate_tf
check_polarity
flip_check!
```

## Transfer function I/O

```@docs
write_edi
read_edi
write_modem
read_modem
plot_tf
export_tf
write_tf_report
```

## Transfer function window

```@docs
TKProc
run_tkproc
```

## Index

```@index
```
