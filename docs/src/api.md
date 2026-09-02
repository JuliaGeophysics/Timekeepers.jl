# API Reference

Everything Timekeepers exports, grouped by what it is for.

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
metronix_site_runs
write_metronix_site
write_metronix_site_masked
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

## Index

```@index
```
