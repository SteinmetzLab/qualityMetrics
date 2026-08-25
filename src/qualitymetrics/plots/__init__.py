"""Figures. Each returns a matplotlib Figure and draws nothing to screen."""
from .channels import band_rms, channel_health, depth_power  # noqa: I001
from .channels import filter_state, plot_channel_health
from .diagnostics import (depth_power_wide, noise_cutoff_diagnostic,
                          units_spanning_noise_cutoff)
from .drift import amplitude_cdf_over_depth, drift_map, unit_drift
from .examples import draw_acg, draw_template_waterfall, example_neurons
from .population import (depth_correlation, depth_profiles, firing_rate_image,
                         lfp_band_profiles)
from .rawviews import draw_waveforms_on_probe, raw_traces, raw_with_spikes
from .rawviews import sample_waveforms, spike_snippets, wall_heatmap
from .units import amp_depth_scatter, template_waterfall
from .units import templates_grid, unit_summary_table

__all__ = [
    "amp_depth_scatter", "amplitude_cdf_over_depth", "band_rms",
    "channel_health", "depth_correlation", "depth_power", "depth_power_wide",
    "depth_profiles", "draw_acg", "draw_template_waterfall",
    "draw_waveforms_on_probe", "drift_map", "example_neurons", "filter_state",
    "firing_rate_image", "lfp_band_profiles", "noise_cutoff_diagnostic",
    "plot_channel_health", "raw_traces", "raw_with_spikes", "sample_waveforms",
    "spike_snippets", "template_waterfall", "templates_grid", "unit_drift",
    "unit_summary_table", "units_spanning_noise_cutoff", "wall_heatmap",
]
