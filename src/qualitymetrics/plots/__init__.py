"""Figures. Each returns a matplotlib Figure and draws nothing to screen."""
from .channels import band_rms, channel_health, depth_power  # noqa: I001
from .channels import filter_state, plot_channel_health
from .drift import amplitude_cdf_over_depth, drift_map, unit_drift
from .rawviews import raw_traces, raw_with_spikes, sample_waveforms, wall_heatmap
from .units import amp_depth_scatter, template_waterfall
from .units import templates_grid, unit_summary_table

__all__ = [
    "amp_depth_scatter", "amplitude_cdf_over_depth", "band_rms",
    "channel_health", "depth_power", "drift_map", "filter_state",
    "plot_channel_health", "raw_traces", "raw_with_spikes",
    "sample_waveforms", "template_waterfall", "templates_grid", "unit_drift",
    "unit_summary_table", "wall_heatmap",
]
