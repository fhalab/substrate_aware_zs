from __future__ import annotations

import os

import pandas as pd

from cairosvg import svg2png

import seaborn as sns

import matplotlib.pyplot as plt

import bokeh
from bokeh.io import export_svg
from bokeh.models import NumeralTickFormatter
from bokeh.themes.theme import Theme

import holoviews as hv

bokeh.io.output_notebook()
hv.extension("bokeh", "matplotlib")


FZL_PALETTES = {
    "blue": "#4bacc6",
    "orange": "#f79646ff",
    "light_orange": "#ffbb78",
    "red": "#ff8888",
    "maroon": "#7A303F",
    "green": "#9bbb59",
    "yellow": "#f9be00",
    "purple": "#8064a2",
    "brown": "#ae682f",
    "dark_brown": "#6e4a2eff",
    "gray": "#666666",
    "light_gray": "#D3D3D3",
    "light_blue":"#849895",
    "light_green": "#9DAE88",
    "light_yellow": "#F1D384",
    "light_brown": "#C7B784",
    "black": "#000000",
}
