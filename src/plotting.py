import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import string
#from matplotlib_venn import venn2, venn2_circles
import matplotlib.patches as mpatches  
from matplotlib import transforms as mtransforms
from matplotlib.colors import to_rgba

import matplotlib as mpl
from matplotlib import font_manager as fm


mpl.rcParams.update({
    "text.usetex": False,              # use MathText
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


class Plot:
    def __init__(self, n_rows=1, n_cols=1, figsize=(10, 5), 
                 title_fontdict=None, label_fontdict=None, tick_fontdict=None,
                 shared_title=None, shared_title_fontdict=None, legend_fontdict=None,
                 share_x=False, share_y=False):
        """
        Initialize the Plot class with optional shared title and shared axis settings.
        """
        self.custom_legends = {}
        self.figure_legends = []
        self.fig = plt.figure(figsize=figsize)
        self.gs = gridspec.GridSpec(n_rows, n_cols, figure=self.fig)
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.shared_title = shared_title
        self.shared_title_fontdict = shared_title_fontdict if shared_title_fontdict else {'fontsize': 18, 'fontweight': 'bold'}
        self.axes = {}  # Store axes
        self.plot_data = {}  # Store plot data
        self.textboxes = {}
        # Default font settings
        self.title_fontdict = title_fontdict if title_fontdict else {'fontsize': 14, 'fontweight': 'bold'}
        self.label_fontdict = label_fontdict if label_fontdict else {'fontsize': 12}
        self.tick_fontdict = tick_fontdict if tick_fontdict else {'fontsize': 10}
        self.legend_fontdict = legend_fontdict if legend_fontdict else {'size': 10}
        

    def _get_axes(self, row, col, colspan=1):
        """Retrieve or create axes, merging subplots if needed."""
        key = (row, col, colspan)
        if key not in self.axes:
            self.axes[key] = self.fig.add_subplot(self.gs[row, col:col+colspan])
        return self.axes[key]

    def _facecolor_with_alpha(self, colour, alpha, n=None):
        """
        Return color(s) with alpha applied to face only.
        - colour: str/tuple or list-like of colors
        - alpha: None, float in [0,1], or list-like (per-bar/per-series)
        - n: optional length hint (e.g., number of bars)
        """
        # No alpha provided -> keep colour as-is
        if alpha is None:
            return colour

        # Normalize list-ness
        is_color_list = isinstance(colour, (list, tuple, np.ndarray))
        is_alpha_list = isinstance(alpha, (list, tuple, np.ndarray))

        # If scalar color and list alpha, expand color
        if not is_color_list and is_alpha_list:
            if n is None:
                # best effort
                return [to_rgba(colour, a) for a in alpha]
            return [to_rgba(colour, alpha[i] if i < len(alpha) else alpha[-1]) for i in range(n)]

        # If list color and scalar alpha
        if is_color_list and not is_alpha_list:
            return [to_rgba(c, alpha) for c in colour]

        # If list color and list alpha
        if is_color_list and is_alpha_list:
            L = min(len(colour), len(alpha))
            return [to_rgba(colour[i], alpha[i]) for i in range(L)] + \
                   [to_rgba(colour[i], alpha[-1]) for i in range(L, len(colour))]

        # Both scalar
        return to_rgba(colour, alpha)

    def add_hist(self, row, col, data, bins=10, title=None, overlay=False, colour=None, alpha=None, label=None, x_label=None, y_label=None, ymax=None, colspan=1, rotation=0, labelpos=(0,1.1), bar_edgecolor="Black", bar_edgewidth=0.5, start_x_at_zero=True):
        """
        Add histogram to a specific subplot at (row, col).
        """
        ax = self._get_axes(row, col, colspan)
        if overlay:
            self.plot_data.setdefault((row, col, colspan), []).append(("hist", data, bins, colour, alpha, label, x_label, y_label, title, ymax, rotation, labelpos, bar_edgecolor, bar_edgewidth, start_x_at_zero))
        else:
            self.plot_data[(row, col, colspan)] = [("hist", data, bins, colour, alpha, label, x_label, y_label, title, ymax, rotation, labelpos, bar_edgecolor, bar_edgewidth, start_x_at_zero)]


    def add_bar(self, row, col, categories, values, bar_width=0.35, title=None, colour=None, alpha=None, bar_labels=None, x_label=None, y_label=None, side_by_side=False, colspan=1, ymax=None, rotation=0, labelpos=(0,1.1), bar_edgecolor="Black", bar_edgewidth=1.5):
        """Add bar chart to the specified subplot with an option for side-by-side bars."""
        ax = self._get_axes(row, col, colspan)

        if side_by_side:
            self.plot_data[(row, col, colspan)] = [("side-by-side-bar", categories, values, colour, alpha, x_label, y_label, title, bar_labels, bar_width, bar_edgecolor, bar_edgewidth, ymax, rotation, labelpos)]
        else:
            self.plot_data[(row, col, colspan)] = [("bar", categories, values, colour, alpha, bar_edgecolor, bar_edgewidth, x_label, y_label, title, bar_labels, ymax, rotation, labelpos)]

    def add_horizontal_bar(self, row, col, categories, values, bar_height=0.35, title=None, colour=None, alpha=None, bar_labels=None, x_label=None, y_label=None, colspan=1, xmax=None, rotation=0, labelpos=(0, 1.1), bar_edgecolor="Black", bar_edgewidth=1.5, side_by_side=False):
        """
        Add a horizontal bar chart to the specified subplot.
        If side_by_side=True:
          - values should be a list of series (one per group), each of length len(categories)
          - colour/alpha/bar_edgecolor/bar_edgewidth/bar_labels can be scalars or lists per series
        """
        ax = self._get_axes(row, col, colspan)
        if side_by_side:
            # Store as a separate dataset type like vertical side-by-side bars
            self.plot_data[(row, col, colspan)] = [(
                "horizontal-side-by-side-bar",
                categories, values, colour, alpha,
                x_label, y_label, title, bar_labels,
                bar_height, bar_edgecolor, bar_edgewidth,
                xmax, rotation, labelpos
            )]
        else:
            self.plot_data[(row, col, colspan)] = [(
                "horizontal-bar",
                categories, values, bar_height, colour, alpha,
                bar_edgecolor, bar_edgewidth, x_label, y_label,
                title, bar_labels, xmax, rotation, labelpos
            )]

    def add_pie(self, row, col, values, labels=None, colours=None, alpha=None, explode=None,
                title=None, autopct="%1.1f%%", startangle=90, legend=False, legend_title=None,
                legend_loc='best', colspan=1,radius=0.5,
                small_slice_threshold=0.04,  # fraction of total below which labels moved outside
                outside_distance=1.08,       # radial distance for outside labels
                outside_alignment='right', pct_fontsize=12):  # alignment for outside labels
        """Add a pie chart to the specified subplot.
        Args:
            values: sequence of numeric values.
            labels: sequence of labels (optional).
            colours: list/tuple or single colour.
            alpha: scalar or list of alphas matching values.
            explode: list of explode fractions per slice.
            autopct: format string or None for slice value labels.
            startangle: rotation of the start angle.
            legend: if True, add legend instead of srelying on labels on wedges.
        """
        ax = self._get_axes(row, col, colspan)
        # Store pie dataset; apply alpha later during render similar to other plots
        self.plot_data[(row, col, colspan)] = [(
            "pie", values, labels, colours, alpha, explode,
            title, autopct, startangle, legend, legend_title, legend_loc,
            small_slice_threshold, outside_distance, outside_alignment, pct_fontsize, radius
        )]

    def add_stacked_bar(self, row, col, categories, values, title=None, colours=None, alpha=None, bar_labels=None,
                    x_label=None, y_label=None, colspan=1, ymax=None, rotation=0, labelpos=(0, 1.1), bar_edgecolor="Black", bar_edgewidth=1.5):
        """
                Add stacked bar chart to the specified subplot.
                values may be provided in either orientation:
                    - layers x categories (preferred): values[layer][category]
                    - categories x layers: values[category][layer] (will be transposed)
                colours and alpha can be scalars or lists with length == n_layers.
        """
        ax = self._get_axes(row, col, colspan)
        self.plot_data[(row, col, colspan)] = [("stacked-bar", categories, values, colours, alpha,
                                                x_label, y_label, title, bar_labels, bar_edgecolor, bar_edgewidth, ymax, rotation, labelpos)]


    def add_venn2(self, row, col, subsets, set_labels=("A", "B"), set_colors=("blue", "#FF7701"),
                  alpha=0.8, title=None, normalize_to=1.0, circle_color="black", circle_alpha=1.0,
                  circle_lw=1.0, set_label_fontsize=None, subset_fontsize=None,
                  set_label_fontweight=None, subset_fontweight=None, colspan=1):
        """
        Add a 2-set Venn diagram.
        :param subsets: tuple (A_only, B_only, A_and_B) or set-like inputs accepted by matplotlib_venn.
        :param set_label_fontsize: font size for the set labels ('A','B').
        :param subset_fontsize: font size for the region labels ('10','01','11').
        """
        ax = self._get_axes(row, col, colspan)
        self.plot_data[(row, col, colspan)] = [(
            "venn2",
            subsets, set_labels, set_colors, alpha, title, normalize_to,
            circle_color, circle_alpha, circle_lw,
            set_label_fontsize, subset_fontsize,
            set_label_fontweight, subset_fontweight
        )]

    def add_textbox(self, row, col, text, xy=(0.5, 0.5), color="black", fontstyle="normal", fontsize=12, ha="center", va="center", colspan=1, fontweight='normal', bbox=None, coord='axes'):
        """
        Add a custom text box to a subplot at (row, col).
        :param xy: (x, y) tuple.
        :param coord: 'axes' (default), 'data', or 'axes,data' (x in axes coords, y in data coords).
        """
        key = (row, col, colspan)
        self.textboxes.setdefault(key, []).append({
            "text": text,
            "fontweight": fontweight,
            "xy": xy,
            "color": color,
            "fontsize": fontsize,
            "fontstyle": fontstyle,
            "ha": ha,
            "va": va,
            "bbox": bbox,
            "coord": coord,  # NEW
        })
   
    def add_custom_legend(
        self,
        row, col,
        entries,
        title=None,
        anchor=(0.95, 0.95),   # axes-fraction coords
        loc='upper right',
        colspan=1,
        frame=True,
        frame_facecolor='white',
        frame_alpha=0.85,
        frame_edgecolor='black',
        ncol=1,
        label_fontsize=None,
        label_fontweight=None,
        label_color='black',
        title_fontsize=None,
        title_fontweight='bold',
        block_edgecolor='black',
        block_linewidth=0.0,
        borderpad=0.4,
        columnspacing=0.8,
        handletextpad=0.6,
    ):
        """
        Add a custom legend composed of colored blocks (with alpha) and labels.
        entries: list of dicts like:
            {
              'label': 'Reaxys',
              'color': '#FF7701',
              'alpha': 0.5
            }
        Position anywhere using 'anchor' (axes fraction) and 'loc'.
        """
        key = (row, col, colspan)
        self.custom_legends.setdefault(key, []).append({
            'entries': entries,
            'title': title,
            'anchor': anchor,
            'loc': loc,
            'frame': frame,
            'frame_facecolor': frame_facecolor,
            'frame_alpha': frame_alpha,
            'frame_edgecolor': frame_edgecolor,
            'ncol': ncol,
            'label_fontsize': label_fontsize,
            'label_fontweight': label_fontweight,
            'label_color': label_color,
            'title_fontsize': title_fontsize,
            'title_fontweight': title_fontweight,
            'block_edgecolor': block_edgecolor,
            'block_linewidth': block_linewidth,
            'borderpad': borderpad,
            'columnspacing': columnspacing,
            'handletextpad': handletextpad,
        })
    
    def add_figure_legend(
        self,
        entries,
        title=None,
        anchor=(0.5, 0.02),     # figure-fraction coordinates (x,y)
        loc='lower center',     # alignment relative to anchor
        frame=True,
        frame_facecolor='white',
        frame_alpha=0.85,
        frame_edgecolor='black',
        ncol=1,
        label_fontsize=None,
        label_fontweight=None,
        label_color='black',
        title_fontsize=None,
        title_fontweight='bold',
        block_edgecolor='black',
        block_linewidth=0.0,
        borderpad=0.4,
        columnspacing=0.8,
        handletextpad=0.6,
        zorder=10,
    ):
        """
        Add a legend at figure level (does not change subplot layout when moved).
        entries: list of dicts:
          {'label': 'Reaxys', 'color': '#FF7701', 'alpha': 0.4}
        """
        self.figure_legends.append({
            'entries': entries,
            'title': title,
            'anchor': anchor,
            'loc': loc,
            'frame': frame,
            'frame_facecolor': frame_facecolor,
            'frame_alpha': frame_alpha,
            'frame_edgecolor': frame_edgecolor,
            'ncol': ncol,
            'label_fontsize': label_fontsize,
            'label_fontweight': label_fontweight,
            'label_color': label_color,
            'title_fontsize': title_fontsize,
            'title_fontweight': title_fontweight,
            'block_edgecolor': block_edgecolor,
            'block_linewidth': block_linewidth,
            'borderpad': borderpad,
            'columnspacing': columnspacing,
            'handletextpad': handletextpad,
            'zorder': zorder,
        })

    def plot(self, savefig=False, filepath=None, x_tick_colors=None, y_tick_colors=None,
             tight=True, tight_pad=0.3, tight_w_pad=0.3, tight_h_pad=0.3, x_tick_mask=1,
             pie_tight=False):
        """
        Render all plots with optional tick label colors.
        :param x_tick_colors: List of colors for x-tick labels (must match the number of x-tick labels).
        :param y_tick_colors: List of colors for y-tick labels (must match the number of y-tick labels).
        :param x_tick_mask: Integer >=1; if >1 only every nth x tick label is shown (others hidden). Default 1 (show all).
        """
        subplot_labels = list(string.ascii_lowercase)  # ['a', 'b', 'c', ..., 'z']
        label_idx = 0  # Counter for subplot labels
        has_pie = False
        pie_axes = []
        pie_axes_pos = {}
        for (row, col, colspan), datasets in self.plot_data.items():
            ax = self._get_axes(row, col, colspan)
            ax.clear()
            title = None
            use_legend = False
            x_label = None         # NEW: safe defaults
            y_label = None         # NEW: safe defaults
            rotation = 0           # NEW: safe default

             # ...existing code in plot(), replace the textbox drawing block...
            for tb in self.textboxes.get((row, col, colspan), []):
                # choose transform
                coord = tb.get("coord", "axes")
                if coord == "axes":
                    T = ax.transAxes
                elif coord == "data":
                    T = ax.transData
                elif coord in ("axes,data", "data,axes"):
                    T = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
                else:
                    T = ax.transAxes
                ax.text(
                    tb["xy"][0], tb["xy"][1], tb["text"],
                    color=tb["color"], fontsize=tb["fontsize"],
                    ha=tb["ha"], va=tb["va"],
                    transform=T,
                    bbox=tb["bbox"],
                    fontweight=tb.get("fontweight", None),
                    fontstyle=tb.get("fontstyle", None)
                )

            for dataset in datasets:
                if dataset[0] == "hist":
                    _, data, bins, colour, alpha, label, x_label, y_label, dataset_title, ymax, rotation, labelpos, bar_edgecolor, bar_edgewidth, start_x_at_zero = dataset
                    # Apply alpha only to the face (edges remain opaque)
                    face = self._facecolor_with_alpha(colour, alpha) if colour is not None else None
                    ax.hist(
                        data,
                        bins=bins,
                        color=face,                 # RGBA face; no global alpha
                        label=label,
                        edgecolor=bar_edgecolor,
                        linewidth=bar_edgewidth
                    )
                    if label:
                        use_legend = True
                    if ymax:
                        ax.set_ylim(ymax=ymax)
                    if start_x_at_zero:
                        # Ensure the left x-limit starts at zero; keep current right limit
                        current_xlim = ax.get_xlim()
                        ax.set_xlim(left=0, right=current_xlim[1])

                elif dataset[0] == "bar":
                    _, categories, values, colour, alpha, bar_edgecolor, bar_edgewidth, x_label, y_label, dataset_title, bar_labels, ymax, rotation, labelpos = dataset
                    # Apply alpha only to face color
                    face = self._facecolor_with_alpha(colour, alpha, n=len(values) if hasattr(values, "__len__") else None)
                    ax.bar(categories, values, color=face, edgecolor=bar_edgecolor, linewidth=bar_edgewidth, label=bar_labels)
                    if bar_labels:
                        use_legend = True
                    if ymax:
                        ax.set_ylim(ymax=ymax)
                    ax.set_xticks(categories)
                    ax.set_xticklabels(categories)

                elif dataset[0] == "horizontal-bar":
                    _, categories, values, bar_height, colour, alpha, bar_edgecolor, bar_edgewidth, x_label, y_label, dataset_title, bar_labels, xmax, rotation, labelpos = dataset
                    y_positions = np.arange(len(categories))
                    face = self._facecolor_with_alpha(colour, alpha, n=len(values) if hasattr(values, "__len__") else None)
                    ax.barh(y_positions, values, height=bar_height, color=face, edgecolor=bar_edgecolor, linewidth=bar_edgewidth)
                    ax.set_yticks(y_positions)
                    ax.set_yticklabels(categories)
                    if xmax:
                        ax.set_xlim(right=xmax)
                    if bar_labels:
                        for i, value in enumerate(values):
                            ax.text(value, y_positions[i], str(value), va='center', ha='left', fontsize=10)

                elif dataset[0] == "horizontal-side-by-side-bar":
                    _, categories, values, colour, alpha, x_label, y_label, dataset_title, bar_labels, bar_height, bar_edgecolor, bar_edgewidth, xmax, rotation, labelpos = dataset
                    n_groups = len(categories)
                    n_bars = len(values)
                    y_index = np.arange(n_groups)

                    # Normalize list/scalar inputs per series
                    if not isinstance(colour, list):
                        colour = [colour] * n_bars
                    if not isinstance(alpha, list):
                        alpha = [alpha] * n_bars
                    if not isinstance(bar_labels, list):
                        bar_labels = [None] * n_bars
                    if not isinstance(bar_edgecolor, list):
                        bar_edgecolor = [bar_edgecolor] * n_bars
                    if not isinstance(bar_edgewidth, list):
                        bar_edgewidth = [bar_edgewidth] * n_bars

                    for i in range(n_bars):
                        offset = i * bar_height
                        face = self._facecolor_with_alpha(colour[i], alpha[i])
                        ax.barh(
                            y_index + offset,
                            values[i],
                            height=bar_height,
                            color=face,
                            edgecolor=bar_edgecolor[i],
                            linewidth=bar_edgewidth[i],
                            label=bar_labels[i]
                        )

                    # Center ticks across grouped bars
                    ax.set_yticks(y_index + (n_bars - 1) * bar_height / 2)
                    ax.set_yticklabels(categories)
                    if xmax:
                        ax.set_xlim(right=xmax)
                    if any(bar_labels):
                        use_legend = True

                elif dataset[0] == "side-by-side-bar":
                    _, categories, values, colour, alpha, x_label, y_label, dataset_title, bar_labels, bar_width, bar_edgecolor, bar_edgewidth, ymax, rotation, labelpos = dataset
                    n_groups = len(categories)
                    n_bars = len(values)
                    bar_index = np.arange(n_groups)

                    # Normalize list/scalar inputs
                    if not isinstance(colour, list):
                        colour = [colour] * n_bars
                    if not isinstance(alpha, list):
                        alpha = [alpha] * n_bars
                    if not isinstance(bar_labels, list):
                        bar_labels = [None] * n_bars
                    if not isinstance(bar_edgecolor, list):
                        bar_edgecolor = [bar_edgecolor] * n_bars
                    if not isinstance(bar_edgewidth, list):
                        bar_edgewidth = [bar_edgewidth] * n_bars

                    for i, (val, col, alph, lab, ecolor, ewidth) in enumerate(zip(values, colour, alpha, bar_labels, bar_edgecolor, bar_edgewidth)):
                        offset = i * bar_width
                        # Series-wide face color; alpha only on face, not edge
                        face = self._facecolor_with_alpha(col, alph)
                        ax.bar(bar_index + offset, val, bar_width, color=face, label=lab, edgecolor=ecolor, linewidth=ewidth)

                    ax.set_xticks(bar_index + (n_bars - 1) * bar_width / 2)
                    ax.set_xticklabels(categories)
                    if ymax:
                        ax.set_ylim(ymax=ymax)
                    if any(bar_labels):
                        use_legend = True

                elif dataset[0] == "stacked-bar":
                    _, categories, values, colours, alpha, x_label, y_label, dataset_title, bar_labels, bar_edgecolor, bar_edgewidth, ymax, rotation, labelpos = dataset
                    bar_index = np.arange(len(categories))

                    # Normalize values orientation to layers x categories
                    # Accept list/tuple/np.ndarray
                    values_list = [list(v) if isinstance(v, (list, tuple, np.ndarray)) else [v] for v in values]
                    n_categories = len(categories)

                    # Detect categories-major (categories x layers) and transpose
                    if len(values_list) == n_categories and all(isinstance(v, (list, tuple, np.ndarray)) for v in values_list):
                        inner_len = len(values_list[0]) if values_list else 0
                        # If inner length is not equal to number of categories, likely categories x layers
                        if inner_len != n_categories:
                            try:
                                values_list = [list(x) for x in zip(*values_list)]  # transpose
                            except Exception:
                                raise ValueError("Could not transpose stacked-bar values; ensure rectangular list-of-lists.")

                    # Validate that each layer has one value per category
                    if any(len(layer) != n_categories for layer in values_list):
                        raise ValueError(
                            "Each stacked layer must have len == len(categories). "
                            f"Got category count {n_categories} and layer lengths {[len(l) for l in values_list]}"
                        )

                    n_stacks = len(values_list)
                    bottom = np.zeros(n_categories)

                    # Normalize list/scalar inputs to n_stacks
                    if not isinstance(colours, list):
                        colours = [colours] * n_stacks
                    if not isinstance(alpha, list):
                        alpha = [alpha] * n_stacks
                    if not isinstance(bar_labels, list):
                        bar_labels = [None] * n_stacks
                    if not isinstance(bar_edgecolor, list):
                        bar_edgecolor = [bar_edgecolor] * n_stacks
                    if not isinstance(bar_edgewidth, list):
                        bar_edgewidth = [bar_edgewidth] * n_stacks

                    for i in range(n_stacks):
                        face = self._facecolor_with_alpha(colours[i], alpha[i])
                        layer_vals = np.array(values_list[i])
                        ax.bar(
                            bar_index,
                            layer_vals,
                            bottom=bottom,
                            color=face,
                            label=bar_labels[i],
                            edgecolor=bar_edgecolor[i],
                            linewidth=bar_edgewidth[i],
                        )
                        bottom += layer_vals

                    ax.set_xticks(bar_index)
                    ax.set_xticklabels(categories)
                    if any(bar_labels):
                        use_legend = True
                    if ymax is not None:
                        ax.set_ylim(top=ymax)

                elif dataset[0] == "pie":
                    has_pie = True
                    (_tag, values, labels, colours, alpha_p, explode, dataset_title,
                     autopct, startangle, pie_legend, legend_title, legend_loc,
                     small_slice_threshold, outside_distance, outside_alignment, pct_fontsize, radius) = dataset
                    n = len(values)
                    # Normalize colours
                    if colours is None:
                        colours_list = None
                    else:
                        if not isinstance(colours, (list, tuple)):
                            colours_list = [colours] * n
                        else:
                            colours_list = list(colours)
                            if len(colours_list) < n:
                                colours_list += [colours_list[-1]] * (n - len(colours_list))
                    # Apply alpha per slice if provided
                    if alpha_p is not None:
                        if not isinstance(alpha_p, (list, tuple)):
                            alpha_list = [alpha_p] * n
                        else:
                            alpha_list = list(alpha_p)
                            if len(alpha_list) < n:
                                alpha_list += [alpha_list[-1]] * (n - len(alpha_list))
                        if colours_list is not None:
                            colours_list = [to_rgba(c, a) for c, a in zip(colours_list, alpha_list)]
                    wedges, texts, autotexts = ax.pie(
                        values,
                        labels=None if pie_legend else labels,
                        colors=colours_list,
                        explode=explode if explode is not None else None,
                        autopct=autopct,
                        startangle=startangle,
                        radius=radius,
                        textprops={"fontsize": pct_fontsize if pct_fontsize else 10},
                    )
                    # Ensure a circular pie and make visual size respond to `radius`
                    # Use datalim so the axes respect the data limits set by radius
                    ax.set_aspect('equal', adjustable='datalim')
                    ax.set_xlim(-radius, radius)
                    ax.set_ylim(-radius, radius)
                    # Track pie axes and their positions to protect from tight_layout if requested
                    if ax not in pie_axes:
                        pie_axes.append(ax)
                        pie_axes_pos[ax] = ax.get_position().frozen()
                    # Reposition small slice percentage labels outside
                    if autotexts and small_slice_threshold is not None and small_slice_threshold > 0:
                        total = float(sum(values)) if sum(values) != 0 else 1.0
                        # Build combined outside labels for small slices (label + percent)
                        for idx, (w, at) in enumerate(zip(wedges, autotexts)):
                            val = values[idx]
                            fraction = val / total if total else 0
                            if fraction < small_slice_threshold:
                                ang = (w.theta2 + w.theta1) / 2.0
                                rad = np.deg2rad(ang)
                                x = np.cos(rad) * outside_distance
                                y = np.sin(rad) * outside_distance
                                # Retrieve percent text
                                pct_text = at.get_text()
                                label_text = labels[idx] if labels and idx < len(labels) and not pie_legend else ''
                                combined = f"{label_text} {pct_text}".strip()
                                at.set_text(combined)
                                at.set_fontsize(pct_fontsize)
                                at.set_position((x, y))
                                ha = 'right' if outside_alignment == 'right' else ('left' if outside_alignment == 'left' else 'center')
                                if outside_alignment == 'right' and x < 0:
                                    ha = 'left'
                                at.set_ha(ha)
                                at.set_va('center')
                                at.set_bbox(dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none', alpha=0.6))
                                # Hide inside label if we had separate text objects (when not using legend)
                                if texts and not pie_legend and idx < len(texts):
                                    texts[idx].set_visible(False)
                    if pie_legend and labels is not None:
                        leg = ax.legend(wedges, labels, title=legend_title, loc=legend_loc, prop=self.legend_fontdict)
                        if leg.get_title() and self.title_fontdict.get('fontsize'):
                            leg.get_title().set_fontsize(self.title_fontdict['fontsize'])
                    ax.axis('equal')
                    # No need to flip use_legend flag because handled internally

                elif dataset[0] == "venn2":
                    (_tag, subsets, set_labels, set_colors, alpha_v, dataset_title, normalize_to,
                     circle_color, circle_alpha, circle_lw,
                     set_label_fontsize, subset_fontsize,
                     set_label_fontweight, subset_fontweight) = dataset

                    v = venn2(
                        subsets=subsets,
                        set_labels=set_labels,
                        set_colors=set_colors,
                        alpha=alpha_v,
                        normalize_to=normalize_to,
                        ax=ax
                    )
                    venn2_circles(subsets=subsets, ax=ax, color=circle_color, alpha=circle_alpha, linewidth=circle_lw)

                    # Apply font sizes/weights to set labels
                    if v.set_labels is not None:
                        for lbl in v.set_labels:
                            if lbl is not None:
                                if set_label_fontsize is not None:
                                    lbl.set_fontsize(set_label_fontsize)
                                if set_label_fontweight is not None:
                                    lbl.set_fontweight(set_label_fontweight)

                    # Apply font sizes/weights to subset labels ('10','01','11')
                    for rid in ("10", "01", "11"):
                        lab = v.get_label_by_id(rid)
                        if lab is not None:
                            if subset_fontsize is not None:
                                lab.set_fontsize(subset_fontsize)
                            if subset_fontweight is not None:
                                lab.set_fontweight(subset_fontweight)

                    ax.set_axis_off()  # cleaner look
                    dataset_title and ax.set_title(dataset_title, fontdict=self.title_fontdict)

                # capture last non-None title
                if 'dataset_title' in locals() and dataset_title:
                    title = dataset_title

            if use_legend:
                ax.legend(prop=self.legend_fontdict)

            # Render custom legends (preserve any existing legend by adding as artist)
            for spec in self.custom_legends.get((row, col, colspan), []):
                handles = []
                labels = []
                for e in spec['entries']:
                    color = e.get('color', 'gray')
                    alpha = e.get('alpha', 1.0)
                    label = e.get('label', '')
                    patch = mpatches.Patch(
                        facecolor=color,
                        alpha=alpha,
                        edgecolor=spec['block_edgecolor'],
                        linewidth=spec['block_linewidth']
                    )
                    handles.append(patch)
                    labels.append(label)

                leg = ax.legend(
                    handles=handles,
                    labels=labels,
                    loc=spec['loc'],
                    bbox_to_anchor=spec['anchor'],
                    frameon=spec['frame'],
                    ncol=spec['ncol'],
                    title=spec['title'],
                    borderpad=spec['borderpad'],
                    columnspacing=spec['columnspacing'],
                    handletextpad=spec['handletextpad'],
                    prop={'size': spec['label_fontsize']} if spec['label_fontsize'] else None,
                )
                # Style legend frame
                if spec['frame']:
                    leg.get_frame().set_facecolor(spec['frame_facecolor'])
                    leg.get_frame().set_alpha(spec['frame_alpha'])
                    leg.get_frame().set_edgecolor(spec['frame_edgecolor'])
                # Style legend texts
                for txt in leg.get_texts():
                    if spec['label_color']:
                        txt.set_color(spec['label_color'])
                    if spec['label_fontweight']:
                        txt.set_fontweight(spec['label_fontweight'])
                # Style legend title
                if leg.get_title():
                    if spec['title_fontsize']:
                        leg.get_title().set_fontsize(spec['title_fontsize'])
                    if spec['title_fontweight']:
                        leg.get_title().set_fontweight(spec['title_fontweight'])
                ax.add_artist(leg)  # keep any existing legend too

            # Apply font settings to title, labels, and ticks
            ax.set_title(title if title else "", fontdict=self.title_fontdict)
            ax.set_xlabel(x_label if x_label else "", fontdict=self.label_fontdict)
            ax.set_ylabel(y_label if y_label else "", fontdict=self.label_fontdict)
            ax.tick_params(axis='both', labelsize=self.tick_fontdict.get('fontsize', 10))
            ax.tick_params(axis='x', labelrotation=rotation)

            # Apply colors to x-tick labels
            if x_tick_colors:
                x_ticks = ax.get_xticklabels()
                if len(x_ticks) == len(x_tick_colors):
                    for tick, color in zip(x_ticks, x_tick_colors):
                        tick.set_color(color)
                else:
                    raise ValueError("Length of x_tick_colors must match the number of x-tick labels.")
 

            # Apply colors to y-tick labels
            if y_tick_colors:
                y_ticks = ax.get_yticklabels()
                if len(y_ticks) == len(y_tick_colors):
                    for tick, color in zip(y_ticks, y_tick_colors):
                        tick.set_color(color)
                else:
                    raise ValueError("Length of y_tick_colors must match the number of y-tick labels.")


            # Align x-tick labels to the left
            for tick in ax.get_xticklabels():
                tick.set_ha('right')

            # Apply x tick masking (hide labels except every nth)
            if isinstance(x_tick_mask, int) and x_tick_mask > 1:
                current_labels = [lbl.get_text() for lbl in ax.get_xticklabels()]
                # If labels are empty (e.g., numeric ticks auto-generated), populate from tick values
                if not any(current_labels):
                    current_labels = [str(tick) for tick in ax.get_xticks()]
                masked_labels = [lab if (i % x_tick_mask == 0) else '' for i, lab in enumerate(current_labels)]
                ax.set_xticklabels(masked_labels)

        # Apply shared title if specified
        if self.shared_title:
            self.fig.suptitle(self.shared_title, **self.shared_title_fontdict)

        # Render figure-level legends
        # Run layout first so figure legends don’t trigger subplot reflow
        if tight:
            # Use figure.tight_layout to control whitespace for the whole figure
            self.fig.tight_layout(pad=tight_pad, w_pad=tight_w_pad, h_pad=tight_h_pad)
            # If pies should not be tightened, restore their original positions
            if has_pie and not pie_tight:
                for pax in pie_axes:
                    if pax in pie_axes_pos:
                        pax.set_position(pie_axes_pos[pax])

        # Render figure-level legends (positions are in figure coords)
        for spec in self.figure_legends:
            handles = []
            labels = []
            for e in spec['entries']:
                color = e.get('color', 'gray')
                alpha = e.get('alpha', 1.0)
                label = e.get('label', '')
                patch = mpatches.Patch(
                    facecolor=color, alpha=alpha,
                    edgecolor=spec['block_edgecolor'],
                    linewidth=spec['block_linewidth']
                )
                handles.append(patch)
                labels.append(label)

            leg = self.fig.legend(
                handles=handles,
                labels=labels,
                loc=spec['loc'],
                bbox_to_anchor=spec['anchor'],
                bbox_transform=self.fig.transFigure,  # explicit figure coords
                frameon=spec['frame'],
                ncol=spec['ncol'],
                title=spec['title'],
                borderpad=spec['borderpad'],
                columnspacing=spec['columnspacing'],
                handletextpad=spec['handletextpad'],
                prop={'size': spec['label_fontsize']} if spec['label_fontsize'] else None,
            )
            # Style frame
            if spec['frame']:
                leg.get_frame().set_facecolor(spec['frame_facecolor'])
                leg.get_frame().set_alpha(spec['frame_alpha'])
                leg.get_frame().set_edgecolor(spec['frame_edgecolor'])
            # Style label texts
            for txt in leg.get_texts():
                if spec['label_color']:
                    txt.set_color(spec['label_color'])
                if spec['label_fontweight']:
                    txt.set_fontweight(spec['label_fontweight'])
            # Style title
            if leg.get_title():
                if spec['title_fontsize']:
                    leg.get_title().set_fontsize(spec['title_fontsize'])
                if spec['title_fontweight']:
                    leg.get_title().set_fontweight(spec['title_fontweight'])
            leg.set_zorder(spec['zorder'])
            self.fig.add_artist(leg)

        # Apply plt.tight_layout as a final pass, then restore pie axes if needed
        if tight:
            plt.tight_layout()
            if has_pie and not pie_tight:
                for pax in pie_axes:
                    if pax in pie_axes_pos:
                        pax.set_position(pie_axes_pos[pax])
        # Save the figure if you want to
        if savefig == True:
            plt.savefig(filepath, bbox_inches="tight", pad_inches=0.1,dpi=500)
        plt.show()