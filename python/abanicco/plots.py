"""Plotting for ABANICCO: matplotlib figures matching the MATLAB output,
plus interactive Plotly polar description (hover = coordinates + shade name),
mirroring the planned FIJI plugin interactions."""

import numpy as np
import matplotlib.pyplot as plt

from .core import load_regions, wrap_to_360


def show_segmentations(an_im, regions, cols=4):
    """Original image + one tile per category segmentation."""
    live = [r for r in regions]
    n = len(live) + 1
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.2 * rows))
    axes = np.atleast_1d(axes).ravel()
    axes[0].imshow(an_im)
    axes[0].set_title("New Image")
    for ax, reg in zip(axes[1:], live):
        ax.imshow(reg.segmented if reg.segmented is not None else np.zeros_like(an_im))
        ax.set_title(f"{reg.name}-Region ({reg.mask_percentage:.1f}%)")
    for ax in axes:
        ax.axis("off")
    fig.tight_layout()
    return fig


def _draw_region_ring(ax, regions):
    """Colored fan ring (rho 100-145) showing each category's angular span."""
    for reg in regions:
        t1, t2 = sorted((reg.theta1, reg.theta2))
        theta = np.linspace(t1, t2, 60)
        ax.fill_between(theta, reg.rho1, reg.rho2, color=reg.code,
                        linewidth=0, zorder=1)


def _draw_guides(ax, regions, lines, lines2):
    for i in range(lines.shape[0]):
        ax.plot(lines[i], lines2[i] * 98, "k-", lw=0.5, zorder=2)
    th1 = np.linspace(0, 2 * np.pi, 50)
    ax.plot(th1, np.full_like(th1, 3), "k-", lw=1, zorder=2)
    brown = regions[-2]
    th2 = np.linspace(brown.theta1, brown.theta2, 50)
    ax.plot(th2, np.full_like(th2, 35), "k-", lw=1, zorder=2)
    ax.scatter([r.base for r in regions], [r.rho1 for r in regions],
               marker="*", c="k", zorder=4)


def polar_description(regions, rows, lines=None, lines2=None, bubble=True):
    """The ABANICCO polar description: fan of categories + detected shades.

    bubble=True sizes markers by percentage of the shade in the image.
    """
    if lines is None or lines2 is None:
        _, lines, lines2, _, _ = load_regions()
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="polar")
    _draw_region_ring(ax, regions)
    _draw_guides(ax, regions, lines, lines2)
    if rows:
        th = [r["theta"] for r in rows]
        rho = [r["rho"] for r in rows]
        col = [r["rgb"] for r in rows]
        if bubble:
            per = np.array([r["percent"] for r in rows])
            size = 40 + 900 * per / max(per.max(), 1e-9)
        else:
            size = 100
        ax.scatter(th, rho, s=size, c=col, edgecolors="k",
                   linewidths=0.4, zorder=3)
    ax.set_rlim(0, 150)
    # radius measurement along the downward axis
    ax.set_rlabel_position(270)
    ax.set_yticks([30, 60, 90, 120, 150])
    ax.set_yticklabels(["30", "60", "90", "120", "150"],
                       fontsize=8, color="#444")
    ax.text(np.radians(270), 175, "radius \u03c1 = \u221a(a\u00b2+b\u00b2)",
            ha="center", va="top", fontsize=9, color="#444")
    fig.tight_layout()
    return fig


def shade_charts(an_im, regions, cols=4):
    """Per-category bubble chart (rho vs angle) next to its segmentation."""
    live = [r for r in regions if r.shades]
    n = 2 * len(live) + 1
    rows_ = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows_, cols, figsize=(4.2 * cols, 3.4 * rows_))
    axes = np.atleast_1d(axes).ravel()
    axes[0].imshow(an_im)
    axes[0].set_title("Image")
    axes[0].axis("off")
    i = 1
    for reg in live:
        ax = axes[i]; i += 1
        per = np.array([s["percent"] for s in reg.shades])
        ax.scatter([s["rho"] for s in reg.shades],
                   [wrap_to_360(np.degrees(s["theta"])) for s in reg.shades],
                   s=60 + 800 * per / max(per.max(), 1e-9),
                   c=[s["rgb"] for s in reg.shades], edgecolors="k")
        ax.set_title(f"The {reg.name}s")
        ax.set_xlabel("Radius"); ax.set_ylabel("Angle (º)")
        ax2 = axes[i]; i += 1
        ax2.imshow(reg.segmented)
        ax2.set_title(f"The Segmented {reg.name}s")
        ax2.axis("off")
    for ax in axes[i:]:
        ax.axis("off")
    fig.tight_layout()
    return fig


def interactive_polar(regions, rows, lines=None, lines2=None, bubble=True):
    """Plotly interactive polar description.

    Hover a shade -> radius, angle, L value and mixing name.
    Hover/click the outer fan -> category name and % of image.
    bubble=True sizes shade markers by their % of the image;
    bubble=False draws fixed-size dots (scatter mode, easier to read
    when many shades cluster densely).
    """
    import plotly.graph_objects as go
    if lines is None or lines2 is None:
        _, lines, lines2, _, _ = load_regions()

    fig = go.Figure()

    # category fan ring
    for reg in regions:
        t1, t2 = sorted((reg.theta1, reg.theta2))
        theta = np.degrees(np.linspace(t1, t2, 40))
        r_lo, r_hi = reg.rho1, reg.rho2
        fig.add_trace(go.Scatterpolar(
            theta=np.concatenate([theta, theta[::-1]]),
            r=np.concatenate([np.full_like(theta, r_lo),
                              np.full_like(theta, r_hi)]),
            fill="toself", mode="lines",
            line=dict(color=reg.code, width=0.5),
            fillcolor=reg.code, name=reg.name,
            hoverinfo="text",
            text=f"{reg.name}: {reg.mask_percentage:.2f}% of image",
        ))

    # boundary guide lines
    for i in range(lines.shape[0]):
        fig.add_trace(go.Scatterpolar(
            theta=np.degrees(lines[i]), r=lines2[i] * 98,
            mode="lines", line=dict(color="black", width=0.6),
            showlegend=False, hoverinfo="skip"))

    # detected shades
    if rows:
        per = np.array([r["percent"] for r in rows])
        if bubble:
            size = 12 + 38 * per / max(per.max(), 1e-9)
        else:
            size = np.full(len(rows), 9.0)
        fig.add_trace(go.Scatterpolar(
            theta=[wrap_to_360(np.degrees(r["theta"])) for r in rows],
            r=[r["rho"] for r in rows],
            mode="markers",
            marker=dict(
                size=size,
                color=["rgb(%d,%d,%d)" % tuple(int(c * 255) for c in r["rgb"])
                       for r in rows],
                line=dict(color="black", width=1)),
            showlegend=False,
            hoverinfo="text",
            text=[(f"<b>{r['region']}</b><br>"
                   f"radius = {r['rho']:.2f}<br>"
                   f"angle = {r.get('angle', 0):.2f}º<br>"
                   f"L = {r['L']:.2f}<br>"
                   f"{r['percent']:.2f}% of image<br>"
                   f"{r.get('shade_name', '')}") for r in rows],
        ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                range=[0, 150], showticklabels=True, angle=270,
                tickangle=0, tickvals=[30, 60, 90, 120, 150],
                tickfont=dict(size=10, color="#444"),
                ticks="outside", showline=True, linecolor="#888",
                gridcolor="#e5e5e5"),
            angularaxis=dict(rotation=0, direction="counterclockwise")),
        width=760, height=700, template="plotly_white",
        title=("ABANICCO polar description — "
               + ("bubble (size = % of image)" if bubble else "scatter (fixed size)")),
        annotations=[dict(
            text="radius \u03c1 = \u221a(a\u00b2+b\u00b2)", showarrow=False,
            xref="paper", yref="paper", x=0.5, y=-0.06,
            font=dict(size=11, color="#444"))],
    )
    return fig


def cluster_polar(regions, rows, labels, lines=None, lines2=None, title=None):
    """Interactive polar description with shades colored by cluster label.

    Same layout as ``interactive_polar`` (category fan ring + guides), but
    every shade marker is painted with its cluster color; noise points
    (label -1, from DBSCAN/OPTICS) are drawn as grey open circles.
    """
    import plotly.graph_objects as go
    import plotly.express as px
    if lines is None or lines2 is None:
        _, lines, lines2, _, _ = load_regions()

    fig = go.Figure()

    # category fan ring (context, low opacity so clusters stand out)
    for reg in regions:
        t1, t2 = sorted((reg.theta1, reg.theta2))
        theta = np.degrees(np.linspace(t1, t2, 40))
        fig.add_trace(go.Scatterpolar(
            theta=np.concatenate([theta, theta[::-1]]),
            r=np.concatenate([np.full_like(theta, reg.rho1),
                              np.full_like(theta, reg.rho2)]),
            fill="toself", mode="lines", opacity=0.45,
            line=dict(color=reg.code, width=0.5),
            fillcolor=reg.code, name=reg.name, showlegend=False,
            hoverinfo="text",
            text=f"{reg.name}: {reg.mask_percentage:.2f}% of image"))

    for i in range(lines.shape[0]):
        fig.add_trace(go.Scatterpolar(
            theta=np.degrees(lines[i]), r=lines2[i] * 98,
            mode="lines", line=dict(color="#999", width=0.5),
            showlegend=False, hoverinfo="skip"))

    labels = np.asarray(labels)
    palette = px.colors.qualitative.Dark24
    for lb in sorted(set(labels.tolist())):
        idx = np.where(labels == lb)[0]
        members = [rows[i] for i in idx]
        noise = lb == -1
        color = "#aaaaaa" if noise else palette[lb % len(palette)]
        fig.add_trace(go.Scatterpolar(
            theta=[wrap_to_360(np.degrees(r["theta"])) for r in members],
            r=[r["rho"] for r in members],
            mode="markers",
            name="noise" if noise else f"cluster {lb}",
            marker=dict(
                size=11, color="rgba(0,0,0,0)" if noise else color,
                line=dict(color=color, width=2)),
            hoverinfo="text",
            text=[(f"<b>{'noise' if noise else f'cluster {lb}'}</b><br>"
                   f"{r['region']}<br>"
                   f"radius = {r['rho']:.2f} | angle = "
                   f"{wrap_to_360(np.degrees(r['theta'])):.2f}\u00ba | "
                   f"L = {r['L']:.2f}<br>"
                   f"{r['percent']:.2f}% of image<br>"
                   f"{r.get('shade_name', '')}") for r in members]))

    n_found = len(set(labels.tolist()) - {-1})
    fig.update_layout(
        polar=dict(
            radialaxis=dict(range=[0, 150], showticklabels=True, angle=270,
                            tickvals=[30, 60, 90, 120, 150],
                            tickfont=dict(size=10, color="#444"),
                            ticks="outside", showline=True, linecolor="#888",
                            gridcolor="#e5e5e5"),
            angularaxis=dict(rotation=0, direction="counterclockwise")),
        width=760, height=700, template="plotly_white",
        title=title or f"Automatic clustering — {n_found} clusters found",
        legend=dict(orientation="h", yanchor="bottom", y=-0.14),
    )
    return fig
