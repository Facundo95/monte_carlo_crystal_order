#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:27:03 2026

@author: facundo
"""
import matplotlib.pyplot as plt
import numpy as np
import argparse
import re


def read_header(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()

    if not first_line:
        raise ValueError("Empty file. Expected a header line.")

    if first_line.startswith("#"):
        first_line = first_line[1:].strip()

    headers = first_line.split()
    if not headers:
        raise ValueError("Could not parse header names from first line.")

    return headers


def get_data(file_path):
    headers = read_header(file_path)
    data = np.genfromtxt(
        file_path,
        names=headers,
        dtype=float,
        skip_header=1,
        invalid_raise=True,
    )

    if data.size == 0:
        raise ValueError("No numeric rows found in input file.")

    # genfromtxt returns a scalar record when there is only one data row.
    if data.shape == ():
        data = np.array([data], dtype=data.dtype)

    return data


def find_column(columns, candidates):
    lowered = {c.lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def mean_by_temperature_and_field(data, temp_col, field_col):
    temp = data[temp_col]
    field = data[field_col]

    keys = np.column_stack((temp, field))
    unique_keys, inverse = np.unique(keys, axis=0, return_inverse=True)
    counts = np.bincount(inverse)

    grouped = {
        temp_col: unique_keys[:, 0],
        field_col: unique_keys[:, 1],
    }

    for col in data.dtype.names:
        if col in (temp_col, field_col):
            continue
        values = data[col]
        sums = np.bincount(inverse, weights=values)
        grouped[col] = sums / counts

    # Sort by field, then by temperature for cleaner lines.
    order = np.lexsort((grouped[temp_col], grouped[field_col]))
    for col in grouped:
        grouped[col] = grouped[col][order]

    return grouped


def maybe_filter_by_field(grouped, field_col, target_field):
    if target_field is None:
        return grouped

    mask = np.isclose(grouped[field_col], target_field)
    if not np.any(mask):
        raise ValueError(f"No rows found for field h={target_field}.")

    filtered = {}
    for col, values in grouped.items():
        filtered[col] = values[mask]
    return filtered


def plot_scalar_vs_temp(grouped, temp_col, field_col, value_col, title, ylabel):
    fig, ax = plt.subplots(figsize=(8, 5))

    fields = np.unique(grouped[field_col])
    for h in fields:
        mask = np.isclose(grouped[field_col], h)
        ax.plot(
            grouped[temp_col][mask],
            grouped[value_col][mask],
            marker="o",
            linewidth=1.5,
            label=f"{value_col} (h={h:g})",
        )

    ax.set_title(title)
    ax.set_xlabel(temp_col)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()


def extract_lro_columns(columns):
    # Matches x_*, y_*, z_* where * identifies the element/site label.
    by_element = {}
    for col in columns:
        match = re.match(r"^([xyz])_(.+)$", col, flags=re.IGNORECASE)
        if not match:
            continue

        axis = match.group(1).lower()
        element = match.group(2)
        by_element.setdefault(element, {})[axis] = col

    return by_element


def plot_lro(grouped, temp_col, field_col, lro_map):
    elements = sorted(lro_map.keys())
    if not elements:
        raise ValueError("No LRO columns found. Expected headers like x_*, y_*, z_*.")

    fig, axes = plt.subplots(
        len(elements),
        1,
        figsize=(10, 3.2 * len(elements)),
        sharex=True,
    )

    if len(elements) == 1:
        axes = [axes]

    fields = np.unique(grouped[field_col])

    for ax, element in zip(axes, elements):
        element_cols = lro_map[element]
        for axis in ("x", "y", "z"):
            if axis not in element_cols:
                continue
            col_name = element_cols[axis]
            for h in fields:
                mask = np.isclose(grouped[field_col], h)
                ax.plot(
                    grouped[temp_col][mask],
                    grouped[col_name][mask],
                    marker="o",
                    linewidth=1.3,
                    label=f"{col_name} (h={h:g})",
                )

        ax.set_title(f"Element/site: {element}")
        ax.set_ylabel("LRO")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncol=2)

    axes[-1].set_xlabel(temp_col)
    fig.suptitle("LRO parameters (mean by temperature and field)", y=1.01)
    fig.tight_layout()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot .out data averaged by temperature and field."
    )
    parser.add_argument(
        "-f",
        "--file",
        type=str,
        required=True,
        help="Path to .out file.",
    )
    parser.add_argument(
        "-p",
        "--plot",
        type=str,
        choices=["lro", "magn", "energy", "all"],
        default="all",
        help="Type of plot to generate.",
    )
    parser.add_argument(
        "--field",
        type=float,
        default=None,
        help="Optional field value h to filter before plotting.",
    )
    args = parser.parse_args()

    data = get_data(args.file)
    columns = list(data.dtype.names)

    field_col = find_column(columns, ["h", "field", "campo"])
    temp_col = find_column(columns, ["temperature", "temp", "t"])
    magn_col = find_column(columns, ["magnetization", "magn", "m"])
    energy_col = find_column(columns, ["etotal", "energy", "e"])

    if field_col is None or temp_col is None:
        raise ValueError(
            "Could not identify temperature/field columns from header. "
            "Expected names like 'temperature' and 'h'."
        )

    grouped = mean_by_temperature_and_field(data, temp_col, field_col)
    grouped = maybe_filter_by_field(grouped, field_col, args.field)

    lro_map = extract_lro_columns(columns)

    if args.plot in ("lro", "all"):
        plot_lro(grouped, temp_col, field_col, lro_map)

    if args.plot in ("magn", "all"):
        if magn_col is None:
            raise ValueError("Magnetization column not found in header.")
        plot_scalar_vs_temp(
            grouped,
            temp_col,
            field_col,
            magn_col,
            "Magnetization (mean by temperature and field)",
            magn_col,
        )

    if args.plot in ("energy", "all"):
        if energy_col is None:
            raise ValueError("Energy column not found in header.")
        plot_scalar_vs_temp(
            grouped,
            temp_col,
            field_col,
            energy_col,
            "Total energy (mean by temperature and field)",
            energy_col,
        )

    plt.show()