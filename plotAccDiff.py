#!/usr/bin/env python3
import argparse
import math
from array import array
import ROOT

# The two plots contain three rapidity slices each.  The values are bin edges.
RAPIDITY_GROUPS = (
    ((0.0, 0.3), (0.3, 0.6), (0.6, 0.9)),
    ((0.9, 1.2), (1.2, 1.8), (1.8, 2.4)),
)

# A distinct color and marker for each y interval.
STYLES = (
    (ROOT.kBlue + 1, 20),
    (ROOT.kRed + 1, 21),
    (ROOT.kGreen + 2, 22),
)


def parse_table(filename):
    """Read the whitespace-separated bin table from Acceptance_reweighted.txt."""
    rows = []
    with open(filename, "r", encoding="utf-8") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            fields = line.split()
            if len(fields) < 10:
                continue

            try:
                values = [float(value) for value in fields[:10]]
            except ValueError:
                # Skip titles, column headers, separators, and summary lines.
                continue

            (pt_low, pt_high, y_low, y_high, denominator, numerator,
             acceptance, weighted_denominator, weighted_numerator,
             weighted_acceptance) = values

            if denominator <= 0:
                raise ValueError(
                    f"Line {line_number}: Denominator must be positive, got {denominator}"
                )
            if numerator < 0 or numerator > denominator:
                raise ValueError(
                    f"Line {line_number}: Numerator={numerator} is inconsistent "
                    f"with Denominator={denominator}"
                )
            if not math.isclose(numerator, round(numerator), abs_tol=1e-6):
                raise ValueError(
                    f"Line {line_number}: Numerator is not an integer count: {numerator}"
                )
            if not math.isclose(denominator, round(denominator), abs_tol=1e-6):
                raise ValueError(
                    f"Line {line_number}: Denominator is not an integer count: {denominator}"
                )
            if pt_high <= pt_low:
                raise ValueError(
                    f"Line {line_number}: invalid pT interval [{pt_low}, {pt_high}]"
                )

            rows.append({
                "pt_low": pt_low,
                "pt_high": pt_high,
                "y_low": y_low,
                "y_high": y_high,
                "denominator": int(round(denominator)),
                "numerator": int(round(numerator)),
                "acceptance_from_file": acceptance,
                "weighted_denominator": weighted_denominator,
                "weighted_numerator": weighted_numerator,
                "weighted_acceptance": weighted_acceptance,
            })

    if not rows:
        raise RuntimeError(f"No bin-by-bin rows were found in '{filename}'.")
    return rows


def same_edge(a, b):
    return math.isclose(a, b, rel_tol=0.0, abs_tol=1e-6)


def select_rapidity_slice(rows, y_low, y_high):
    """Return pT-sorted rows for one rapidity interval."""
    selected = [
        row for row in rows
        if same_edge(row["y_low"], y_low) and same_edge(row["y_high"], y_high)
    ]
    selected.sort(key=lambda row: row["pt_low"])
    if not selected:
        raise RuntimeError(f"No rows found for y interval [{y_low}, {y_high}].")
    return selected


def build_tefficiency(rows, name):
    """Build a TEfficiency object for one y slice, using unweighted counts."""
    rows = sorted(rows, key=lambda row: row["pt_low"])
    edges = [rows[0]["pt_low"]] + [row["pt_high"] for row in rows]

    # Ensure pT bins form one contiguous variable-width histogram.
    for index in range(1, len(rows)):
        if not same_edge(rows[index - 1]["pt_high"], rows[index]["pt_low"]):
            raise RuntimeError(
                f"Non-contiguous pT bins in y slice [{rows[0]['y_low']}, "
                f"{rows[0]['y_high']}]."
            )

    n_bins = len(rows)
    edge_array = array("d", edges)
    h_total = ROOT.TH1D(f"h_total_{name}", "", n_bins, edge_array)
    h_passed = ROOT.TH1D(f"h_passed_{name}", "", n_bins, edge_array)
    h_total.SetDirectory(0)
    h_passed.SetDirectory(0)

    for bin_index, row in enumerate(rows, start=1):
        h_total.SetBinContent(bin_index, row["denominator"])
        h_passed.SetBinContent(bin_index, row["numerator"])

    if not ROOT.TEfficiency.CheckConsistency(h_passed, h_total):
        raise RuntimeError(f"TEfficiency consistency check failed for {name}.")

    efficiency = ROOT.TEfficiency(h_passed, h_total)
    efficiency.SetName(f"efficiency_{name}")
    efficiency.SetStatisticOption(ROOT.TEfficiency.kFCP)  # Clopper-Pearson

    # Keep the source histograms alive for the lifetime of the TEfficiency.
    return efficiency, h_passed, h_total, rows


def propagated_relative_error(acceptance, weighted_acceptance,
                              err_a_low, err_a_up,
                              err_aw_low, err_aw_up):
    """Return R, its lower error, and its upper error.

    First-order propagation for R = 100*abs(A - A_w)/A.

    The two estimates are assigned the same TEfficiency uncertainty, but are
    treated as uncorrelated here.  For asymmetric input intervals, the correct
    lower/upper side is selected according to the sign of each derivative.
    """
    if acceptance <= 0:
        raise ValueError("Acceptance must be positive to form the relative difference.")

    difference = acceptance - weighted_acceptance
    relative_difference = 100.0 * abs(difference) / acceptance

    # At exactly zero difference the absolute-value function has a cusp.  The
    # positive branch is used for the linearized uncertainty; the lower error
    # is subsequently limited so the plotted observable cannot go below zero.
    sign = 1.0 if difference >= 0.0 else -1.0
    derivative_a = 100.0 * sign * weighted_acceptance / (acceptance * acceptance)
    derivative_aw = -100.0 * sign / acceptance

    # Independent-error propagation with asymmetric intervals:
    # positive derivative: lower uses input lower error, upper uses upper;
    # negative derivative: the choices are reversed.
    low_variance = 0.0
    up_variance = 0.0
    for derivative, error_low, error_up in (
        (derivative_a, err_a_low, err_a_up),
        (derivative_aw, err_aw_low, err_aw_up),
    ):
        if derivative >= 0.0:
            low_variance += (derivative * error_low) ** 2
            up_variance += (derivative * error_up) ** 2
        else:
            low_variance += (derivative * error_up) ** 2
            up_variance += (derivative * error_low) ** 2

    error_low = math.sqrt(low_variance)
    error_up = math.sqrt(up_variance)

    # R is non-negative, so its lower error cannot exceed its central value.
    error_low = min(error_low, relative_difference)
    return relative_difference, error_low, error_up


def make_graph(rows, efficiency, name, color, marker_style):
    """Create a TGraphAsymmErrors for one rapidity interval."""
    graph = ROOT.TGraphAsymmErrors(len(rows))
    graph.SetName(f"graph_{name}")
    graph.SetTitle("")
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(1.05)
    graph.SetLineWidth(2)
    max_y_with_error = 0.0

    for point_index, row in enumerate(rows):
        root_bin = point_index + 1

        # TEfficiency supplies the unweighted acceptance and its binomial errors.
        acceptance = efficiency.GetEfficiency(root_bin)
        err_a_low = efficiency.GetEfficiencyErrorLow(root_bin)
        err_a_up = efficiency.GetEfficiencyErrorUp(root_bin)

        # First-order approximation requested: assign the same errors to A_w.
        weighted_acceptance = row["weighted_acceptance"]
        err_aw_low = err_a_low
        err_aw_up = err_a_up

        value, error_low, error_up = propagated_relative_error(
            acceptance, weighted_acceptance,
            err_a_low, err_a_up, err_aw_low, err_aw_up,
        )

        pt_center = 0.5 * (row["pt_low"] + row["pt_high"])
        pt_error = 0.5 * (row["pt_high"] - row["pt_low"])
        graph.SetPoint(point_index, pt_center, value)
        graph.SetPointError(point_index, pt_error, pt_error, error_low, error_up)
        max_y_with_error = max(max_y_with_error, value + error_up)

    return graph, max_y_with_error


def draw_group(rows, y_intervals, output_filename, canvas_name, error):
    """Draw three rapidity slices on one canvas and save as a PNG."""
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetEndErrorSize(4)

    slices = []
    all_y_max = 0.0
    for index, (y_low, y_high) in enumerate(y_intervals):
        slice_rows = select_rapidity_slice(rows, y_low, y_high)
        teff, h_passed, h_total, slice_rows = build_tefficiency(
            slice_rows, f"{canvas_name}_{index}"
        )
        color, marker_style = STYLES[index]
        graph, graph_y_max = make_graph(
            slice_rows, teff, f"{canvas_name}_{index}", color, marker_style
        )
        slices.append((graph, teff, h_passed, h_total, y_low, y_high))
        all_y_max = max(all_y_max, graph_y_max)

    pt_min = min(row["pt_low"] for _, _, _, _, y1, y2 in slices
                 for row in select_rapidity_slice(rows, y1, y2))
    pt_max = max(row["pt_high"] for _, _, _, _, y1, y2 in slices
                 for row in select_rapidity_slice(rows, y1, y2))
    y_max = max(1.0, 1.25 * all_y_max)

    canvas = ROOT.TCanvas(canvas_name, "", 800, 600)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)
    canvas.SetLogx()

    frame = canvas.DrawFrame(pt_min, -4.0, pt_max, y_max)
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("p_{T}(J/#psi) [GeV]")
    frame.GetYaxis().SetTitle("|A - A_{w}| / A  [%]")
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleOffset(1.15)
    frame.GetYaxis().SetTitleOffset(1.30)
    frame.GetXaxis().SetMoreLogLabels()

    legend = ROOT.TLegend(0.60, 0.69, 0.91, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)

    for graph, _, _, _, y_low, y_high in slices:
        graph.Draw("P" + error + "SAME")
        legend.AddEntry(graph, f"{y_low:.1f} #leq y < {y_high:.1f}", "lp")

    legend.Draw()

    note = ROOT.TLatex()
    note.SetNDC(True)
    note.SetTextSize(0.027)
    note.SetTextColor(ROOT.kGray + 2)
    note.DrawLatex(
        0.14, 0.93,
        "#sigma(A_{w}) = #sigma(A); zero covariance assumed",
    )

    canvas.RedrawAxis()
    canvas.SaveAs(output_filename)
    return canvas, slices


def main():
    parser = argparse.ArgumentParser(
        description="Plot relative acceptance differences with PyROOT."
    )
    parser.add_argument(
        "input_txt",
        help="Input acceptance table, e.g. Acceptance_reweighted.txt",
    )
    parser.add_argument(
        "--error",
        default="X",
        help="error bar option (default: no error bar)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display canvases interactively after saving the PNG files.",
    )
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(not args.show)
    rows = parse_table(args.input_txt)
    error = args.error

    output_files = (
        (RAPIDITY_GROUPS[0], "acc_y0p0_0p9.png", "c_y0p0_0p9"),
        (RAPIDITY_GROUPS[1], "acc_y0p9_2p4.png", "c_y0p9_2p4"),
    )

    # Keep canvases and graph/efficiency objects alive until the script exits.
    keep_alive = []
    for y_intervals, output_filename, canvas_name in output_files:
        canvas, slices = draw_group(
            rows, y_intervals, output_filename, canvas_name, error
        )
        keep_alive.append((canvas, slices))
        print(f"Saved: {output_filename}")

    if args.show:
        input("Press Enter to close the canvases... ")


if __name__ == "__main__":
    main()
