import os
import pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from pathogenicity_benchmark import (
    read_clingen,
    run_competitor_comparison,
    compare_seq_vs_clingen,
)


def prepare_sankey_data(data_dict, mapping1, mapping2, node_colors):
    sources = []
    targets = []
    counts = []
    colors = []

    for key, val in data_dict.items():
        sources.append(mapping1[key[0]])
        targets.append(mapping2[key[1]])
        counts.append(val)
        colors.append(f"rgba({node_colors[mapping1[key[0]]]}, 0.5)")

    return sources, targets, counts, colors


def create_sankey_subplot(
        fig,
        sources,
        targets,
        counts,
        colors,
        node_labels,
        node_colors,
        titles,
        subplot_index,
):
    # convert to RGB
    node_colors = [f"rgb({color})" for color in node_colors] * 2

    fig.add_trace(
        go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=node_labels,
                color=node_colors,
            ),
            link=dict(
                source=sources,
                target=targets,
                value=counts,
                color=colors,
            ),
            domain=dict(
                x=[subplot_index * (0.3 + 0.05), subplot_index * (0.3 + 0.05) + 0.3],
                y=[0, 1],
            ),
        )
    )

    # https://stackoverflow.com/questions/67540925/plotly-how-to-write-a-text-over-my-sankey-diagram-columns
    cols = titles
    for x_coordinate, column_name in enumerate(cols):
        fig.add_annotation(
            x=subplot_index * (0.3 + 0.08) + 0.25 * x_coordinate,
            y=1.05,
            xref="paper",
            yref="paper",
            text=column_name,
            showarrow=False,
            font=dict(
                family="Courier New, monospace",
                size=16,
            ),
        )


def create_sankey_figure(fig, data_dict, clingen_mapping, competitor_mapping, node_colors, node_labels, titles,
                         subplot_index):
    sources, targets, counts, colors = prepare_sankey_data(
        data_dict,
        clingen_mapping,
        competitor_mapping,
        node_colors,
    )
    create_sankey_subplot(
        fig,
        sources,
        targets,
        counts,
        colors,
        node_labels,
        node_colors,
        titles,
        subplot_index,
    )


def dict_to_tsv(comparison_dict, index_mapping, filename_output):
    # Create a reverse mapping dictionary
    reverse_competitor_mapping = {v: k for k, v in index_mapping.items()}

    # Create empty DataFrame
    df = pd.DataFrame(
        columns=["Clingen pathogenicity", "Predicted pathogenicity", "Variants counts"]
    )

    # Process and fill the DataFrame
    for key, value in comparison_dict.items():
        merged_key = (index_mapping[key[0]], index_mapping[key[1]])
        clingen_name = reverse_competitor_mapping[merged_key[0]]
        predicted_name = reverse_competitor_mapping[merged_key[1]]
        # check if clingen_name, predicted_name pair already exists
        if df[
            (df["Clingen pathogenicity"] == clingen_name)
            & (df["Predicted pathogenicity"] == predicted_name)
        ].empty:
            df.loc[len(df.index)] = [clingen_name, predicted_name, value]
        else:
            index = df[
                (df["Clingen pathogenicity"] == clingen_name)
                & (df["Predicted pathogenicity"] == predicted_name)
                ].index[0]
            df.at[index, "Variants counts"] += value

    # Save as TSV
    df.to_csv(os.path.join("data", "output", filename_output), sep="\t", index=False)


def calculate_statistics_from_tsv(tsv_file, filename_output, merged=False):
    df = pd.read_csv(os.path.join("data", "output", tsv_file), sep="\t")
    result_df = pd.DataFrame(columns=["F1", "Precision", "Recall"])
    pathogenicity_values = ["P", "LP", "VUS", "LB", "B"]
    if merged:
        pathogenicity_values = ["P", "VUS", "B"]

    for patho in pathogenicity_values:
        relevant_data = df

        true_positives = relevant_data[
            (relevant_data["Predicted pathogenicity"] == patho)
            & (relevant_data["Clingen pathogenicity"] == patho)
            ]["Variants counts"].sum()

        false_positives = relevant_data[
            (relevant_data["Predicted pathogenicity"] == patho)
            & (relevant_data["Clingen pathogenicity"] != patho)
            ]["Variants counts"].sum()

        false_negatives = relevant_data[
            (relevant_data["Predicted pathogenicity"] != patho)
            & (relevant_data["Clingen pathogenicity"] == patho)
            ]["Variants counts"].sum()

        # Calculation of metrics (assuming you already have the above)
        if true_positives + false_positives > 0:
            precision = true_positives / (true_positives + false_positives)
            precision = round(precision, 3)
        else:
            precision = 0  # Handle division by zero

        if true_positives + false_negatives > 0:
            recall = true_positives / (true_positives + false_negatives)
            recall = round(recall, 3)
        else:
            recall = 0

        f1 = 2 * (precision * recall) / (precision + recall)
        f1 = round(f1, 3)

        result_df.loc[patho] = [f1, precision, recall]

    # save to data folder
    result_df.to_csv(os.path.join("data", "output", filename_output), sep="\t")

    return result_df


def radar_chart(seq_ensembl_stat, seq_refseq_stat, competitor_stat, filename, merged=False):
    # set same grid for every subplot
    fig_radar = make_subplots(rows=1, cols=3, specs=[[{"type": "polar"}] * 3] * 1)

    categories = ["P", "LP", "VUS", "LB", "B"]
    if merged:
        categories = ["P", "VUS", "B"]

    # Define colors for each group
    color_map = {"SEQ-Ensembl": "blue", "SEQ-RefSeq": "lightgreen", "Competitor": "salmon"}

    for group, color in color_map.items():
        fig_radar.add_trace(
            go.Scatterpolar(
                r=(
                    seq_ensembl_stat.loc[categories]["F1"].values
                    if group == "SEQ-Ensembl"
                    else (
                        seq_refseq_stat.loc[categories]["F1"].values
                        if group == "SEQ-RefSeq"
                        else competitor_stat.loc[categories]["F1"].values
                    )
                ),
                theta=categories,
                fill="toself",
                name=group,
                marker=dict(color=color),  # Assign color based on the group
                legendgroup=group,
                showlegend=True,
            ),
            row=1,
            col=1,
        )

        fig_radar.add_trace(
            go.Scatterpolar(
                r=(
                    seq_ensembl_stat.loc[categories]["Recall"].values
                    if group == "SEQ-Ensembl"
                    else (
                        seq_refseq_stat.loc[categories]["Recall"].values
                        if group == "SEQ-RefSeq"
                        else competitor_stat.loc[categories]["Recall"].values
                    )
                ),
                theta=categories,
                fill="toself",
                name=group,
                marker=dict(color=color),  # Assign color based on the group
                legendgroup=group,
                showlegend=False,
            ),
            row=1,
            col=2,
        )

        fig_radar.add_trace(
            go.Scatterpolar(
                r=(
                    seq_ensembl_stat.loc[categories]["Precision"].values
                    if group == "SEQ-Ensembl"
                    else (
                        seq_refseq_stat.loc[categories]["Precision"].values
                        if group == "SEQ-RefSeq"
                        else competitor_stat.loc[categories]["Precision"].values
                    )
                ),
                theta=categories,
                fill="toself",
                name=group,
                marker=dict(color=color),  # Assign color based on the group
                legendgroup=group,
                showlegend=False,
            ),
            row=1,
            col=3,
        )

    fig_radar.update_layout(
        autosize=False,
        width=1280,
        height=720,
        polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
        polar2=dict(radialaxis=dict(visible=True, range=[0, 1])),
        polar3=dict(radialaxis=dict(visible=True, range=[0, 1])),
        showlegend=True,
    )

    fig_radar.add_annotation(
        x=0.12,
        y=1.1,
        xref="paper",
        yref="paper",
        text="F1",
        showarrow=False,
        font=dict(
            family="Courier New, monospace",
            size=20,
        ),
    )
    fig_radar.add_annotation(
        x=0.48,
        y=1.1,
        xref="paper",
        yref="paper",
        text="Recall",
        showarrow=False,
        font=dict(
            family="Courier New, monospace",
            size=20,
        ),
    )
    fig_radar.add_annotation(
        x=0.85,
        y=1.1,
        xref="paper",
        yref="paper",
        text="Precision",
        showarrow=False,
        font=dict(
            family="Courier New, monospace",
            size=20,
        ),
    )

    fig_radar.write_image(os.path.join("data", "output", filename), width=1280, height=720)


def main(
        clingen_json_file,
        seq_refseq_annotation_file,
        seq_ensembl_annotation_file,
        competitor_annotation_file,
):
    # check if output folder exists
    if not os.path.exists(os.path.join("data", "output")):
        os.makedirs(os.path.join("data", "output"))

    # process data
    clingen_truthset_dict = read_clingen(clingen_json_file)
    (
        seq_ensembl_pathogenicity_comparison_dict,
        seq_ensembl_evidence_code_comparison_dict,
    ) = compare_seq_vs_clingen(seq_ensembl_annotation_file)
    (
        seq_refseq_pathogenicity_comparison_dict,
        seq_refseq_evidence_code_comparison_dict,
    ) = compare_seq_vs_clingen(seq_refseq_annotation_file)
    (
        competitor_pathogenicity_comparison_dict,
        competitor_evidence_code_comparison_dict,
        competitor_double_counting_dict,
        missing,
    ) = run_competitor_comparison(
        competitor_annotation_file, clingen_truthset_dict, merge_vus=True
    )

    # define index mappings
    clingen_index_mapping = {"P": 0, "LP": 1, "VUS": 2, "LB": 3, "B": 4}
    competitor_index_mapping = {
        "P": 5,
        "LP": 6,
        "VUS++": 7,
        "VUS+": 7,
        "VUS": 7,
        "LB": 8,
        "B": 9,
    }
    clingen_index_mapping_merged = {
        "LP": 0,
        "P": 0,
        "VUS": 1,
        "LB": 2,
        "B": 2,
    }  # LP should be before P
    competitor_index_mapping_merged = {
        "LP": 3,
        "P": 3,
        "VUS++": 4,
        "VUS+": 4,
        "VUS": 4,
        "LB": 5,
        "B": 5,
    }  # LP should be before P

    # define node colors
    node_colors = [
        "245, 132, 98",
        "248, 182, 165",
        "90, 169, 218",
        "151, 210, 178",
        "92, 189, 123",
    ]  # "#f58462", "#f8b6a5", "#5aa9da", "#97d2b2", "#5cbd7b"
    node_colors_merged = [
        "245, 132, 98",
        "90, 169, 218",
        "92, 189, 123",
    ]  # "#f58462", "#5aa9da", "#5cbd7b"

    # define graph labels
    node_labels = ["P", "LP", "VUS", "LB", "B", "P", "LP", "VUS", "LB", "B"]
    node_labels_merged = ["P/LP", "VUS", "B/LB", "P/LP", "VUS", "B/LB"]

    # Create a figure with subplots
    fig_comparison = go.Figure()

    create_sankey_figure(
        fig_comparison,
        seq_ensembl_pathogenicity_comparison_dict,
        clingen_index_mapping,
        competitor_index_mapping,
        node_colors,
        node_labels,
        ["Clingen", "SEQ-Ensembl"],
        0,
    )

    create_sankey_figure(
        fig_comparison,
        seq_refseq_pathogenicity_comparison_dict,
        clingen_index_mapping,
        competitor_index_mapping,
        node_colors,
        node_labels,
        ["Clingen", "SEQ-RefSeq"],
        1,
    )

    create_sankey_figure(
        fig_comparison,
        competitor_pathogenicity_comparison_dict,
        clingen_index_mapping,
        competitor_index_mapping,
        node_colors,
        node_labels,
        ["Clingen", "Competitor"],
        2,
    )

    fig_comparison.add_annotation(
        x=0.5,
        y=1.1,
        xref="paper",
        yref="paper",
        text="Comparison of Pathogenicity Prediction",
        showarrow=False,
        font=dict(
            family="Courier New, monospace",
            size=20,
        ),
    )
    fig_comparison.write_image(os.path.join("data", "output", "sankey_diagram.pdf"), width=1280, height=720)

    # Create a figure with merged subplots
    fig_comparison_merged = go.Figure()

    create_sankey_figure(
        fig_comparison_merged,
        seq_ensembl_pathogenicity_comparison_dict,
        clingen_index_mapping_merged,
        competitor_index_mapping_merged,
        node_colors_merged,
        node_labels_merged,
        ["Clingen", "SEQ-Ensembl"],
        0,
    )

    create_sankey_figure(
        fig_comparison_merged,
        seq_refseq_pathogenicity_comparison_dict,
        clingen_index_mapping_merged,
        competitor_index_mapping_merged,
        node_colors_merged,
        node_labels_merged,
        ["Clingen", "SEQ-RefSeq"],
        1,
    )

    create_sankey_figure(
        fig_comparison_merged,
        competitor_pathogenicity_comparison_dict,
        clingen_index_mapping_merged,
        competitor_index_mapping_merged,
        node_colors_merged,
        node_labels_merged,
        ["Clingen", "Competitor"],
        2,
    )

    fig_comparison_merged.add_annotation(
        x=0.5,
        y=1.1,
        xref="paper",
        yref="paper",
        text="Comparison of Pathogenicity Prediction",
        showarrow=False,
        font=dict(
            family="Courier New, monospace",
            size=20,
        ),
    )
    fig_comparison_merged.write_image(os.path.join("data", "output", "sankey_diagram_merged.pdf"), width=1280, height=720)

    # Save the comparison data as TSV
    dict_to_tsv(
        seq_ensembl_pathogenicity_comparison_dict,
        competitor_index_mapping,
        "seq_ensembl.tsv",
    )
    dict_to_tsv(
        seq_refseq_pathogenicity_comparison_dict,
        competitor_index_mapping,
        "seq_refseq.tsv",
    )
    dict_to_tsv(
        competitor_pathogenicity_comparison_dict,
        clingen_index_mapping,
        "competitor.tsv"
    )
    dict_to_tsv(
        seq_ensembl_pathogenicity_comparison_dict,
        competitor_index_mapping_merged,
        "seq_ensembl_merged.tsv",
    )
    dict_to_tsv(
        seq_refseq_pathogenicity_comparison_dict,
        competitor_index_mapping_merged,
        "seq_refseq_merged.tsv",
    )
    dict_to_tsv(
        competitor_pathogenicity_comparison_dict,
        clingen_index_mapping_merged,
        "competitor_merged.tsv",
    )

    # Calculate statistics
    seq_ensembl_stat = calculate_statistics_from_tsv(
        "seq_ensembl.tsv",
        "seq_ensembl_statistics.tsv"
    )
    seq_refseq_stat = calculate_statistics_from_tsv(
        "seq_refseq.tsv",
        "seq_refseq_statistics.tsv"
    )
    competitor_stat = calculate_statistics_from_tsv(
        "competitor.tsv",
        "competitor_statistics.tsv"
    )
    seq_ensembl_stat_merged = calculate_statistics_from_tsv(
        "seq_ensembl_merged.tsv",
        "seq_ensembl_merged_statistics.tsv",
        merged=True,
    )
    seq_refseq_stat_merged = calculate_statistics_from_tsv(
        "seq_refseq_merged.tsv",
        "seq_refseq_merged_statistics.tsv",
        merged=True,
    )
    competitor_stat_merged = calculate_statistics_from_tsv(
        "competitor_merged.tsv",
        "competitor_merged_statistics.tsv",
        merged=True
    )

    radar_chart(
        seq_ensembl_stat,
        seq_refseq_stat,
        competitor_stat,
        "radar_chart.pdf"
    )
    radar_chart(
        seq_ensembl_stat_merged,
        seq_refseq_stat_merged,
        competitor_stat_merged,
        "radar_chart_merged.pdf",
        merged=True,
    )


if __name__ == "__main__":
    clingen_json_file_path = os.path.join("data", "clingen", "clingen_variant_hg38.json.gz")
    seq_ensembl_annotation_file_path = os.path.join("data", "seq", "clingen_annotation_2024feb_ensembl.json.gz")
    seq_refseq_annotation_file_path = os.path.join("data", "seq", "clingen_annotation_2024feb_refseq.json.gz")
    competitor_annotation_file_path = (
        os.path.join("data", "competitor", "competitor_clingen2024feb_hg38_variant_annotation.tsv.gz")
    )
    main(
        clingen_json_file_path,
        seq_refseq_annotation_file_path,
        seq_ensembl_annotation_file_path,
        competitor_annotation_file_path,
    )
