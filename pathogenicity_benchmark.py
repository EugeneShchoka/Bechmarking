from collections import defaultdict

from lib.annotation_lib import Annotation
from lib.baseutils import load_json, parse_json_lines, open_func
from lib.unchangable_variables import (
    EVIDENCE_CODE_LIST,
    PATHOGENICITY_MAPPING_SHRINKAGE,
)


def read_clingen(clingen_json_file):
    clingen_truthset_dict = {}
    indata = load_json(clingen_json_file)
    for entry in indata["data"]:
        clingen_truthset_dict[entry["identifier"]] = {
            "pathogenicity": entry["pathogenicity"],
            "evidence_codes": entry["evidence_codes"],
            "unmet_evidence_codes": entry["unmet_evidence_codes"],
        }
    return clingen_truthset_dict


def run_competitor_comparison(
        competitor_annotation_file,
        clingen_truthset_dict,
        merge_vus=False,
):
    pathogenicity_mapping = PATHOGENICITY_MAPPING_SHRINKAGE.copy()
    if merge_vus:
        pathogenicity_mapping["Uncertain significance P"] = "VUS"
        pathogenicity_mapping["Uncertain significance LP"] = "VUS"

    pathogenicity_compare_dict = defaultdict(int)
    evidence_code_dict_competitor = {
        ec: {
            "tp": 0,
            "fp": 0,
            "tn": 0,
            "fn": 0,
        }
        for ec in EVIDENCE_CODE_LIST
    }
    missing = 0
    double_counting_dict = defaultdict(int)
    for entry in open_func(competitor_annotation_file, to_dict=True):
        # possible string is: {'Variant': 'chr1:171636330 C⇒T', 'Chromosome': 'chr1', 'Position': '171636330', 'RS ID': 'rs149881467', 'Ref seq': 'C' etc.
        if entry["Ref seq"] == "":
            entry["Ref seq"] = "."
        if entry["Var seq"] == "":
            entry["Var seq"] = "."
        if entry["Chromosome"].startswith("chr"):
            entry["Chromosome"] = entry["Chromosome"][3:]
        identifier = "{}-{}-{}-{}".format(
            entry["Chromosome"], entry["Position"], entry["Ref seq"], entry["Var seq"]
        )
        if identifier in clingen_truthset_dict:
            compare_id = (
                clingen_truthset_dict[identifier]["pathogenicity"],
                pathogenicity_mapping[entry["Germline Class"]],
            )
            pathogenicity_compare_dict[compare_id] += 1

            evidence_codes_competitor = [
                ec.split("_")[0] for ec in entry["Germline rules"].split(",")
            ]

            clingen_evidence_codes = clingen_truthset_dict[identifier]["evidence_codes"]
            clingen_evidence_codes_unmet = clingen_truthset_dict[identifier][
                "unmet_evidence_codes"
            ]

            for ec in evidence_code_dict_competitor:
                if ec in evidence_codes_competitor and ec in clingen_evidence_codes:
                    evidence_code_dict_competitor[ec]["tp"] += 1
                elif (
                        ec in evidence_codes_competitor and ec in clingen_evidence_codes_unmet
                ):
                    evidence_code_dict_competitor[ec]["fp"] += 1
                elif ec not in evidence_codes_competitor and ec in clingen_evidence_codes:
                    evidence_code_dict_competitor[ec]["fn"] += 1
                elif (
                        ec not in evidence_codes_competitor
                        and ec in clingen_evidence_codes_unmet
                ):
                    evidence_code_dict_competitor[ec]["tn"] += 1

            # check whether the seq_evidence_codes_competitor contains at least 2 of elemtnis in PS1, PP5 and PM5
            clingen_related_codes = sorted(
                {"PP5", "PS1", "PM5"}.intersection(evidence_codes_competitor)
            )
            if len(clingen_related_codes) >= 2:
                double_counting_dict[tuple(clingen_related_codes)] += 1
        else:
            missing += 1

    return (
        pathogenicity_compare_dict,
        evidence_code_dict_competitor,
        double_counting_dict,
        missing,
    )


def compare_seq_vs_clingen(seq_annotation_json, merge_vus=False):
    pathogenicity_compare_dict = defaultdict(int)
    evidence_code_dict = {
        ec: {
            "tp": 0,
            "fp": 0,
            "tn": 0,
            "fn": 0,
        }
        for ec in EVIDENCE_CODE_LIST
    }
    for entry in parse_json_lines(seq_annotation_json):
        annot = Annotation(entry)
        if not annot.clingen_entries:
            continue
        clingen_pathogeniciy = annot.clingen_entries[0].pathogenicity
        clingen_evidence_codes = [
            ec.split("_")[0] for ec in annot.clingen_entries[0].evidence_codes
        ]
        clingen_evidence_codes_unmet = [
            ec.split("_")[0] for ec in annot.clingen_entries[0].unmet_evidence_codes
        ]
        seq_evidence_codes = [
            ec.replace("+", "").replace("-", "")
            for ec in entry["annotations"]["transcript"]["acmg_evidence_codes"]
            if ec != "BA1-"
        ]

        for ec in EVIDENCE_CODE_LIST:
            if ec in seq_evidence_codes and ec in clingen_evidence_codes:
                evidence_code_dict[ec]["tp"] += 1
            elif ec in seq_evidence_codes and ec in clingen_evidence_codes_unmet:
                evidence_code_dict[ec]["fp"] += 1
            elif ec not in seq_evidence_codes and ec in clingen_evidence_codes:
                evidence_code_dict[ec]["fn"] += 1
            elif ec not in seq_evidence_codes and ec in clingen_evidence_codes_unmet:
                evidence_code_dict[ec]["tn"] += 1
        autopat_code = annot.autopat_code.replace("-", "")
        if merge_vus:
            autopat_code = autopat_code.replace("+", "")
        compare_id = (
            clingen_pathogeniciy,
            autopat_code,
        )
        pathogenicity_compare_dict[compare_id] += 1
    return pathogenicity_compare_dict, evidence_code_dict
