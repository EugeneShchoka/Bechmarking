EVIDENCE_CODE_LIST = [
    "PVS1",
    "PS1",
    "PM1",
    "PM2",
    "PM4",
    "PM5",
    "PP2",
    "PP3",
    "BA1",
    "BS1",
    "BS2",
    "BP3",
    "BP4",
    "BP7",
]

PATHOGENICITY_MAPPING_SHRINKAGE = {
    "Pathogenic": "P",
    "Likely pathogenic": "LP",
    "Uncertain significance P": "VUS++",
    "Uncertain significance LP": "VUS+",
    "Uncertain significance": "VUS",
    "Likely benign": "LB",
    "Benign": "B",
}

PATHOGENICITY_MAPPING_EXTENDED = {
    "P": "P",
    "LP": "LP",
    "VUS": "VUS",
    "LB": "LB",
    "B": "B",
    "P-": "P",
    "P--": "P",
    "P---": "P",
    "LP-": "LP",
    "LP--": "LP",
    "LP---": "LP",
    "VUS++": "VUS++",
    "VUS++-": "VUS++",
    "VUS++--": "VUS++",
    "VUS++---": "VUS++",
    "VUS+": "VUS+",
    "VUS+-": "VUS+",
    "VUS+--": "VUS+",
    "VUS+---": "VUS+",
    "VUS-": "VUS",
    "VUS--": "VUS",
}