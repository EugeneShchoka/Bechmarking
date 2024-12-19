from lib.unchangable_variables import PATHOGENICITY_MAPPING_EXTENDED


class ClingenVariant:
    def __init__(self, entry):
        self.mondo_id = entry["mondo_id"]
        self.disease_id = f'MONDO:{entry["mondo_id"]}'
        self.disease_name = entry["disease_name"]
        self.mode_of_inheritances = entry["mode_of_inheritances"]
        self.pathogenicity = entry["pathogenicity"]
        self.evidence_codes = entry["evidence_codes"]
        self.unmet_evidence_codes = entry.get("unmet_evidence_codes", [])
        self.variant_id = entry["variant_id"]

        self.gene_ids = entry.get("gene_ids")
        self.transcript_ids = entry.get("transcript_ids", [])
        self.aa_change = entry.get("aa_change")
        self.identifier = entry.get("identifier")
        self.clinvar_id = entry.get("clinvar_id")
        self.hgvsg = entry.get("hgvsg")
        self.hgvs_all = entry.get("hgvs_all")
        self.gene_symbol = entry.get("gene_symbol")
        self.gene_unique_ids = entry.get("unique_gene_ids")


class Annotation:
    def __init__(self, data):
        self.data = data
        self.clingen_data = data["annotations"]["variant"].get("clingen", [])
        self.clingen_entries = [ClingenVariant(entry) for entry in self.clingen_data]
        self.autopat_code = PATHOGENICITY_MAPPING_EXTENDED[
            data["annotations"]["transcript"]["auto_pathogenicity"]
        ]
