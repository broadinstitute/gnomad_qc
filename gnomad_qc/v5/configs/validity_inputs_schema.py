"""Script defining schema for json config that is input to the federated validity checks."""

schema = {
    "type": "object",
    "properties": {
        "missingness_threshold": {"type": "number"},
        "struct_annotations_for_missingness": {
            "type": "array",
            "items": {"type": "string"},
        },
        "freq_fields": {
            "type": "object",
            "properties": {
                "freq": {"type": "string"},
                "freq_meta": {"type": "string"},
                "freq_index_dict": {"type": "string"},
                "freq_meta_sample_count": {"type": "string"},
            },
            "required": [
                "freq",
                "freq_meta",
                "freq_meta_sample_count",
            ],
            "additionalProperties": False,
        },
        "faf_fields": {
            "type": "object",
            "properties": {
                "faf": {"type": "string"},
                "faf_meta": {"type": "string"},
                "faf_index_dict": {"type": "string"},
            },
            "additionalProperties": False,
        },
        "freq_annotations_to_sum": {
            "type": "array",
            "items": {"type": "string"},
        },
        "sort_order": {
            "type": "array",
            "items": {"type": "string"},
        },
        "nhomalt_metric": {"type": "string"},
        "subsets": {
            "type": "array",
            "items": {"type": "string"},
        },
        "variant_filter_field": {"type": "string"},
    },
    "required": [
        "missingness_threshold",
        "freq_fields",
        "freq_annotations_to_sum",
        "sort_order",
        "nhomalt_metric",
        "subsets",
        "variant_filter_field",
    ],
    "additionalProperties": False,
}
