"""Script defining schema for federated validty check input config file."""

schema = {
    "type": "object",
    "properties": {
        "missingness_threshold": {"type": "number"},
        "struct_annotations_for_missingness": {
            "type": "array",
            "items": {"type": "string"},
        },
        "freq_names": {
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
                "freq_index_dict",
                "freq_meta_sample_count",
            ],
            "additionalProperties": False,
        },
        "faf_names": {
            "type": "object",
            "properties": {
                "faf": {"type": "string"},
                "faf_meta": {"type": "string"},
                "faf_index_dict": {"type": "string"},
            },
            "required": ["faf", "faf_meta", "faf_index_dict"],
            "additionalProperties": False,
        },
        "freq_annotations_to_sum": {
            "type": "array",
            "items": {"type": "string"},
        },
        "freq_sort_order": {
            "type": "array",
            "items": {"type": "string"},
        },
        "nhomalt_metric": {"type": "string"},
    },
    "required": [
        "missingness_threshold",
        "struct_annotations_for_missingness",
        "freq_names",
        "faf_names",
        "freq_annotations_to_sum",
        "freq_sort_order",
        "nhomalt_metric",
    ],
    "additionalProperties": False,
}
