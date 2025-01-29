schema = {
    "type": "object",
    "properties": {
        "missingness_threshold": {"type": "number"},
        "indexed_array_annotations": {
            "type": "object",
            "additionalProperties": {"type": "string"}  # Allow any string keys with string values
        },
        "struct_annotations_for_missingness": {
            "type": "array",
            "items": {"type": "string"}
        },
        "freq_meta_expr": {"type": "string"},
        "freq_annotations_to_sum": {
            "type": "array",
            "items": {"type": "string"}
        },
        "freq_sort_order": {
            "type": "array",
            "items": {"type": "string"}
        },
        "nhomalt_metric": {"type": "string"},
    },
    "required": [
        "missingness_threshold",
        "indexed_array_annotations",
        "struct_annotations_for_missingness",
        "freq_meta_expr",
        "freq_annotations_to_sum",
        "freq_sort_order",
        "nhomalt_metric"
    ],
    "additionalProperties": False
}