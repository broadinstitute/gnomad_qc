import re
import hail as hl


def parse_log_file(log_file):
    """Parses a log file and categorizes messages for formatting, extracting function names and sources."""
    parsed_logs = []
    log_pattern = re.compile(
        r"(\w+) \(([^)]+)\.(\w+) (\d+)\): (.*)"
    )  # Extract log level, module, function, line number, and message

    function_mapping = {
        "validate_config": "general info",
        "validate_ht_fields": "general info",
        "main": "general info",
        "validate_federated_data": "general info",
        "summarize_variants": "variant summary",
        "sum_group_callstats": "group summations",
        "make_group_sum_expr_dict": "group summations",
        "check_sex_chr_metrics": "XY check",
        "generic_field_check": "raw/adj check",
    }

    with hl.hadoop_open(log_file, "r") as f:
        for line in f:
            match = log_pattern.match(line)
            if match:
                level, module, function_name, line_number, message = match.groups()
                source = (
                    f"{module}.{function_name} {line_number}"  # Create a source column
                )
                validity_check = function_mapping.get(
                    function_name, function_name
                )  # Map function names

                # Determine category based on log message
                if "PASSED" in message:
                    category = "pass"
                elif "FAILED" in message:
                    category = "fail"
                elif "WARNING" in message:
                    category = "warn"
                else:
                    category = "info"

                parsed_logs.append((validity_check, category, source, message))

    return parsed_logs


def generate_html_report(parsed_logs, output_file):
    """Generates an HTML report with sortable and filterable columns for validity check, status, source, and message."""
    html_template = """
    <html>
    <head>
        <style>
            body { font-family: Arial, sans-serif; }
            table { width: 100%; border-collapse: collapse; }
            th, td { border: 1px solid black; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; cursor: pointer; }
            .pass { color: #6d1faa; }
            .fail { color: #D42736; font-weight: bold; }
            .warn { color: #DAA520; }
            .info { color: black; }
        </style>
        <script>
            function sortTable(n) {
                var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
                table = document.getElementById("logTable");
                switching = true;
                dir = "asc";
                while (switching) {
                    switching = false;
                    rows = table.rows;
                    for (i = 1; i < (rows.length - 1); i++) {
                        shouldSwitch = false;
                        x = rows[i].getElementsByTagName("TD")[n].innerHTML.toLowerCase();
                        y = rows[i + 1].getElementsByTagName("TD")[n].innerHTML.toLowerCase();
                        if ((dir == "asc" && x > y) || (dir == "desc" && x < y)) {
                            shouldSwitch = true;
                            break;
                        }
                    }
                    if (shouldSwitch) {
                        rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                        switching = true;
                        switchcount++;
                    } else {
                        if (switchcount === 0 && dir === "asc") {
                            dir = "desc";
                            switching = true;
                        }
                    }
                }
            }
            
            function filterTable(column, value) {
                var table, tr, td, i;
                table = document.getElementById("logTable");
                tr = table.getElementsByTagName("tr");
                for (i = 1; i < tr.length; i++) {
                    td = tr[i].getElementsByTagName("td")[column];
                    if (td) {
                        if (value === "all" || td.innerHTML.toLowerCase() === value.toLowerCase()) {
                            tr[i].style.display = "";
                        } else {
                            tr[i].style.display = "none";
                        }
                    } 
                }
            }
        </script>
    </head>
    <body>
        <h2>Log Report</h2>
        <label for="functionFilter">Filter by Validity Check:</label>
        <select id="functionFilter" onchange="filterTable(0, this.value)">
            <option value="all">All</option>
    """

    validity_checks = set()
    statuses = set()
    for validity_check, category, source, message in parsed_logs:
        validity_checks.add(validity_check)
        statuses.add(category)

    for validity_check in sorted(validity_checks):
        html_template += f'<option value="{validity_check}">{validity_check}</option>'

    html_template += """
        </select>
        <label for="statusFilter">Filter by Status:</label>
        <select id="statusFilter" onchange="filterTable(1, this.value)">
            <option value="all">All</option>
    """

    for status in sorted(statuses):
        html_template += f'<option value="{status}">{status.upper()}</option>'

    html_template += """
        </select>
        <table id="logTable">
            <tr>
                <th onclick="sortTable(0)">Validity Check</th>
                <th onclick="sortTable(1)">Status</th>
                <th onclick="sortTable(2)">Source</th>
                <th onclick="sortTable(3)">Message</th>
            </tr>
    """

    for validity_check, category, source, message in parsed_logs:
        html_template += f'<tr class="{category}"><td>{validity_check}</td><td class="{category}">{category.upper()}</td><td>{source}</td><td>{message}</td></tr>'

    html_template += """
        </table>
    </body>
    </html>
    """

    with hl.hadoop_open(output_file, "w") as f:
        f.write(html_template)
