import re
import hail as hl


def parse_log_file(log_file):
    """Parses a log file and categorizes messages for formatting, extracting function names."""
    parsed_logs = []
    log_pattern = re.compile(
        r"\(([^)]+)\.(\w+) (\d+)\)"
    )  # Extract module and function name

    with hl.hadoop_open(log_file, "r") as f:
        for line in f:
            match = log_pattern.search(line)
            function_name = match.group(2) if match else "Unknown"

            # Determine the category based on log message content
            if "PASSED" in line:
                category = "pass"
            elif "FAILED" in line:
                category = "fail"
            elif "WARNING" in line:
                category = "warn"
            else:
                category = "info"

            parsed_logs.append((function_name, category, line.strip()))

    return parsed_logs


def generate_html_report(parsed_logs, output_file):
    """Generates an HTML report with sortable and filterable columns for function name, status, and message."""
    html_template = """
    <html>
    <head>
        <style>
            body { font-family: Arial, sans-serif; }
            table { width: 100%; border-collapse: collapse; }
            th, td { border: 1px solid black; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; cursor: pointer; }
            .pass { color: green; }
            .fail { color: red; font-weight: bold; }
            .warn { color: orange; }
            .info { color: blue; }
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
        <label for="functionFilter">Filter by Function Name:</label>
        <select id="functionFilter" onchange="filterTable(0, this.value)">
            <option value="all">All</option>
    """

    function_names = set()
    statuses = set()
    for function_name, category, message in parsed_logs:
        function_names.add(function_name)
        statuses.add(category)

    for function_name in sorted(function_names):
        html_template += f'<option value="{function_name}">{function_name}</option>'

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
                <th onclick="sortTable(0)">Function Name</th>
                <th onclick="sortTable(1)">Status</th>
                <th onclick="sortTable(2)">Message</th>
            </tr>
    """

    for function_name, category, message in parsed_logs:
        html_template += f'<tr class="{category}"><td>{function_name}</td><td class="{category}">{category.upper()}</td><td>{message}</td></tr>'

    html_template += """
        </table>
    </body>
    </html>
    """

    with hl.hadoop_open(output_file, "w") as f:
        f.write(html_template)


if __name__ == "__main__":
    log_file = "log_output.txt"  # Change this to your log file
    output_html = "log_report.html"

    parsed_logs = parse_log_file(log_file)
    generate_html_report(parsed_logs, output_html)

    print(f"HTML report generated: {output_html}")
