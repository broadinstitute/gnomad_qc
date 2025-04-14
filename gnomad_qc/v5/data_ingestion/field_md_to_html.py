import markdown
from bs4 import BeautifulSoup

# Read markdown content.
with open("field_requirements.md", "r") as f:
    markdown_text = f.read()

# Convert markdown to HTML.
html = markdown.markdown(markdown_text, extensions=["extra"])
soup = BeautifulSoup(html, "html.parser")

# Color map based on Field Necessity.
color_map = {
    "required": "#fff0f0",
    "optional": "#f0faff",
    "not_needed": "#dedbe4",
}

final_sections = []
first_table_processed = False
second_table_processed = False

elements = soup.find_all(["p", "table"])

paragraphs_before_first_table = []
paragraphs_before_second_table = []

for element in elements:
    if element.name == "p":
        if not first_table_processed:
            paragraphs_before_first_table.append(str(element))
        elif not second_table_processed:
            paragraphs_before_second_table.append(str(element))

    elif element.name == "table":
        if not first_table_processed:
            first_table_processed = True
            final_sections.append("<h1>Global Fields Specification</h1>")
            final_sections.extend(paragraphs_before_first_table)

        elif not second_table_processed:
            second_table_processed = True
            final_sections.append("<h1>Row Fields Specification</h1>")
            final_sections.extend(paragraphs_before_second_table)

        # Build a styled table
        styled_table = BeautifulSoup(
            "<table style='border-collapse: collapse; border: 1px solid #000; width: 100%;'><thead></thead><tbody></tbody></table>",
            "html.parser",
        )
        thead = styled_table.find("thead")
        tbody = styled_table.find("tbody")

        rows = element.find_all("tr")
        if not rows:
            continue

        # Header row
        header_cells = rows[0].find_all(["th", "td"])
        header_row = styled_table.new_tag("tr", style="background-color:#f8f8f8;")
        for cell in header_cells:
            th = styled_table.new_tag(
                "th", style="border: 1px solid #000; padding: 5px; text-align: left;"
            )
            th.append(BeautifulSoup(cell.decode_contents(), "html.parser"))
            header_row.append(th)
        thead.append(header_row)

        # Data rows.
        for row in rows[1:]:
            cells = row.find_all("td")
            if len(cells) < 5:
                continue

            necessity = cells[4].get_text(strip=True).lower()
            bg_color = color_map.get(necessity, "#ffffff")
            tr = styled_table.new_tag("tr", style=f"background-color:{bg_color};")

            for cell in cells:
                td = styled_table.new_tag(
                    "td", style="border: 1px solid #000; padding: 5px;"
                )

                # Check if cell content is not empty.
                if cell.decode_contents():
                    td.append(BeautifulSoup(cell.decode_contents(), "html.parser"))
                else:
                    # Handle empty cells by appending a blank space or placeholder.
                    td.append(" ")

                tr.append(td)

            tbody.append(tr)

        final_sections.append(str(styled_table))

final_html = """
<style>
  table {
    border-collapse: collapse;
  }
  th, td {
    border: 1px solid #000;
    padding: 5px;
  }
</style>
""" + "\n".join(
    final_sections
)

# Write final output.
with open("field_requirements.html", "w") as out:
    out.write(final_html)

print("✅ HTML file saved with correct spans and borders.")
