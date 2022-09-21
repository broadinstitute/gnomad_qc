"""
Script that pulls Terra workspace owners.

Script reads through input file (with one Terra workspace per line),
pulls owners associated with workspace, and writes workspace + owners
to output TSV.
"""
#!/usr/bin/env python3
​
# Install dependencies with:
# pip install google-auth requests
#
# Set up Application Default Credentials with:
# gcloud auth application-default login
​
​
import argparse
import csv
import fileinput
import sys
​
import google.auth
from google.auth.transport.requests import AuthorizedSession
​
​
parser = argparse.ArgumentParser()
parser.add_argument("files", metavar="file", nargs="+", help="File containing workspaces (one per line, formatted namespace/name).")
args = parser.parse_args()
​
workspaces = list(line.strip() for line in fileinput.input(args.files))
​
credentials, _ = google.auth.default()
session = AuthorizedSession(credentials)
​
writer = csv.writer(sys.stdout, delimiter="\t")
​
for workspace in workspaces:
    owners = session.get(
        f"https://rawls.dsde-prod.broadinstitute.org/api/workspaces/{workspace}",
        params={"fields": "owners"}
    ).json().get("owners")
    
    writer.writerow([workspace, ",".join(owners)])
    