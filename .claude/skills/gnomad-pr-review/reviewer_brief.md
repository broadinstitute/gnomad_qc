# PR review: shared brief for reviewer agents

<!--
Coordinator: copy this file to {REVIEW_DIR}/brief.md and fill in every {PLACEHOLDER}
before launching agents. Delete sections that don't apply (e.g. "Tests" if the
repo has none). Agents read the filled-in copy, not this template.
-->

## The PR
- {OWNER/REPO} PR #{PR}, "{TITLE}", author {AUTHOR}.
- Repo: `{REPO_PATH}`. The working tree is checked out at the PR head, commit `{HEAD_SHA}`. Base is `origin/{BASE}`.
- See the full diff of a file with `git diff origin/{BASE}...HEAD -- <path>`. See the old version with `git show origin/{BASE}:<path>`.
- PR description: `gh pr view {PR}` (read-only). It lists the intended behavior changes. Anything that changes behavior but is NOT in the description is worth flagging.
- Project conventions: `{REPO_PATH}/CLAUDE.md`.

## Hard rules
- READ-ONLY. Do not edit, create, or delete files in any repo. Do not commit, push, or post anything to GitHub (no `gh pr comment`, `gh pr review`, `gh api -X POST/PATCH`).
- Do NOT run Dataproc, Hail Batch, hailctl, gcloud, or gsutil. Do not read from GCS.
- Never read All of Us or other individual-level data, in any form.
- Do NOT start Hail, Spark, or pytest. Several reviewers share one machine, and parallel JVMs exhaust its memory. The coordinator runs the tests and any repro scripts. If a bug can be proven with a small local Hail snippet that uses only hand-built Tables (`hl.Table.parallelize`), write the snippet in your report under "Repro" and the coordinator will run it.
- You may run plain Python that doesn't start Hail, plus git, grep, and similar read-only tools. The review venv is at `{VENV}` ({VENV_DESCRIPTION}).

## Code outside this repo that the PR depends on (read it, don't guess)
- gnomad_methods: `{GNOMAD_METHODS_PATH}` (check it is at the commit the PR needs; installed copy: `{VENV}/lib/python3.*/site-packages/gnomad/`).
- gnomad_qc: `{GNOMAD_QC_PATH}` (e.g. `resource_utils.py` for `PipelineResourceCollection`).
- Hail source, for checking semantics: `{VENV}/lib/python3.*/site-packages/hail/`.

When your code calls a gnomad_methods (or other external) function, open it and check what it really does with the arguments passed: argument order, defaults, how it handles missing values, which index it uses, what fields it returns.

## Cross-file call map
`{REVIEW_DIR}/callmap.txt` maps each function to its callers in other areas, including tests. Large files are split into lettered areas: {SPLIT_DESCRIPTION}.

## Reviewer roster (who owns what)
The coordinator sends you the agent IDs once everyone is running. Areas:
{ROSTER_AREAS}

## Talking to other reviewers
- Use the SendMessage tool. Your plain text is not visible to other agents.
- When your code depends on a function in someone else's area, read it yourself first. If the hand-off looks wrong (argument order, field names, schema, missing-value handling, units, index), send the owner one specific question with file:line, then keep working. Don't stop and wait for the answer.
- When you find a likely bug in someone else's area, send it to the owner. Also list it in your report under "Cross-area findings".
- When a peer asks you something, answer briefly and specifically with file:line, then go back to your work.
- Message `main` (the coordinator) only for something that blocks you.
- **Once a peer has delivered its final report, don't message it again.** A message to a finished agent restarts it, and its report then goes to you instead of the coordinator.

## What to look for, in priority order
1. **Wrong numbers.** Wrong math or sign, wrong filter, off-by-one, wrong array index (freq index, downsampling index, genetic ancestry group), wrong join key, double counting, wrong grouping, wrong chrX/chrY/PAR handling. Missing values in Hail: aggregators skip missing; `NA == x` is NA, so `filter()` drops the row; arithmetic with NA gives NA; `hl.or_else` and `hl.or_missing` change counts; a dict `.get()` returns missing silently while `dict[key]` errors.
2. **Pipeline breaks.** Crashes, wrong resource path, a flag wired wrong, a step reading the wrong step's output, a mismatch between the schema one step writes and the schema the next step reads, broken `--test`/`--overwrite`/version handling, a `--test` run that can write to a production path.
3. **Silent behavior changes compared with `{BASE}`** that the PR description doesn't mention, including changed CLI defaults.
4. **Cost on full data.** Extra shuffles, `.count()` used for logging, collecting big data to the driver, a missing checkpoint that makes Hail recompute a big pipeline, explode blowups, `order_by` on huge tables, no `naive_coalesce` after a heavy filter, hard-coded partition counts, `repartition()` on a large dataset.
5. **Tests** that don't test what they claim, assert a wrong expected value, or pin a bug as correct. Ask "would this test still pass if the code were broken?" Common holes: every group has one row, so a scaling step is never exercised; symmetric values that hide a sign or mirroring bug; the same value in every group, so a wrong index passes; a fixture with a field the pipeline never produces.
6. **Docs** that disagree with the code, or that describe code that isn't on this branch.
7. **Maintainability.** Dead code, duplicated logic, misleading names, docstrings, or comments, wrong type annotations, violations of the CLAUDE.md conventions.

Do NOT report: formatting that black/isort handle, pure style taste with no real effect, or anything already raised in the existing PR comments at `{REVIEW_DIR}/existing_comments.md` (unless you have something materially new to add; then say so).

## Priority levels
- **Blocker**: wrong numbers in a released metric on the default path, or the pipeline crashes or writes to the wrong place on a normal run.
- **High**: wrong results in a supported non-default path, wrong released output that users can't detect, overwrite or data-loss risk, a big cost problem, a CI break, or a test that pins a real bug as correct.
- **Medium**: edge-case bug, a misleading doc that could cause misuse, a missing test for risky logic, an undocumented behavior change.
- **Low**: cleanup, clarity, a small cost improvement.
- **Nit**: tiny.

## Confidence
- **Confirmed**: you traced the full code path, including the called functions in other files and repos.
- **Likely**: strong evidence, but one step is unverified (say which).
- **Question**: depends on the author's intent. Phrase it as a question.

Be honest. A short list of real problems beats a long list of guesses. Don't inflate priority.

## Line numbers
- Use line numbers in the working tree (the PR head). Name the function each line is in.
- Say whether the line is part of the PR diff. GitHub only accepts inline comments on lines in the diff. Check with `git diff origin/{BASE}...HEAD -U0 -- <path>` and read the `@@ -a,b +c,d @@` hunks: new lines c..c+d-1 are in the diff. If the line isn't in the diff, give the nearest line that is. Files not touched by the PR need a general PR comment.

## Report format (your final message)
Start with a one-line summary: how many findings at each priority.

Then one block per finding, sorted by priority:

```
### [Priority] Short title
- Where: path:LINE (function_name), in diff: yes/no (nearest diff line: N)
- Confidence: Confirmed / Likely / Question
- What: 1–3 plain sentences.
- Why it matters: 1 sentence (what goes wrong, for whom).
- Evidence: the path you traced (file:line → file:line). Quote at most 5 lines.
- Repro (optional): a small local Hail snippet using only hand-built Tables.
- Draft PR comment: 1–4 plain, direct sentences to the author, with a suggested fix. Written as a human reviewer. Common language, no jargon padding, no hedging filler.
```

After the findings:
- **Interface notes**: for each function in your area that another area calls, what it expects as input, what it guarantees as output (fields and types), and any assumption the caller must meet. Keep it short.
- **Cross-area findings**: suspected problems in other areas that you sent to the owner.
- **Open questions** for the author.
- **Checked and fine**: one line per risky area you checked that turned out OK, so the coordinator knows it was covered.
