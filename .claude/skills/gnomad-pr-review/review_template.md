# PR #{PR} review: ranked findings with draft comments

- Line numbers are at the PR head (`{HEAD_SHA_SHORT}`), so they match the "Files changed" tab on GitHub.
- Models: coordinator {SESSION_MODEL}; reviewers {e.g. "Opus: core logic A/B/C, pipeline; Sonnet: resources, tests, docs/config"}; second opinion on High findings: {"yes, other model" / "no"}.
- "Not in diff" means GitHub won't take an inline comment on that line. Use the line given, or leave a general PR comment.
- Tags:
  - **[ran]**: reproduced locally on Hail {HAIL_VERSION}.
  - **[code]**: confirmed by reading the code path end to end.
  - **[Q]**: depends on what the author intended, so it's phrased as a question.
- Already raised by others, so not repeated here: {EXISTING_COMMENT_TOPICS}.

---

## HIGH

### H1. {One-line claim: what is wrong and what it causes} [ran|code] [Q?]
- **Where:** `{path}:{line}` (`{function}`)
  - Related lines, if any: `:{line}` (`{function}`), with one phrase each on what happens there.
- **What:**
  - 2–5 bullets: the mechanism, and what `{BASE}` did instead.
  - **Verified end to end** ({the real functions chained, and the input used}): {what came out}.
  - **Second opinion** (only if run): {model} {agrees / disagrees}: {one-line reason}.
  - **Correction** (only if verification narrowed the claim): {what the earlier version overstated}.
  - **Severity note** (if it's a non-default path or depends on intent): {why it's still ranked here}.
- **Draft comment:**
  > {1–5 plain sentences to the author: what's wrong, the evidence in one clause, a suggested fix, and a test that would pin it.}

---

## MEDIUM

### M1. {One-line claim} [ran|code] [Q?]
- **Where:** `{path}:{line}` (`{function}`)
- **What:** {1–4 bullets}
- **Draft comment:**
  > {plain sentences}

---

## LOW

<!-- Group by area. One bullet per finding, with the draft inline. -->
### {Area}
- **L1. `{path}:{line}`: {one-line claim}.** {Optional extra sentence.}
  - Draft: "{short draft comment}"

---

## NIT
- `{path}:{line}`: {one line}.

---

## Open questions to ask {AUTHOR} (general comment)
1. {question}

## Checked and fine
- **{Area}:** {what was checked and found correct}.
- **Checks run:** {tests: N passed / skipped / xfail on Hail X with Y GB heap}; {lint commands}.
