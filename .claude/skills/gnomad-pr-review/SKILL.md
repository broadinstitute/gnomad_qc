---
name: gnomad-pr-review
description: Extremely thorough, multi-agent review of a large gnomAD Hail pipeline PR (gnomad_qc, gnomad-constraint, gnomad_methods). Splits the PR by file across reviewer agents that message each other about cross-file calls, runs the tests and local Hail repros, re-checks every High finding end to end, and writes a ranked REVIEW.md with file:line, function, and an editable draft PR comment per finding. Read-only - never posts to GitHub, never submits to Dataproc or Hail Batch.
argument-hint: "[PR number, URL, or branch] [model=mixed|opus|sonnet] [second-opinion=yes|no]"
disable-model-invocation: true
---

# gnomAD PR review

Review a pull request the way a careful senior reviewer would, using one agent per file or area, and hand the user a ranked list of findings they can edit and paste into GitHub.

The target is the PR the user named in the invocation (number, URL, or branch). If they didn't name one, use the PR for the current branch (`gh pr view`). The PR can be in a sibling repo (e.g. `../gnomad-constraint`); work in that repo's directory.

The invocation can also carry options as `key=value` words after the PR:
- `model=mixed|opus|sonnet`: which model the reviewer agents use (see "Model choice").
- `second-opinion=yes|no`: whether to get an independent check on each High finding from the other model (Phase 6).

Any option given here answers the matching Phase 1 question, so don't ask it again. If a value isn't recognized, ask.

## How the user wants this done

- **Plan first.** Size up the PR, then tell the user your plan and wait for a go-ahead before starting the review. If the plan changes along the way, tell them.
- **Read-only.** Never post comments or reviews, never commit or push, never submit to Dataproc or Hail Batch, never run `hailctl`, and never write to GCS. Reading GCS (even small reads) also needs the user's OK first. All of Us or other individual-level data must never be read locally, ever.
- **One agent per file, and the agents talk.** When one file uses another file's functions, the owning agents check the hand-off with each other.
- **Ranked output with locations.** Every finding gets a priority, a file:line (at the PR head), the function name, and a draft PR comment.
- **Draft comments are for the user to edit**, not final text. Write them in plain, direct language, as a human reviewer would, with no filler.
- **Be direct and concise, in common language,** in chat and in the file.
- **Check the High findings end to end** before calling the review done (Phase 6). The user expects this every time.
- When the user asks to see the review file, **print its contents in your reply.** Tool output isn't shown to them.

## Model choice

Set the model for each reviewer agent with the Agent tool's `model` parameter. Without it, agents inherit the session's model.

| Choice | Reviewer agents | When to use |
|---|---|---|
| **mixed** (default, recommended) | **Opus** for core-logic areas (code that computes the numbers) and pipeline wiring (the main CLI script). **Sonnet** for tests, docs/config, and resources/constants. | Most PRs. Opus goes where the hard cross-function tracing is; Sonnet takes the text-vs-code checks. |
| **opus** | Opus for every area. | Very large or high-stakes PRs that change released numbers. Most thorough, most expensive. |
| **sonnet** | Sonnet for every area. | Smaller PRs or a quick first pass. Cheaper and faster. |

- **Reference point:** with all-Opus on a ~9k-line PR (8 agents), each agent used about 210–270k tokens and took 12–19 minutes. Check current per-token pricing before quoting cost; don't quote prices from memory.
- **The coordinator (you) stays on the session's model.** Don't pin a model in this skill. Checking findings, repros and the write-up happen in the main conversation and are the hardest part. In the plan, tell the user that running the session on Opus gives the best results.
- **Don't use Haiku for reviewers;** it's too shallow for this work. The Agent tool also accepts `fable`. It hasn't been compared with Opus for this skill, so use it only if the user asks for it.
- **If an area mixes core logic with simple code,** give it the stronger model.

## Phase 0: Size up the PR (before the plan)

1. Get the PR facts. Run:
   ```bash
   gh pr view <PR> --json number,title,author,headRefName,baseRefName,additions,deletions,changedFiles,url
   gh pr view <PR> --json body -q .body
   ```
2. Collect what others have already said, so agents don't repeat it. Save it to `{REVIEW_DIR}/existing_comments.md`, where `REVIEW_DIR=/tmp/pr<N>_review`:
   ```bash
   gh api repos/<owner>/<repo>/pulls/<PR>/comments --paginate -q '.[] | "=== \(.path):\(.line // .original_line) [\(.user.login)]\n\(.body)\n"'
   gh pr view <PR> --json comments -q '.comments[] | "[\(.author.login)] \(.body)"'
   ```
   Read the author's own comments closely. They often document known problems, like a CI failure with a fix that was never committed.
3. Size it up:
   - `git fetch`, then run `git diff --stat origin/<base>...origin/<head>`.
   - Get the line count of each changed file at the PR head.
   - For big files, list the functions with `grep -n "^def "`.
4. Check whether the local branch matches the PR head. Line numbers must match GitHub's.
5. Check the local environment: the Python and Hail version, and the gnomad_methods version and commit. Compare them with what the PR needs (CI config, requirements, any functions it imports from gnomad_methods).

## Phase 1: Present the plan and ask

Tell the user:
- The PR, its size, and where the code is (source vs. tests vs. docs).
- The agent split: one per file. Split any file over about 1,000 lines into stages at function boundaries. Group small related files together. Give the tests one or two agents, and docs plus config one agent. Keep it under 10 agents.
- How the agents will talk to each other (call map, roster, peer messages, then a second round run by you).
- How findings get checked (Phases 5 and 6), and the output format.
- The ground rules: read-only, no Dataproc or Batch.
- The cost: many agents, heavy on tokens, and it takes a while.

Then ask, giving your recommended answer for each. Skip any question the invocation options already answered.
1. **Checkout.** Fast-forward the local branch to the PR head, or use a separate worktree?
2. **Tests.** Static review only, or set up a separate venv and run the tests and local repros (local Spark, hand-built Tables only)?
3. **Scope.** Everything, including tests, docs and config, or source code only?
4. **Reviewer model.** Mixed (recommended), all Opus, or all Sonnet? In the plan, show which areas get which model under mixed.
5. **Second opinion on High findings.** Yes or no? Yes means one extra agent per High finding, on the other model. Recommend yes for PRs that change released numbers.

AskUserQuestion takes at most 4 questions per call. If all 5 are still open, put the second-opinion question in a short follow-up call.

Wait for answers.

## Phase 2: Setup

1. **Check out the PR head** the way the user chose. Record the short SHA; every line number refers to it.
2. **Test environment**, if approved:
   - Create the venv in `{REVIEW_DIR}/venv` with `uv venv --python <CI's version>`.
   - **Mirror CI's install order exactly.** Installing the requirements and Hail in one resolve can quietly pick an older Hail: a `dill==0.3.7` pin resolved to Hail 0.2.135, which has no `hl.qgamma`, so the gamma tests skipped without failing.
   - Print `hl.version()` and check for any feature the tests need.
   - Check whether gnomad_methods excludes the Hail version CI installs (e.g. `hail!=0.2.139` in its `requirements.txt`). That's a finding by itself.
3. **Run the tests**, if the repo has any (gnomad_qc itself has no unit tests; gnomad-constraint does). Run them in the background, and force the local backend:
   ```bash
   HAIL_QUERY_BACKEND=spark \
   PYSPARK_SUBMIT_ARGS="--conf spark.driver.bindAddress=127.0.0.1 --conf spark.driver.host=127.0.0.1 --driver-memory 4g pyspark-shell" \
   python -m pytest -q -rfEsxX --durations=10 -p no:cacheprovider
   ```
   Report passed, failed, skipped, xfail and xpass counts. Zero xfail means any bug found is either untested or locked in by a test.
4. **Run CI's lint commands**, copied from `.github/workflows/ci.yml`. Also run pylint once with warnings enabled, to spot dead code.
5. **Build the call map.** This skill's `callmap.py` does it; run it from the PR repo root:
   ```bash
   python <skill_dir>/callmap.py <src files...> --split <big_file>:<line>,<line> > {REVIEW_DIR}/callmap.txt
   ```
   It lists every function used across areas, with callers as `file:line(function)`, including tests.
6. **Write the brief.** Copy `reviewer_brief.md` from this skill to `{REVIEW_DIR}/brief.md` and fill in every placeholder.

## Phase 3: Reviewer agents

- **Launch all reviewers in one message**, so they run in parallel. Use `subagent_type: general-purpose`, which has SendMessage.
- **Set `model` on every Agent call to match the user's choice** (`opus` or `sonnet`; see "Model choice"). Add each agent's model to `roster.md`.
- **Each prompt contains:**
  - the agent's role name, and "read `{REVIEW_DIR}/brief.md` and follow it exactly";
  - its area, with exact line ranges;
  - its callers and callees from the call map, with the owning area of each;
  - a focus list written for that area.
- **Focus lists should name concrete risks**, not generic ones:
  - preprocessing: missing-value handling, freq/meta index alignment, chrX/chrY/PAR, filters;
  - model application: coefficient order, which index each model uses, equivalence between alternative paths;
  - metrics: CI parameterization, z-score sign, which transcripts are ranked;
  - pipeline: every call's signature, step-to-step resource wiring, `--test`/`--overwrite` safety, changed defaults;
  - resources: path collisions, test paths, stale checkpoints;
  - tests: "would it still pass if the code were broken?";
  - docs: check against the code, not against other docs.
- **Ask each agent to open the external functions it calls** (gnomad_methods at the right commit) rather than assume what they do.
- **The Agent tool has no name parameter.** After launch, write `{REVIEW_DIR}/roster.md` with each area and its agent ID, then send every agent a message pointing to it, along with the test results once you have them.
- **Agents must not run Hail, Spark or pytest.** On a 16 GB laptop, 8 JVMs thrash. Agents write repro snippets using hand-built Tables, and you run them.

## Phase 4: Collect

- As each report arrives, add a condensed entry to `{REVIEW_DIR}/ledger.md`: priority, file:line, one-line claim, confidence, and whether a repro is needed.
- Merge duplicates. A finding reached independently by two agents is stronger; note that.
- **Watch for duplicate or misrouted reports.** A message sent to a finished agent restarts it, and its report may go to the sender instead of you. Agents sometimes tell you about this; check the ledger instead of resending.
- Run a second round only for cross-area questions that are still open, by messaging the owners. Most get settled agent to agent.

## Phase 5: Check every High and Medium yourself

Don't trust agent output by itself; mistakes come back sounding just as confident as real findings.

- Re-read the code at every cited line. Confirm line numbers with `grep -n`, and confirm in-diff status from the `git diff -U0` hunks.
- **Run repros one JVM at a time.** Batch several into one script so Hail starts once. Print `RESULT ...` lines and grep for them.
- For CI claims, read the logs with `gh run view <id> --log-failed` (read-only). Check that the run's branch has the same tests and config as the PR.
- **For "this test can't catch X" claims, break the code on purpose.** Monkeypatch the function with the logic removed (`inspect.getsource` → edit → `exec`), and show that the test still passes.
- Drop or downgrade anything that doesn't hold. Fix wrong line numbers.

## Phase 6: Check every High finding end to end (always)

The user asks for this after every review, so do it before presenting results. For each High finding:

1. **Look for reasons it's wrong.**
   - Does other code (upstream, downstream, or in gnomad_methods) already handle it?
   - What did the base branch do? Is this actually new?
   - Does the draft comment overstate the problem?

   Example from a past review: a finding said zero-observed transcripts newly get ranked. Main ranked them too; what was actually new was the missing flag and their inclusion in published thresholds. The draft had to be corrected.
2. **Reproduce it through the real public functions, chained the way the pipeline chains them**, not only the helper where the bug lives.
   - Reuse test fixtures by importing the test module (add `tests/...` to `sys.path`).
   - Monkeypatch fixture constants to add the case you need.
   - Checkpoint intermediate Tables to `{REVIEW_DIR}/repros/` so a failed later step doesn't redo everything.
3. **Test CI and runtime claims with controlled comparisons**, e.g. the same suite at two heap sizes, or on two Hail versions.
4. **Second opinion (only if the user said yes).** For each High finding, start one fresh agent on the *other* model from the one that found it: Sonnet if an Opus agent found it, Opus if a Sonnet agent did.
   - Launch them all in one message, with `subagent_type: general-purpose`.
   - **Give it the code, not the conclusion.** Include the repo path, the PR head SHA, the file:line and function, and a neutral question about behavior. Don't include the finding, the draft comment, or what you expect. For example:
     > In `get_exomes_observed_and_possible` (`constraint.py:434-541`), what do `observed_variants` and `possible_variants` come out as for an exome variant that has a non-empty filter set and AF 0.2, with `max_af=0.001`? Compare that with an unfiltered variant at the same AF. Trace through any gnomad_methods functions it calls, and answer with file:line evidence.
   - Point it to the "Hard rules" section of `{REVIEW_DIR}/brief.md`: read-only, no Hail, no GCS.
   - **If it disagrees, re-check the code yourself.** Report the disagreement to the user rather than quietly choosing a side.
5. **Update `REVIEW.md` under each High finding:**
   - a "Verified end to end" bullet with the input and the output;
   - a "Second opinion" bullet, if one was run: the model, and whether it agrees, with its one-line reason;
   - a "Correction" bullet if the claim got narrower, along with a rewritten draft comment;
   - a "Severity note" if it's a non-default path or depends on intent.
6. **Say what you couldn't check**, e.g. how many real transcripts are affected, which would need GCS reads. Offer to do it, and ask first.

## Phase 7: Write up

- Write `{REVIEW_DIR}/REVIEW.md` using this skill's `review_template.md`, with sections High, Medium, Low (grouped by area), Nit, open questions for the author, and checked-and-fine.
- Every finding gets a location, the function name, a short "what", and a draft comment.
- Mark lines that aren't in the diff, and files the PR doesn't touch (those need a general comment).
- Fill in the "Models" line in the header: the coordinator model, which reviewers ran on which model, and whether a second opinion was run.
- **Chat summary:**
  - one opening line with the counts and the file path;
  - what was run, including the models used;
  - the High findings with full draft comments;
  - the Medium findings as a table (location, issue);
  - the Low and Nit findings as counts plus the main themes;
  - decisions that belong to the author;
  - a cleanup note for the venvs in `/tmp`.

## Gotchas

- Foreground `sleep` is blocked. Run long jobs with `run_in_background` and wait for the completion notice.
- Hail `Table.distinct()` needs a keyed Table. Use `group_by(...).aggregate(...)` to build small lookup Tables in repros.
- The test suite can pass locally at a 1 GB heap and still run out of memory on a GitHub runner. Test heap claims at 512 MB.
- `gh` comment line numbers (`.line // .original_line`) can point at an older commit than the PR head. Re-check them before calling something a duplicate.
- The base branch can be far behind. Compare behavior against `origin/<base>` and also against the last released state, if the PR says releases were cut from the branch.
