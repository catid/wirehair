# Agent Instructions

This project uses **bd** (beads) for issue tracking. Run `bd prime` for full workflow context.

> **Architecture in one line:** Issues live in a local Dolt database
> (`.beads/dolt/`); cross-machine sync uses `bd dolt push/pull` (a
> git-compatible protocol), stored under `refs/dolt/data` on your git
> remote — separate from `refs/heads/*` where your code lives.
> `.beads/issues.jsonl` is a passive export, not the wire protocol.
>
> See [SYNC_CONCEPTS.md](https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md)
> for the one-screen overview and anti-patterns (don't treat JSONL as the
> source of truth; don't `bd import` during normal operation; don't
> reach for third-party Dolt hosting before trying the default).

## Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work atomically
bd close <id>         # Complete work
bd dolt push          # Push beads data to remote
```

## Non-Interactive Shell Commands

**ALWAYS use non-interactive flags** with file operations to avoid hanging on confirmation prompts.

Shell commands like `cp`, `mv`, and `rm` may be aliased to include `-i` (interactive) mode on some systems, causing the agent to hang indefinitely waiting for y/n input.

**Use these forms instead:**
```bash
# Force overwrite without prompting
cp -f source dest           # NOT: cp source dest
mv -f source dest           # NOT: mv source dest
rm -f file                  # NOT: rm file

# For recursive operations
rm -rf directory            # NOT: rm -r directory
cp -rf source dest          # NOT: cp -r source dest
```

**Other commands that may prompt:**
- `scp` - use `-o BatchMode=yes` for non-interactive
- `ssh` - use `-o BatchMode=yes` to fail instead of prompting
- `apt-get` - use `-y` flag
- `brew` - use `HOMEBREW_NO_AUTO_UPDATE=1` env var

## Bug Fix Passes

Bug fix passes are not just running unit tests. They require actually reading
the code thoroughly, especially newly edited areas, and looking for potential
bugs, edge cases, and incorrect assumptions before declaring the pass complete.
After making code edits, always do repeated bug fix passes until a full pass
finds no bugs.

## External Max-Effort Planning and Review

When the user authorizes Claude/Fable, use it as a read-only planning and
adversarial-review aid:

```bash
env -u ANTHROPIC_API_KEY -u OPENAI_API_KEY -u GOOGLE_API_KEY \
  -u GEMINI_API_KEY -u EXA_API_KEY -u SLACK_APP_TOKEN -u SLACK_BOT_TOKEN \
  claude -p --model fable --effort max \
  --permission-mode plan --tools 'Read,Grep,Glob,Bash' \
  --no-session-persistence --max-budget-usd 20 '<focused review prompt>'
```

- Prefer the signed-in Claude Max session (`ANTHROPIC_API_KEY` unset) when a
  separately configured API key is quota-limited.  Strip unrelated secrets
  from the subprocess environment.
- Scope the prompt to an exact commit/diff and subsystem.  State that the
  review is read-only, forbid edits and long workloads, enumerate the risks to
  inspect, and request prioritized findings with exact file:line evidence and
  minimal fixes.
- For complex work, use two narrow passes when worthwhile: before coding, ask
  for architecture alternatives, failure modes, and a validation plan; after
  the diff is coherent, ask for the adversarial file:line review.  Reusing one
  giant session is slower and makes budget failures less actionable.
- A budget-exhausted run that returns no report is not a review.  Narrow the
  scope or, when authorized, use a sufficient cap; the first successful
  repository review here needed more than a $5 cap and completed under $20.
- Fable is advisory, not an authority.  Reproduce every finding locally,
  reject unsound suggestions explicitly, apply fixes ourselves, and rerun the
  relevant gates.  Never expose private code to an external reviewer without
  the user's authorization.

## Web Research With Exa

Use the Exa MCP as a fast discovery-and-fetch pipeline when web research is
needed:

1. Call `web_search_exa` with a semantically rich description of the ideal
   page (including technology/version and the evidence sought), normally with
   3-5 results.  Do not submit a loose keyword pile.
2. Prefer primary papers, standards, official documentation, and upstream
   repositories from the result set.  Follow the best URLs with
   `web_fetch_exa`; batch related URLs and set a bounded character limit.
3. Verify material claims against the fetched primary source and cite that
   source, not the search highlight.  Use the normal browser/search path when
   Exa is unavailable, when an authoritative site needs direct inspection, or
   when current facts require an additional check.

Exa is for external discovery, not facts already available in the repository.
The validated pattern in this project was a targeted technical search that
surfaced the primary inactivation-decoding paper, followed by a focused fetch
of its authoritative arXiv page.

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:7510c1e2 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
