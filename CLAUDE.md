# Project Instructions for AI Agents

This file provides instructions and context for AI coding agents working on this project.

## Web Research With Exa

Use the Exa MCP as the default discovery-and-fetch pipeline for web searches
and other external research:

1. Call `web_search_exa` with a semantically rich description of the ideal
   page (including technology/version and the evidence sought), normally with
   3-5 results. Do not submit a loose keyword pile.
2. Prefer primary papers, standards, official documentation, and upstream
   repositories from the result set. Follow the best URLs with
   `web_fetch_exa`; batch related URLs and set a bounded character limit.
3. Verify material claims against the fetched primary source and cite that
   source, not the search highlight. Use the normal browser/search path when
   Exa is unavailable, when an authoritative site needs direct inspection, or
   when current facts require an additional check.

Exa is for external discovery, not facts already available in the repository.
The validated pattern in this project was a targeted technical search that
surfaced the primary inactivation-decoding paper, followed by a focused fetch
of its authoritative arXiv page.

## Long-Running Claude CLI Reviews

Claude CLI reviews can legitimately run far longer than 300 seconds. Keep the
subprocess in a resumable session and poll it carefully; do not kill or
classify an otherwise-live review as failed merely because a five-minute
client/tool wait expired.

If the most recent Claude CLI attempt produced no report because the client
timed out, retry it once with the same scoped prompt and a sufficiently long
session/wall-time allowance. Do not treat authentication, quota, model, or
substantive command failures as timeout-only failures.

## Forwarded SSH Agent Refresh

The client SSH session cycles frequently, so an inherited `SSH_AUTH_SOCK` is
usually stale by the time a push is needed. Before Git/Beads SSH operations,
refresh the latest shell environment and select the newest live forwarded
`/tmp/ssh-*/agent.*` socket. Verify it non-interactively with `ssh-add -l`
and `ssh -o BatchMode=yes -T git@github.com`; do not treat the inherited socket
as authoritative.

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


## Build & Test

_Add your build and test commands here_

```bash
# Example:
# npm install
# npm test
```

## Architecture Overview

_Add a brief overview of your project architecture_

## Conventions & Patterns

_Add your project-specific conventions here_
