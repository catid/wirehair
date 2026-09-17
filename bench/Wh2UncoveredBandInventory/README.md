# Uncovered-band recovery inventory — disabled draft

This directory checkpoints unfinished work on `main`. It is not qualified,
not connected to the production build, and provides no recovery or speed result.
The controller and native inventory entry points fail closed before loading a
codec, reading qualification evidence, or creating output. There is no command-line
override. Before removing the guard, resolve the review findings and pass the
prerequisite synthetic tests and source review tracked in Beads issue
`wirehair-sxvz.16.1.20.82.3.2`. Actual neutral qualification must then pass before
any scientific launch.

The intended diagnostic compares unchanged ordinary WH2 and owned WH1 at
K7/9/12/16, widths 2/64/1280, using actual decoder outcomes and independent
GF(256) packet/rank checks. The prospective roster has 3,072 loss cases and
24 separate low-repair cases. It does not evaluate a candidate, timing,
holdout performance, or all-K success. Widths share roots and are not independent
seed samples.

The sole scientific namespace `/var/tmp/wh2-uncovered-band-inventory-r0` has
not been launched. Native/portable neutral qualification has not been run.
Historical codec paths and hashes in the draft are proposed inputs, not evidence
that this harness is qualified. Its source manifest intentionally still names
support/launcher test files that have not yet been implemented.

The initial review found gaps in neutral-success enforcement, producing-source
binding, backend parity, recursive type equality, and neutral resource bounds.
Those findings and remaining qualification work are recorded in the Beads issue;
this merge does not resolve or waive them.

Synthetic-only checks (no native codec calls):

```bash
python3 -B -m unittest discover -s bench/Wh2UncoveredBandInventory -p 'test_*.py'
```

Passing these checks is not permission to launch or evidence of codec performance.
