# Deviations from the registration

ASCII, no em dashes.

None from the registered analysis: the data, graph, sets, statistics, thresholds, seeds, command and
single-run rule were followed as written in prereg/PREREG.md.

For the record, four things that the registration does not describe:
- 2026-11-14T16:03Z. The registered command was started through a wrapper, results/run_once.sh, which
  only writes the time stamps, the library versions and the file hashes into results/run.log.
- 2026-11-14T16:02Z. The metadata table was downloaded again for the run; its hash equals the one read
  at the names-only check.
- The Python patch version was 3.13.18 (the registration says 3.13).
- results/run_once.sh writes the exit code line after an echo, so the line exit code: 0 in
  results/run.log is the status of that echo, not of the run. Completion is shown by the line
  saved results_banc.json in the log and by the complete result file. run.log is left as written.
