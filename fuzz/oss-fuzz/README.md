# OSS-Fuzz integration scaffolding

Files in this directory are staged for submission to
[`google/oss-fuzz`](https://github.com/google/oss-fuzz) under
`projects/pfloat/`. They are not used by pfloat's own CI; the
in-tree fuzz lane (`fuzz/fuzz_targets/`) runs under the project's
`Cargo.toml` per ADR-0013.

## Files

- **`Dockerfile`**: pulls pfloat from its primary repo and stages
  `build.sh` into `$SRC/`.
- **`build.sh`**: compiles each in-tree fuzz target with the
  libFuzzer sanitizer and copies the resulting binaries into
  `$OUT/`. Seven targets at slice 6g close: `parse`, `arith`,
  `exp_log_family`, `trig`, `hyperbolic`, `specials`, `fmt`.
- **`project.yaml`**: project metadata. Sets language to `rust`,
  primary_contact to Parnell's address, sanitizer to `address`.

## Submission

The OSS-Fuzz onboarding flow:

1. Fork `google/oss-fuzz`.
2. `git checkout -b add-pfloat`.
3. Create `projects/pfloat/` containing copies of
   `Dockerfile`, `build.sh`, and `project.yaml` from this
   directory.
4. Open a PR with the standard "Add new project" title.
5. Address any review feedback; merge.

Per the `feedback_no_prs_solo.md` engineering memory, pfloat
itself does not use the GitHub PR flow for its own merges; the
OSS-Fuzz upstream PR is the one cross-repo exception. ADR-0013
names the upstream PR as the closer for slice 6g.

## Local verification

The in-tree fuzz lane runs the same targets without the
OSS-Fuzz wrapper:

```sh
cargo +nightly fuzz run <target> -- -max_total_time=60
```

See `fuzz/Cargo.toml` for the target list.
