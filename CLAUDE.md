# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working in this repository.

## Reference registry (docs/references/)

This crate implements an external standard; the reference registry is part of the deliverable, not an afterthought. After any slice that cites a source (the standard, a paper, a vector set, an erratum, a reference implementation's behavior), append or update `docs/references/<slug>.md` with the global accretion schema: citation with edition, canonical URL or document number, Wayback archived URL saved at citation time, retrieval date, sha256 for binaries, license, vendor status (vendored, pointer only, legally cannot), and the code or test paths it grounds. Conformance vector sets always get entries, including their license and their coverage gaps; the gaps feed the README disclosure's named failure mode. Paywalled standards are pointer only: cite clause numbers and the free proxy literature (for 754: Goldberg, Muller et al; for 1788: the Pryce and Kulisch papers).
