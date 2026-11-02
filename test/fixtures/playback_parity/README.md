# Playback parity corpus

`loopback_demo.json` is a verbatim copy of strumento's parity corpus
(`devices/_fixtures/parity/loopback_demo.json`, schema `strumento.parity-corpus/1`, qick 0.2.367),
sha256 `f3fbafc9061aca41f333658ab72700f42748a1e02dfdd17b5b7f40e89512d241`, from
harmoniqs/strumento-dev commit `42be0977569733ea297c333af580c8e17cf4e635` (PR #244, merged as
`35c5267186a1a8444c1dc06da6590eaecbf7c99b`). Do not edit it; refresh by copying a newer corpus and
updating this line. `test/test_playback_parity.jl` plays every pulse from the corpus spec and
asserts exact int16 code and gain-code equality (recipe: the `Playback` module docstring).
