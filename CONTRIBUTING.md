# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

Keep in mind as you contribute, that code, docs and other material submitted to this projects are considered licensed under MIT license.

## Setup developement environment

We recommand to install rust with [rustup](https://rustup.rs/).
If you want perform a documentation contribution install [mdbook](https://rust-lang.github.io/mdBook/guide/installation.html).

## Contribution

Before start any modification please create a specific branch:
```bash
git switch -c fix_11         # branch create to fix issue 11
git switch -c feat_index_rc  # branch to add a new index reverse complement method
```

## Code contribution

Before submit pull request make sure you run:

```bash
cargo fmt --all
cargo clippy --all-targets
cargo test --all-targets
```

You can check your new code are covered by run:
```bash
cargo tarpaulin
```
And open `target/coverage/tarpaulin-report.html`

### Documentation pull request

After change you can run:
```
cargo doc
```
And open `target/doc/kff/index.html` to check effect of your change.

### Website pull request

After change you can run:
```
mdbook serve
```
To check effect of your change.
