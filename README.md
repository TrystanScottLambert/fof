Rust crate for handling the friends-of-friends algorithm functionality for finding galaxy-groups in redshift surveys.

Tests can be run using:

```rust
cargo test
```

A small benchmarking test is included which can be run as
```rust
cargo test --release -- --nocapture tests/integration_bench.rs
```
