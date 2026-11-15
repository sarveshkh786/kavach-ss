#!/usr/bin/env bash
# Compile the field arithmetic in isolation and count conditional branches
# in every secret-handling primitive. Expected: 0 everywhere except the fixed
# 64-iteration loop counter in pow.
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd); tmp=$(mktemp -d)
cat > "$tmp/a.rs" <<RS
#![crate_type="lib"]
#[path="$here/src/field.rs"] mod field;
use field::Fp;
#[no_mangle] pub fn k_mul(a:Fp,b:Fp)->Fp{a.mul(b)}
#[no_mangle] pub fn k_add(a:Fp,b:Fp)->Fp{a.add(b)}
#[no_mangle] pub fn k_sub(a:Fp,b:Fp)->Fp{a.sub(b)}
#[no_mangle] pub fn k_ct_eq(a:Fp,b:Fp)->u64{a.ct_eq(b)}
#[no_mangle] pub fn k_reduce(a:u64)->Fp{Fp::new(a)}
#[no_mangle] pub fn k_pow(a:Fp,e:u64)->Fp{a.pow(e)}
RS
(cd "$tmp" && rustc --edition 2021 -O --emit asm a.rs 2>/dev/null)
for f in k_mul k_add k_sub k_ct_eq k_reduce k_pow; do
  body=$(awk "/^$f:/,/\.cfi_endproc/" "$tmp/a.s")
  printf "%-10s instructions=%-3s conditional_branches=%s\n" "$f" \
    "$(echo "$body" | grep -cE '^\s+[a-z]')" "$(echo "$body" | grep -E '^\s+j[a-z]+' | grep -vc jmp || true)"
done
rm -rf "$tmp"
