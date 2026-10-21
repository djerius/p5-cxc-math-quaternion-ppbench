# Rotating one vector by many quaternions: stacking matrices and `splitdim`

## The statement

> Multiplying many quaternions by one vector works if you vertically
> stack the quaternions' rotation matrices, then do a `splitdim` on the
> output extremely tall vector.

In other words, to rotate one vector $`v`$ by $`N`$ quaternions:

1. Turn each quaternion into its 3×3 rotation matrix $`R_i`$.
2. Stack the $`N`$ matrices vertically into one $`3N \times 3`$ matrix.
3. Multiply that matrix by $`v`$, giving one $`3N`$-long column.
4. Split the column back into $`N`$ rotated vectors.

## Verdict

**Correct, but unnecessary in PDL, and slow.**

- **Correct.** A stacked matrix times $`v`$ is the $`N`$ products
  $`R_i v`$ laid end to end, so splitting the result gives the rotated
  vectors. In PDL each step does what the statement says, and the
  result matches the PP kernel `rotq` to about $`4 \times 10^{-16}`$.
- **Unnecessary.** PDL broadcasts the matrix product over extra
  dimensions, so the `(3,3,N)` matrices can be multiplied by the vector
  directly, with no stacking or splitting. Simpler still, `rotq`
  broadcasts one `(3)` vector against `(4,N)` quaternions, with no
  matrices at all.
- **Slow.** Even leaving out the cost of building the matrices, the
  stacked product is 19–68× slower than `rotq`, and 2.3–6.3× slower
  than broadcasting the same matrices.

## How it looks in PDL

The matrices are in PDL's usual layout: dims `(3, 3, N)`, with the
element in row $`r`$, column $`c`$ of matrix $`i`$ at `$m->at($c, $r, $i)`.

```perl
my $stack = $m->clump(1, 2);                     # (3, 3N): matrices stacked vertically, a view
my $tall  = $stack x $v->dummy(0);               # (3, 3N) x (1, 3) column  ->  (1, 3N)
my $r     = $tall->slice('(0)')->splitdim(0, 3); # (3N)  ->  (3, N)
```

- **Stacking:** `clump(1, 2)` merges the row and matrix dimensions, so
  the rows of all $`N`$ matrices follow each other. It is a view, not a
  copy.
- **The vector:** `$v->dummy(0)` turns the `(3)` vector into a `(1, 3)`
  column for `x`.
- **Splitting:** the product is a `(1, 3N)` column. `slice('(0)')`
  drops the unit dimension, and `splitdim(0, 3)` splits the $`3N`$
  values into `(3, N)`: one rotated vector per quaternion.

The two alternatives need no reshaping:

```perl
my $r = ( $m x $v->dummy(0) )->slice('(0)');     # broadcast: (3,3,N) x (1,3) -> (1,3,N) -> (3,N)
my $r = rotq( $q, $v );                          # PP kernel: (4,N) quaternions, (3) vector -> (3,N)
```

## Test ([splitdim.pl](splitdim.pl))

[splitdim.pl](splitdim.pl) builds `(3,3,N)` rotation matrices from
random unit quaternions. It does so with ordinary PDL operations, one
assignment per matrix element. It then:

1. checks the dims at each step, and the stacked result against `rotq`,
   with $`N = 5`$;
2. times four things, per quaternion, for one vector and $`N`$
   quaternions:
   - `rotq($q, $v)`;
   - building the matrices alone;
   - stacking, `x` and `splitdim`, on matrices already built;
   - `x` broadcasting over the matrices already built.

It was run 3 times on the same laptop as the other benchmarks
(Intel i7-8650U, perl 5.42.3, PDL 2.106). The table gives the per-cell
minimum: [results/splitdim-min.txt](results/splitdim-min.txt). The
load average was about 2.6 during the runs.

## Results

**Check** ($`N = 5`$): the stack is `(3, 15)` and a view. The product
is `(1, 15)`, and the split result is `(3, 5)`. It agrees with `rotq`,
as does broadcasting the matrices, to $`4.4 \times 10^{-16}`$.

**Time per quaternion**, rotating one vector:

|   N | `rotq` | building the matrices | stack + `x` + `splitdim` | broadcast `x` |
|----:|-------:|----------------------:|-------------------------:|--------------:|
| 1e2 | 28.5 ns | 1798.1 ns | 531.8 ns | 232.8 ns |
| 1e4 |  6.6 ns |   92.0 ns | 325.3 ns |  79.1 ns |
| 1e6 |  5.0 ns |  194.0 ns | 339.9 ns |  54.2 ns |

## Analysis

- **The stacked product is 19–68× slower than `rotq`** (325–532 ns
  against 5.0–28.5 ns per quaternion), before counting the matrices.
  It is also 2.3–6.3× slower than broadcasting `x` over the same
  matrices. Collapsing $`N`$ small products into one tall one does not
  help in PDL; the broadcast product is faster at every $`N`$ measured.
- **Building the matrices is an extra pass.** Matrices hold 9 numbers
  per rotation instead of 4, so they are 2.25× the data. Built with
  ordinary PDL operations, as here, they cost 92–1798 ns per quaternion,
  more than everything else. A PP kernel could build them much faster,
  but it would still be an extra pass over more data, and `rotq` needs
  none.
- **`rotq` is fastest** because it reads the 4 quaternion elements
  once, does the rotation in one pass, and needs no intermediate
  matrices. This matches the main benchmarks: rotating with the
  quaternion directly beats going through a matrix
  ([C2](README.md#c2)), and the PP kernel beats `x`
  ([C1](README.md#c1)).

## Pitfalls if the method is used anyway

- **Matrix layout.** `clump(1, 2)` stacks correctly only if the
  matrices are in PDL's `(column, row, N)` layout. With transposed
  matrices, `(row, column, N)`, it stacks columns instead of rows, and
  the result is silently wrong.
- **What `x` does.** Loading PDL::LinearAlgebra replaces PDL's `x` for
  the whole program ([C30](README.md#c30)), so the speed of this method
  depends on what else is loaded.
- **Only one vector.** Stacking handles many quaternions with *one*
  vector. With one vector per quaternion, the stacked matrix would have
  to be block-diagonal, $`3N \times 3N`$ and almost all zeros.
  Broadcasting, and `rotq`, handle that case directly.

## Conclusion

The statement is true, but the method shouldn't be used. To rotate one
vector by many quaternions in PDL, call `rotq($q, $v)`, or the API's
`quat_rotate_vec`, with the `(4, N)` quaternions and the `(3)` vector;
broadcasting does the rest. If the rotations are already matrices,
broadcast `x` over them rather than stacking them.

## Caveats

- These figures come from one machine, which was loaded during the runs.
  The ratios are more reliable than the absolute times.
- The matrices were built with ordinary PDL operations, which are slow.
  The comparison of the methods themselves (columns 2, 4 and 5)
  doesn't depend on that.
