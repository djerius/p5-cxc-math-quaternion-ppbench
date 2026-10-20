# Rotating vectors with PDL: PDL::PP vs. PDL::Transform vs. plain PDL

Claude did this.

These benchmarks were run to help design `CXC::Math::Quaternion`. The
question was how the vectorised quaternion code should be built:

- as PDL::PP kernels,
- on top of `PDL::Transform` (`t_linear`),
- or as compositions of ordinary PDL operations.

A second question came up while interpreting the results: why does the
PP kernel get about 3–4× slower per vector once $`N`$ passes roughly
$`10^6`$?

## Building and running

### Requirements

- Perl with PDL. The results here used perl 5.42.3 and PDL 2.106. The
  build uses PDL::PP and `PDL::Core::Dev`, which come with PDL.
- A C compiler.
- Linux with glibc for the allocation experiments:
  - [faults.pl](faults.pl) reads `/proc/self/stat`.
  - The `MALLOC_*` environment variables below are glibc-specific.

  The other scripts run anywhere PDL does.

### Build

```sh
perl Makefile.PL
make
```

This compiles the PP module `QBench` (`rotq`, `rotm`, `rotq3`) from
[qbench.pd](qbench.pd) into `blib/`.
- **No `make install` is needed.** Each script finds the module in the
  `blib/` next to it, using `FindBin`, so the scripts can be run from
  any directory.
- **After changing [qbench.pd](qbench.pd), re-run
  `perl Makefile.PL`, not just `make`.** The generated Makefile lists
  one `pp-<name>.c` file per kernel, so a new kernel is otherwise
  silently left out.
- MakeMaker prints `Setting VERSION via file 'qbench.pd' failed`. This
  is harmless.

### Run

Each script prints its own tables to stdout. All timings are in ns
(per vector, or per call where labelled).

```sh
perl bench.pl     # PP vs. matrix vs. x vs. t_linear vs. PDL ops   (~35 s)
perl alloc.pl     # fresh vs. preallocated output, threading       (~20 s)
perl mmap.pl      # fresh-output timing alone                      (~5 s)
perl faults.pl    # page faults per call                           (~10 s)
perl cols.pl      # three-column input vs. (3,N)                   (~40 s)
perl nfs.pl       # zeroes vs. new_from_specification outputs      (~35 s)
perl ovh.pl       # fixed per-call costs at N = 1 and 1000         (~20 s)
perl ovh2.pl      # ndarray construction vs. destruction           (~2 s)
```

To see the effect of glibc's mmap threshold, run [mmap.pl](mmap.pl),
[faults.pl](faults.pl), [cols.pl](cols.pl) or [nfs.pl](nfs.pl) a
second time with large blocks kept on the reusable heap:

```sh
MALLOC_MMAP_THRESHOLD_=4294967296 MALLOC_TRIM_THRESHOLD_=4294967296 perl mmap.pl
```

The column, allocation and overhead tables in this document are the
**per-cell minimum over 3 runs**. For `cols.pl` and `nfs.pl` there are
3 runs with each malloc setting, alternating between the settings. The
raw outputs and the minima are in [results/](results/), and the minima
were produced with [results/min.pl](results/min.pl):

```sh
for i in 1 2 3; do
  perl cols.pl > results/cols-default-$i.txt
  MALLOC_MMAP_THRESHOLD_=4294967296 MALLOC_TRIM_THRESHOLD_=4294967296 \
    perl cols.pl > results/cols-nommap-$i.txt
done
perl results/min.pl results/cols-default-?.txt > results/cols-default-min.txt
perl results/min.pl results/cols-nommap-?.txt  > results/cols-nommap-min.txt
```

`results/ovh-min.txt` came from a similar three-run loop over
[ovh.pl](ovh.pl), with the minimum taken separately. `min.pl` only
recognises rows that start with an `N` such as `1e+03`, so it does not
handle `ovh.pl`'s rows.

Timings on a laptop or a shared machine vary by 10–30% from run to
run, so compare ratios rather than absolute numbers.

## Hypotheses

1. **H1: a dedicated PP kernel beats the generic alternatives.** A PP
   kernel does the whole rotation in one pass over the data with no
   temporaries. It should beat:
   - plain PDL operations, which allocate a temporary for every
     intermediate;
   - `x` (matrix multiply);
   - `PDL::Transform::t_linear`, which adds Perl-level dispatch and
     generic machinery.
2. **H2: rotating with the quaternion directly is at least as fast as
   building a 3×3 matrix and using that.** The quaternion form,
   $`v' = v + 2w(u \times v) + 2u \times (u \times v)`$ with
   $`u = (x,y,z)`$, needs about the same number of flops as a
   matrix–vector product. It also avoids the matrix construction and
   reads only 4 numbers instead of 9.
3. **H3 (added after the first run): the slowdown at large $`N`$ comes
   from the kernel or from memory bandwidth.**
   - **H3a:** it is from memory allocation, not from the kernel. glibc
     `malloc` serves blocks above the dynamic mmap threshold with fresh
     `mmap` calls. Its default ceiling is 32 MB on 64-bit systems. A
     fresh mapping is page-faulted in one 4 KiB page at a time when the
     kernel first writes to it, and it is unmapped again when the
     output piddle is freed. So every call pays that cost again.
   - **H3b (alternative):** it is the memory bandwidth limit, reached
     once the working set falls out of the 8 MB L3 cache.

## Setup

| Item         | Value |
|--------------|-------|
| CPU          | Intel Core i7-8650U, 4 cores / 8 threads, 1.9 GHz base, laptop (frequency scaling active) |
| Caches       | L1d 32 KiB/core, L2 256 KiB/core, L3 8 MiB shared |
| Memory       | 15 GiB |
| OS           | Linux 6.12 (Debian kernel), benchmarks run inside a yolobox container (Ubuntu 24.04 userland, glibc 2.39) |
| THP          | `transparent_hugepage=always` |
| Perl         | plenv perl 5.42.3, `usemymalloc=n` (system glibc malloc) |
| PDL          | 2.106, autopthread target 8 threads, size threshold 1 M elements |
| Compiler     | gcc 13.3.0, `-O2` |

### Code

- [qbench.pd](qbench.pd) is a minimal PP module with two kernels, both
  `GenericTypes => ['D']`:
  - `rotq(q(n=4); v(m=3); [o] r(m=3))` rotates using a scalar-last
    quaternion $`(x,y,z,w)`$ and Hamilton convention, $`v' = q v q^*`$.
    It uses the cross-product form above (≈ 18 multiplies, 12 adds).
  - `rotm(mat(i=3,j=3); v(i); [o] r(j))` is a plain 3×3 matrix–vector
    product.
- [Makefile.PL](Makefile.PL) builds it with `pdlpp_stdargs`.
- [bench.pl](bench.pl) is the main comparison. It runs five
  implementations of the same rotation:

  | Label       | Implementation |
  |-------------|----------------|
  | `rotq(PP)`  | `rotq($q, $v)` |
  | `rotm(PP)`  | `rotm($R, $v)`, with $`R`$ precomputed from $`q`$ |
  | `matmult`   | `$v x $R->transpose` (PDL's built-in matrix multiply) |
  | `t_linear`  | `t_linear(matrix => $R->transpose)->apply($v)` |
  | `PDL ops`   | `rotq_ops`: the same formula written with `crossp`, `slice`, `dummy`, `*`, `+` |

  The first table rotates $`N`$ vectors by one quaternion. The second
  rotates $`N`$ vectors by $`N`$ different quaternions, one each, using
  broadcasting.
- [alloc.pl](alloc.pl) compares three runs of `rotq`:
  - a freshly allocated output (`new out`);
  - a preallocated output passed as the `[o]` argument (`prealloc`);
  - a preallocated output with autopthreading forced to 1 thread.
- [mmap.pl](mmap.pl) times `new out` alone. It was run with default
  glibc settings and then with `MALLOC_MMAP_THRESHOLD_` /
  `MALLOC_TRIM_THRESHOLD_` raised to 4 GiB, so that large blocks come
  from the reusable heap instead of fresh mmaps.
- [faults.pl](faults.pl) counts minor page faults per `rotq` call by
  reading `/proc/self/stat` field 10. It compares the count with the
  number of 4 KiB pages in the output.
- [cols.pl](cols.pl) benchmarks vectors that arrive as three separate
  `(N)` columns. It uses a third kernel added to
  [qbench.pd](qbench.pd), `rotq3(q(n=4); a(); b(); c(); [o] ra();
  [o] rb(); [o] rc())`, which has the same arithmetic as `rotq` but
  reads and writes separate columns. See
  [Column-layout input](#column-layout-input).
- [nfs.pl](nfs.pl) compares allocating an output with `zeroes` against
  `PDL->new_from_specification` (plus `make_physical`), and measures
  passing either one, freshly allocated on every call, to
  `rotq`/`rotq3`. See
  [Allocating outputs without zero-filling](#allocating-outputs-without-zero-filling).
- [ovh.pl](ovh.pl) and [ovh2.pl](ovh2.pl) measure fixed per-call
  costs at $`N = 1`$: creating ndarrays from Perl, versus outputs that
  the PP wrapper creates in C.

### Method

- Each measurement does one warm-up call. It then calls the code
  repeatedly for at least 0.5 s ([bench.pl](bench.pl)) or 1 s (the
  others) and reports the mean wall-clock time per call, divided by
  $`N`$. All numbers below are **ns per vector**.
- Input vectors are `random(3, N) - 0.5` and the quaternion is a fixed
  normalised $`(0.1,-0.3,0.5,0.8)`$.
- Correctness check: all five implementations agree with `rotq` to
  ≤ 2.2e-16 (max abs difference over 1000 vectors).
- This is a laptop inside a container with frequency scaling. Repeated
  runs differ by about 10–20% in the memory-bound region. For example,
  `rotq` at $`N=10^6`$ measured 4.6, 5.3 and 5.7 ns in three runs. Only
  differences of about 1.5× or more should be read as meaningful.

## Results

### One quaternion, $`N`$ vectors ([bench.pl](bench.pl))

|     N | rotq(PP) | rotm(PP) | matmult | t_linear | PDL ops |
|------:|---------:|---------:|--------:|---------:|--------:|
|     1 |   1633.7 |   1615.3 |  6220.0 |  27285.8 | 24136.8 |
|  1e2  |     23.7 |     29.5 |    85.4 |    423.9 |   297.2 |
|  1e3  |      7.7 |     12.3 |    25.2 |    128.6 |    51.5 |
|  1e4  |      6.2 |     10.6 |    19.0 |    119.9 |    70.2 |
|  1e5  |      6.2 |     10.5 |    18.2 |    150.2 |    77.6 |
|  1e6  |      4.6 |      5.1 |    18.5 |    153.4 |    75.9 |
|  1e7  |     16.2 |     17.1 |    29.1 |    154.0 |   105.8 |

### $`N`$ quaternions, $`N`$ vectors (one each)

|     N | rotq(PP) | PDL ops |
|------:|---------:|--------:|
|  1e3  |      7.8 |    52.6 |
|  1e5  |      6.8 |    31.6 |
|  1e7  |     20.5 |   118.8 |

### Allocation and threading ([alloc.pl](alloc.pl))

|     N | new out | prealloc | 1 thread, prealloc |
|------:|--------:|---------:|-------------------:|
|  1e5  |     6.2 |      6.1 |                6.1 |
|  1e6  |     5.3 |      3.8 |                6.6 |
|  3e6  |    17.7 |      4.0 |                7.5 |
|  1e7  |    22.9 |      4.2 |                7.9 |
|  3e7  |    22.6 |      3.8 |                7.0 |

### glibc mmap threshold ([mmap.pl](mmap.pl), `new out` only)

|     N | output size | default glibc | mmap/trim threshold raised |
|------:|------------:|--------------:|---------------------------:|
|  1e6  |       24 MB |           5.7 |                        7.5 |
|  3e6  |       72 MB |          20.8 |                        6.4 |
|  1e7  |      240 MB |          28.2 |                        6.1 |

### Page faults per call ([faults.pl](faults.pl))

The "pages in output" column is $`3 \cdot 8 \cdot N / 4096`$.

|     N | new out, default | new out, threshold raised | pages in output | prealloc |
|------:|-----------------:|--------------------------:|----------------:|---------:|
|  1e6  |                8 |                         8 |            5859 |        8 |
|  3e6  |           17 587 |                         8 |          17 578 |        8 |
|  1e7  |           58 602 |                         8 |          58 594 |        8 |

## Analysis

### H1: the PP kernel wins (confirmed)

- `rotq(PP)` is the fastest at every $`N`$.
- Plain matrix multiply (`x`) is 2–4× slower.
- The plain-PDL-ops version is 5–15× slower. It creates around ten
  full-size temporaries (two `crossp` results, scalings, sums), and
  each one is a separate pass over memory and a separate allocation.
- `t_linear` is the slowest of all for $`N \ge 10^2`$: 8–30× slower
  than `rotq`, and 5–8× slower than `x`, which it ultimately uses.
  - Its per-vector cost *rises* from about 120 to about 150 ns as $`N`$
    grows. That points to extra full-array copies inside
    `PDL::Transform::apply` (it handles generic input/output dims,
    `clump`/reshape, and so on), not to dispatch overhead.
  - `PDL::Transform` is a good tool for resampling, projections and
    composing transforms, but not for a rotation hot path.
- At $`N = 1`$, all PDL approaches cost about 1.6 µs or more per call,
  and that is purely fixed overhead:
  - creating the trans and piddles, and dispatching through XS;
  - `t_linear` and `PDL ops` cost 24–27 µs.

  So scalar code paths, single quaternions and single vectors, will be
  dominated by PDL call overhead.

### H2: quaternion kernel vs. matrix kernel (confirmed)

- While the data fits in cache ($`N \le 10^5`$), `rotq` is about 1.6–1.7×
  faster than `rotm`: 6.2 vs. 10.5 ns. `rotm`'s generic `loop(j)
  loop(i)` uses indexed access through increment variables. `rotq`'s
  fully unrolled straight-line code compiles better.
- Once the problem is memory-bound ($`N \ge 10^6`$), the two converge
  (4.6 vs. 5.1, then 16.2 vs. 17.1 ns), because both move the same 48
  bytes per vector.
- So there is no reason to convert a quaternion to a matrix just to
  rotate vectors. The matrix path is also slightly worse in the
  many-quaternion case, because it would have to build $`N`$ matrices
  first.

### H3: the large-$`N`$ slowdown is allocation, not the kernel (H3a confirmed, H3b rejected)

The first table shows a jump in `rotq` between $`10^6`$ (≈ 5 ns) and
$`10^7`$ (≈ 16–20 ns). Three measurements isolate the cause:

1. **Preallocating the output removes the jump entirely.** With
   preallocation, `rotq` stays at 3.8–4.2 ns/vector all the way to
   $`3\times10^7`$ vectors (720 MB of input plus output), well past L3.
   So the kernel and memory bandwidth are not the problem, and H3b is
   rejected.
2. **The jump follows the glibc mmap ceiling, not the cache size.**
   - At $`N=10^6`$ the output is 24 MB. That is far larger than L3,
     but below glibc's 32 MB dynamic mmap-threshold ceiling. After the
     first free, the threshold adapts upward and the block is reused
     from the heap, so the call is fast.
   - At $`3\times10^6`$ the output is 72 MB, which is above the
     ceiling. Every call gets a fresh `mmap` and gives it back with
     `munmap` on free, so the call is slow.
   - Raising `MALLOC_MMAP_THRESHOLD_` and `MALLOC_TRIM_THRESHOLD_` makes
     `new out` as fast as preallocated output: 6.1–6.4 ns.
3. **Page faults match exactly.**
   - With default glibc settings, each `rotq` call at $`N \ge 3\times10^6`$
     takes one minor page fault per 4 KiB page of output: 17 587 faults
     for 17 578 pages, and 58 602 for 58 594.
   - With preallocation, or with the threshold raised, each call takes
     about 8 faults regardless of $`N`$.

   At roughly 1 µs or more per fault (the kernel zeroes the page and
   updates the page tables), this accounts for the 12–18 ns/vector
   penalty: 24 bytes per vector is about 0.006 pages per vector.
   Transparent huge pages (`always`) evidently do not help this
   allocation pattern.

### Threading

- The autopthread threshold is 1 M elements. $`N=10^5`$ vectors is
  3×10⁵ elements, so it runs on one thread, and all three columns
  agree.
- From $`N=10^6`$ up, 8 threads give only about 1.7–1.9× over 1 thread
  (3.8–4.2 vs. 6.6–7.9 ns).
- The kernel moves about 72 B per vector: 24 read, 24 written, plus 24
  read for write-allocate. That is about 18–19 GB/s with 8 threads and
  about 10 GB/s with one. That is close to what this laptop's
  dual-channel DDR4 delivers in practice, so the multithreaded kernel
  is memory-bandwidth-bound.
- The in-cache single-thread figure, about 6 ns or roughly 12 cycles
  per vector for about 30 flops, suggests the generic strided indexing
  prevents SIMD vectorisation. That is a possible later optimisation,
  but it is not worth doing while the problem is bandwidth-bound.

## Conclusions

1. **Write the `CXC::Math::Quaternion` vector core in PDL::PP.** It is
   2–4× faster than `x` and 8–30× faster than `PDL::Transform`. It also
   avoids every temporary that a composition of PDL operations
   allocates.
2. **Rotate directly with the quaternion formula; do not go through a
   rotation matrix.** It is 1.7× faster in cache and equal when
   memory-bound. Matrices remain useful as an output (`transform`), not
   as the rotation path.
3. **Use `PDL::Transform` only for what it uniquely provides**:
   resampling/`map`, projections and composition. A thin
   `to_transform` bridge is enough. Do not use it for bulk rotation.
4. **Let PP functions accept a preallocated `[o]` output.** PDL::PP
   already allows this. Document it as the way to avoid the mmap and
   page-fault penalty for large repeated calls (outputs over 32 MB,
   i.e. more than about 1.4 M 3-vectors or 1 M quaternions). Users who
   cannot preallocate can raise `MALLOC_MMAP_THRESHOLD_`. That is an
   environment-level setting, and the module should not set it.
5. **Single-element operations are dominated by about 1.6 µs of PDL
   call overhead.** That supports keeping a pure-Perl or scalar path
   for the scalar API (single quaternion, single vector). How it
   compares with pure Perl was not measured here and should be checked
   before deciding.
6. **Autopthreading helps less than 2× on this machine,** because the
   kernel is bandwidth-bound. No custom threading work is warranted.

## Column-layout input

### Question

The vector data used above is `(3, N)`, so each vector's x, y and z
are next to each other. Real input often arrives as **three
independent `(N)` ndarrays**, for example table columns read with
`rcols`. Turning those into a `(3, N)` ndarray costs extra work. Two
things need measuring:

1. What does it cost to combine three `(N)` columns into a `(3, N)`
   ndarray?
2. How does "combine, then `rotq`" compare with a kernel that takes
   the three columns directly?

### Hypotheses

- **H4: building `(3, N)` costs at least as much as the rotation.**
  Building means one full read and write of the data. That is the same
  memory traffic as `rotq` itself, so it should at least double the
  total time. The usual idioms (`cat`, `glue`) should cost more,
  because they create intermediate ndarrays.
- **H5: a column kernel is about as fast as `rotq` on `(3, N)`.** It
  does the same arithmetic over the same 48 bytes per vector, just
  spread over six streams instead of two.
- **H6: the stride inside a vector does not matter.** PP kernels read
  through arbitrary strides, so `rotq` on a transposed view of an
  `(N, 3)` ndarray (`$v->xchg(0,1)`, stride $`N`$ between components)
  should be about as fast as on `(3, N)`. This would make a 2-D column
  store usable for free. Only three *separate* ndarrays force a copy,
  because PDL cannot present three separate buffers as one ndarray.

### Benchmarks ([cols.pl](cols.pl))

`@c = ($x, $y, $z)` are three independent `random(N) - 0.5`
ndarrays.

The first table measures building a physical `(3, N)` ndarray:

| Label           | Code |
|-----------------|------|
| `cat+xchg+copy` | `cat(@c)->xchg(0,1)->copy`. `cat` makes `(N, 3)`, then the transposed view is copied. |
| `glue`          | `$x->dummy(0)->glue(0, $y->dummy(0), $z->dummy(0))` |
| `fill prealloc` | `$v->slice("($_)") .= $c[$_] for 0..2`, into an existing `zeroes(3, N)` |

The second table measures getting the rotated result:

| Label           | Code |
|-----------------|------|
| `rotq(3,N)`     | Baseline. Data is already `(3, N)`, so nothing is built. |
| `rotq(N,3)^T`   | Data already stored as `(N, 3)`. Calls `rotq($q, $vt->xchg(0,1))` on a view, with no copy. |
| `cat+xchg`      | `rotq($q, cat(@c)->xchg(0,1))`. One copy (`cat`), then `rotq` reads the transposed view. |
| `cat+xchg+copy`, `glue` | Build as in the first table, then `rotq`. |
| `fill prealloc` | Fill a preallocated `(3, N)` ndarray, then `rotq` into a preallocated output. |
| `rotq3`         | `rotq3($q, @c)`, which allocates three output columns. |
| `rotq3 prealloc`| `rotq3($q, @c, @out)`, with three preallocated output columns. |

All variants agree exactly with each other (max $`\lvert \Delta\rvert = 0`$).

I ran each configuration 3 times, alternating between default glibc
malloc and `MALLOC_MMAP_THRESHOLD_` / `MALLOC_TRIM_THRESHOLD_` = 4 GiB
(which removes the page-fault penalty described under H3). Each run
used ≥ 0.5 s per cell. The tables give the **per-cell minimum over the
3 runs**, in ns per vector. The machine was busy: load average 3–5,
with the benchmark itself using up to 8 threads.
- Raw output: `results/cols-{default,nommap}-{1,2,3}.txt`.
- Minima: [results/cols-default-min.txt](results/cols-default-min.txt)
  and [results/cols-nommap-min.txt](results/cols-nommap-min.txt),
  produced by [results/min.pl](results/min.pl).

The large-$`N`$ results from the main benchmark held up: `rotq(3,N)`
measured 4.2 ns at $`10^6`$ and 16.8 ns at $`10^7`$, against 4.6 and
16.2 earlier. So these runs are comparable with the earlier tables.

### Results

**(1) Building `(3, N)` from three columns**, default malloc (and with
the mmap threshold raised):

|     N | cat+xchg+copy |        glue | fill prealloc |
|------:|--------------:|------------:|--------------:|
|  1e2  |  278.7 (267.4) | 338.2 (317.6) |  138.3 (141.0) |
|  1e3  |   62.9  (60.9) |  41.7  (43.2) |   16.9  (18.2) |
|  1e4  |   60.9  (42.5) |  27.5  (12.1) |    4.4   (4.4) |
|  1e5  |   67.1  (37.7) |   9.1   (9.0) |    3.1   (3.0) |
|  1e6  |   69.9  (43.0) |   9.4   (9.6) |    6.3   (6.1) |
|  1e7  |   72.0  (37.5) |  30.0   (6.8) |    7.1   (7.1) |

For comparison, `rotq` on data that is already `(3, N)` costs 4–7 ns
per vector for $`N \ge 10^4`$.

**(2) Rotating, default malloc:**

|     N | rotq(3,N) | rotq(N,3)^T | cat+xchg | cat+xchg+copy |  glue | fill prealloc | rotq3 | rotq3 prealloc |
|------:|----------:|------------:|---------:|--------------:|------:|--------------:|------:|---------------:|
|  1e2  |      27.9 |        40.4 |    232.1 |         300.4 | 349.6 |         161.9 |  47.1 |           23.6 |
|  1e3  |       8.5 |         9.8 |     31.2 |          70.9 |  51.8 |          25.8 |  11.4 |            7.8 |
|  1e4  |       6.3 |         6.6 |     11.5 |          47.2 |  17.8 |          10.6 |   7.3 |            6.3 |
|  1e5  |       6.8 |         6.1 |     10.3 |          49.4 |  16.7 |           9.1 |   7.1 |            6.1 |
|  1e6  |       4.2 |         4.4 |     37.6 |          83.6 |  49.4 |          10.1 |   7.7 |            6.3 |
|  1e7  |      16.8 |        17.4 |     34.2 |          82.6 |  49.6 |          10.6 |  19.4 |            3.7 |

**(2') Rotating, mmap threshold raised** (allocation page faults
removed):

|     N | rotq(3,N) | rotq(N,3)^T | cat+xchg | cat+xchg+copy |  glue | fill prealloc | rotq3 | rotq3 prealloc |
|------:|----------:|------------:|---------:|--------------:|------:|--------------:|------:|---------------:|
|  1e2  |      25.4 |        37.6 |    215.6 |         281.6 | 326.7 |         154.2 |  43.8 |           21.6 |
|  1e3  |       7.6 |         8.9 |     28.6 |          65.8 |  45.2 |          22.7 |  10.2 |            7.3 |
|  1e4  |       5.8 |         6.0 |     11.9 |          50.6 |  20.1 |          10.8 |   7.6 |            6.9 |
|  1e5  |       6.5 |         6.1 |      9.4 |          45.3 |  16.2 |           9.6 |   7.1 |            6.1 |
|  1e6  |       5.2 |         4.9 |     10.1 |          49.4 |  14.7 |           9.5 |   7.6 |            6.6 |
|  1e7  |       4.6 |         4.5 |      9.1 |          47.6 |  12.3 |          10.5 |   4.8 |            3.7 |

### Analysis

**H4: building `(3, N)` is expensive (confirmed).**
- Even the cheapest method, filling a preallocated ndarray through
  three slice assignments, costs 3–7 ns per vector for
  $`N \ge 10^4`$. That is about as much as the rotation, so
  "fill + rotate" takes 9–11 ns against 4–7 ns for the rotation alone.
- The idiomatic one-liners are far worse:
  - `cat(...)->xchg(0,1)->copy` costs 38–72 ns per vector, 6–15× the
    rotation. It copies twice, and the second copy gathers through the
    transposed view with stride $`N`$.
  - `glue` costs about 9 ns in the middle range. At $`N \le 10^3`$ it
    is slow (fixed Perl-level overhead from `dummy` and `glue`), and at
    $`N = 10^7`$ it is slow (allocation page faults).
- At small $`N`$, every extra PDL operation adds about 1–3 µs of fixed
  overhead per call, so building dominates completely. At $`N = 100`$,
  the cheapest build-then-rotate path (161.9 ns) is 7× slower than
  `rotq3 prealloc` (23.6 ns).

**H5: the column kernel matches `rotq` (confirmed).**
- With preallocated outputs, `rotq3` is as fast as `rotq` on `(3, N)`
  at every $`N`$. At $`10^7`$ it is slightly faster: 3.7 against 4.6 ns
  with the threshold raised.
- With freshly allocated outputs, it is 1.1–1.7× slower at small and
  moderate $`N`$. It creates three output ndarrays instead of one, so
  it pays the fixed per-ndarray cost three times: 47 against 28 ns at
  $`N = 100`$.
- At $`10^7`$ with default malloc, each 80 MB output column is mmapped
  fresh, so `rotq3` gets the same page-fault penalty as `rotq` (19.4
  against 16.8 ns). Raising the threshold removes it (4.8 ns).
- Against the fastest build-then-rotate path at each $`N`$
  (`fill prealloc` or `cat+xchg`):
  - With preallocated outputs on both sides, `rotq3` is 1.5–6.9×
    faster.
  - With fresh outputs, it is 1.3–3.4× faster up to $`10^6`$.
  - The one exception is $`10^7`$ with default malloc. There,
    `fill prealloc` (10.6 ns) beats fresh-output `rotq3` (19.4 ns),
    but only because it reuses a preallocated output and so avoids the
    page faults.

**H6: the stride inside a vector does not matter (confirmed).**
- `rotq` on the transposed view of an `(N, 3)` ndarray runs within
  about 10% of `rotq` on `(3, N)` for $`N \ge 10^3`$.
- The 12 ns extra at $`N = 100`$ is the fixed cost of creating the
  `xchg` view, about 1.2 µs per call.
- So the layout itself costs nothing. The cost comes only from
  *copying* separate columns into one ndarray.
- `cat+xchg` is therefore the best build-then-rotate option: one copy
  (`cat`) and no transposition. It costs 9–10 ns per vector once
  allocation is not a factor.
- With default malloc, `cat+xchg` jumps to about 35 ns at
  $`N \ge 10^6`$. The raised-threshold run brings it back to 10 ns, so
  this is again the allocation page-fault effect. Two 24 MB blocks are
  alive at the same time (the `cat` result and the output), and they
  are not reused from the heap.

### Conclusions for column-layout data

1. **Offer an API that takes the three columns directly, backed by
   column-signature PP kernels like `rotq3`.** It avoids a build step
   that costs at least as much as the rotation itself, and 2–15× more
   with the obvious `cat`/`glue` idioms. Example names:
   `rotate_xyz($q, $x, $y, $z)`, returning `($x', $y', $z')`.
   - It should also accept preallocated output columns.
   - Its outputs are columns, so callers who want columns back pay
     nothing on the way out either.
2. **Keep the `(3, N)` kernels as well.** They are the natural form
   for the rest of the API (vectors as one ndarray, broadcasting with
   `(4, N)` quaternions).
   - An `(N, 3)` ndarray can be passed as `$v->xchg(0,1)` at no extra
     cost, so 2-D column stores are already handled and need only a
     note in the documentation.
   - The API should not guess the layout from the shape, because
     $`N = 3`$ is ambiguous. Use separate entry points or an explicit
     option.
3. **If someone must build `(3, N)`, document the cheap options:**
   - `cat(@c)->xchg(0,1)`, as a view passed straight to the kernel;
   - slice-assigning into a preallocated ndarray.

   Avoid `->copy` on a transposed view. It was the slowest method
   measured.
4. **A column kernel triples the per-call ndarray overhead for outputs
   it allocates itself.** At small $`N`$ that is about 2 µs more per
   call than `rotq`. This only matters for many small calls, and
   passing preallocated outputs avoids it.
5. **Writing a column kernel costs little.** It is the same `Code`
   block with `$a()`/`$b()`/`$c()` in place of `$v(m=>0..2)`. If many
   operations need both forms, they could be generated from one
   template in the `.pd` file.

## Allocating outputs without zero-filling

### Question and hypothesis

`zeroes` writes zeros over the whole new ndarray. That write is wasted
when a PP kernel is about to overwrite every element.

- **H7:** `PDL->new_from_specification(double, @dims)` allocates
  without filling. It should therefore save one pass over the output
  and be faster than `zeroes`, both on its own and as a freshly
  allocated output passed to `rotq`/`rotq3`.

(The preallocated cases in the earlier tables are not affected either
way. Those outputs were allocated once, outside the timed loop.)

### Benchmark ([nfs.pl](nfs.pl))

- `nfs(@dims)` is `PDL->new_from_specification(double, @dims)` followed
  by `make_physical`, which forces the data to be allocated.
- Columns:
  - `rotq`: the kernel allocates its own output.
  - `rotq zeroes`, `rotq nfs`: an output is allocated by that method
    on every call and passed in.
  - `rotq pre`: one output allocated outside the loop.
  - The `rotq3` columns are the same, with three `(N)` output columns.
- Page faults are counted per call, as in [faults.pl](faults.pl).
- Same protocol as the column benchmarks: 3 runs per malloc setting,
  per-cell minimum.
  - Raw output: `results/nfs-{default,nommap}-{1,2,3}.txt`.
  - Minima: [results/nfs-default-min.txt](results/nfs-default-min.txt)
    and [results/nfs-nommap-min.txt](results/nfs-nommap-min.txt).

### Results

**Allocation alone**, ns per vector (3 doubles), and minor page faults
per call. Each cell gives default malloc, then the raised mmap
threshold in parentheses:

|     N | zeroes      | nfs         | faults, zeroes | faults, nfs |
|------:|------------:|------------:|---------------:|------------:|
|  1e3  |  2.2  (2.4) |  2.4  (2.7) |       0  (0)   |     0  (0)  |
|  1e5  |  0.6  (0.6) |  0.6  (0.6) |       0  (0)   |     0  (0)  |
|  1e6  |  1.1  (1.1) |  1.1  (1.1) |       0  (0)   |     0  (0)  |
|  3e6  | 15.2  (1.1) | 14.9  (1.1) |  17 579  (0)   | 17 579  (0) |
|  1e7  | 15.6  (1.1) | 15.5  (1.1) |  58 594  (0)   | 58 594  (0) |

**Rotating with a freshly allocated output on every call**, default
malloc:

|     N | rotq | rotq zeroes | rotq nfs | rotq pre | rotq3 | rotq3 zeroes | rotq3 nfs | rotq3 pre |
|------:|-----:|------------:|---------:|---------:|------:|-------------:|----------:|----------:|
|  1e3  |  8.6 |        10.0 |     10.4 |      7.4 |  11.4 |         13.7 |      14.5 |       8.2 |
|  1e5  |  6.4 |         6.5 |      6.4 |      5.8 |   7.1 |          6.8 |       7.0 |       6.6 |
|  1e6  |  5.1 |         5.2 |      5.2 |      3.8 |   8.0 |          8.2 |       8.1 |       6.7 |
|  3e6  | 18.8 |        19.4 |     18.9 |      3.8 |  15.1 |         14.7 |      14.9 |       3.8 |
|  1e7  | 19.6 |        19.6 |     19.4 |      3.6 |  20.7 |         19.3 |      19.5 |       3.8 |

**The same, mmap threshold raised:**

|     N | rotq | rotq zeroes | rotq nfs | rotq pre | rotq3 | rotq3 zeroes | rotq3 nfs | rotq3 pre |
|------:|-----:|------------:|---------:|---------:|------:|-------------:|----------:|----------:|
|  1e3  |  8.9 |         9.7 |      9.7 |      7.0 |  11.2 |         14.0 |      14.3 |       8.0 |
|  1e5  |  6.5 |         6.5 |      6.5 |      6.3 |   7.2 |          6.9 |       7.0 |       6.5 |
|  1e6  |  5.2 |         5.3 |      5.3 |      3.9 |   7.8 |          8.0 |       8.1 |       7.0 |
|  3e6  |  4.9 |         5.0 |      4.9 |      3.6 |   5.4 |          5.3 |       5.4 |       3.9 |
|  1e7  |  4.9 |         4.8 |      4.9 |      3.5 |   5.2 |          5.2 |       5.1 |       3.7 |

### Analysis

**H7 is rejected for PDL 2.106: `new_from_specification` +
`make_physical` costs exactly as much as `zeroes`.**
- Both cost the same at every $`N`$, with or without the mmap
  threshold.
- Both take one page fault for every 4 KiB page of the new ndarray. So
  allocating already *writes every page*.
- With the threshold raised, nothing faults, yet both still cost about
  1.1 ns per vector. That is 24 bytes per 1.1 ns, about 22 GB/s, which
  is a full write pass at this machine's memory bandwidth.
- The memory from `new_from_specification` also comes back all zero
  (max $`\lvert x\rvert = 0`$ at $`N = 10^6`$), even though PDL does not promise
  that.

So **PDL's own data allocation fills new ndarrays with zeros**, and
`zeroes` adds nothing on top of it. The PDL 2.106 source confirms this. The links below go to the
[`2.106` tag](https://github.com/PDLPorters/pdl/tree/2.106) on GitHub,
which is identical to the CPAN release for the files cited.

- **`zeroes` *is* `new_from_specification`.** In
  [Core.xs:198-202](https://github.com/PDLPorters/pdl/blob/2.106/lib/PDL/Core.xs#L198-L202),
  `PDL::zeroes` is an XS `ALIAS` of `new_from_specification`, with
  `ix = 1`. The only difference between them is one line,
  [Core.xs:249](https://github.com/PDLPorters/pdl/blob/2.106/lib/PDL/Core.xs#L249):
  `if (ix) pdl_make_physical(p)`.
  - So `new_from_specification` on its own allocates nothing: it only
    sets the dims.
  - `zeroes` allocates immediately.
  - There is no separate "write zeros" step in `zeroes`. The Changes
    entry "zeroes now quicker as only does memset, not assgn operation"
    records when that step was removed.
- **Every allocation clears memory.** `pdl_make_physical` ends up in
  `pdl_allocdata`
  ([pdlapi.c:172-208](https://github.com/PDLPorters/pdl/blob/2.106/lib/PDL/Core/pdlapi.c#L172-L208)).
  That function grows the Perl SV that holds the data (`SvGROW`, i.e.
  `malloc`). It then **always** runs
  `memset(it->data + ncurr, 0, nbytes - ncurr)` over the newly
  allocated bytes
  ([pdlapi.c:204](https://github.com/PDLPorters/pdl/blob/2.106/lib/PDL/Core/pdlapi.c#L204)).
  The Changes entries "calloc->malloc+memset, as calloc on Win32 = free
  to wrong pool" and "memset after malloc in pdlapi.c" explain why it
  is not `calloc`. With `calloc`, glibc would skip clearing fresh mmap
  pages.
- **No route avoids it.** A `new_from_specification` ndarray passed as
  an `[o]` argument without `make_physical` is still unallocated. The
  PP code then calls `pdl_make_physical` → `pdl_allocdata` → `memset`
  before the kernel runs. The same happens for outputs the kernel
  creates itself. So `zeroes`, `new_from_specification` (with or
  without `make_physical`) and kernel-allocated outputs all pay the
  same one clearing pass. That is exactly what the timings show.
  - The `memset` is also what triggers the page faults counted above:
    it is the first touch of each fresh mmap page, before the kernel
    writes anything.

The outputs that `rotq` and `rotq3` allocate themselves go through the
same path. Their cost matches the `zeroes`/`nfs` columns at large
$`N`$ (for example 19.6 / 19.6 / 19.4 ns at $`10^7`$). The pass you
wanted to remove therefore happens in every fresh-output case,
whichever way the output is created.

**What the extra pass costs.** With page faults out of the picture
(threshold raised, $`N \ge 10^6`$), a fresh output costs 1.3–1.5 ns per
vector more than a preallocated one: for example 4.9 against 3.5 ns at
$`10^7`$. That is roughly 25–40% of the kernel's time. With default
malloc and outputs over 32 MB, the page faults dominate (about 15 ns
per vector) and the fill is a small part of the total.

**At small $`N`$**, an output the caller creates on every call is
*slower* than one the kernel creates. At $`10^3`$ the figures are 10.0
against 8.6 ns for `rotq`, and 13.7 against 11.4 ns for `rotq3`. The
*data* allocation is identical in both cases: the same `pdl_allocdata`,
`malloc` and `memset`. The difference lies in **creating the ndarray
object around it**.

I measured this per call with [ovh.pl](ovh.pl) and
[ovh2.pl](ovh2.pl), at $`N = 1`$ where the data work is negligible.
Results are the per-cell minimum of 3 runs,
[results/ovh-min.txt](results/ovh-min.txt):

| Operation (ns per call, includes destroying the result)  |  N = 1 | N = 1000 |
|-----------------------------------------------------------|-------:|---------:|
| `zeroes(3,N)` (Core.pm wrapper)                            |   1255 |     1998 |
| `PDL->zeroes(3,N)` (XS directly)                           |   1168 |     1871 |
| `PDL->new_from_specification(double,3,N)`, no data         |   1043 |     1042 |
| … then `->make_physical`                                   |   1158 |     1865 |
| `PDL->null`                                                | ≈ 850 |        – |
| `rotq($q,$v)`, kernel creates output                       |   1580 |     7864 |
| `rotq($q,$v,$r)`, output reused                            |   1136 |     6490 |
| `rotq($q,$v,zeroes(3,N))`                                  |   2413 |     9073 |
| `rotq($q,$v,PDL->new_from_specification(double,3,N))`      |   2169 |     8895 |

- **The fixed cost of `rotq` itself is about 1.1 µs.** That is the
  "output reused" case. It covers creating the trans, sizing dims and
  broadcast setup, and the XS call.
- **The kernel adds about 0.44 µs to create its own output** (1580
  against 1136).
- **A caller-made output adds about 1.0–1.3 µs** (2169–2413 against
  1136). That is the measured cost of
  `new_from_specification`/`zeroes` alone.
- **At $`N = 1`$, almost none of that is data.**
  `new_from_specification` without any data allocation costs 1043 ns,
  and allocating the data (`make_physical`) adds only about 100 ns.
  Even `PDL->null` costs about 850 ns.
- [ovh2.pl](ovh2.pl) splits `new_from_specification` into construction
  (about 650–980 ns) and destruction (about 270–330 ns). Destruction
  happens in both paths, so it is not part of the difference.

The construction paths differ as follows:

- **Caller-made output.** The caller goes through:
  - Perl sub dispatch (the `zeroes` wrapper in Core.pm, then the
    `PDL->` class-method lookup);
  - argument checks and the `PDL::Type` object in
    `new_from_specification`
    ([Core.xs:199-226](https://github.com/PDLPorters/pdl/blob/2.106/lib/PDL/Core.xs#L199-L226));
  - building a Perl array of dims from the arguments, then unpacking it
    to C with `_dims_from_args` and `pdl_packdims`
    ([Core.xs:229-245](https://github.com/PDLPorters/pdl/blob/2.106/lib/PDL/Core.xs#L229-L245));
  - `pdl_setdims`.

  Then `rotq` must check the supplied output's type and dims against
  its signature.
- **Kernel-made output.** It is created in C by `PDL_XS_pdlinit`
  ([pdlperl.h:24-34](https://github.com/PDLPorters/pdl/blob/2.106/lib/PDL/Core/pdlperl.h#L24-L34)).
  For a plain `PDL` class this is a shortcut: `pdlnew()`, then a
  mortal SV, `SetSV_PDL` and `bless`. There is no Perl call, and no
  dims parsing, because the dims come straight from the trans in
  `redodims`.

So per output ndarray, the caller path pays about 0.6–0.9 µs more
object construction than the kernel path. Those few hundred ns to 1 µs
are all fixed Perl/XS bookkeeping. At $`N = 10^3`$ they show up as
1–2 ns per vector, with three times the effect for `rotq3`'s three
outputs. By $`N = 10^5`$ they are lost in the noise. The saving from
caller allocation only appears when **one output is reused across
calls**, so that the construction and the `memset` are paid once.

### Conclusions

1. **Replacing `zeroes` with `new_from_specification` gains nothing in
   PDL 2.106.** `zeroes` is `new_from_specification` plus an immediate
   `make_physical`. The clearing pass is the unconditional `memset` in
   `pdl_allocdata`, and every path that allocates data goes through
   it. Inside `CXC::Math::Quaternion`, use whichever reads better.
2. **The only way to avoid that pass today is to reuse an output
   ndarray across calls.** That is the "preallocated `[o]` argument"
   recommendation from earlier, and it removes the page faults too. It
   is worth about 25–40% at large $`N`$ without page faults, and 4–5×
   with them.
3. **Skipping the zero fill would need a change in PDL itself,** at
   [pdlapi.c:204](https://github.com/PDLPorters/pdl/blob/2.106/lib/PDL/Core/pdlapi.c#L204).
   For example, a state flag could let `pdl_allocdata` skip the
   `memset` for output ndarrays that a PP kernel fully overwrites (that
   is, a non-`[t]`, non-in-place `[o]`). That would be an upstream
   request. Based on the measurements here, it would save about 1.1–1.5
   ns per vector, roughly 25–40% of `rotq`, for fresh outputs at large
   $`N`$ without page faults. With default malloc and outputs over
   32 MB, the page faults would move from the `memset` to the kernel's
   first write. They would not go away, so the 4–5× penalty would
   mostly remain.

## Caveats and open questions

- These results come from a single machine: a laptop CPU inside a
  container, with frequency scaling. Absolute numbers will differ on
  servers. The ratios and the mmap behaviour, which is glibc-specific,
  should carry over to Linux/glibc. Other platforms (macOS, musl,
  Windows) have different allocators.
- Only `double` was benchmarked.
- The column-layout runs used only the one-quaternion case. Quaternions
  supplied as four columns, `(N)` each, were not measured. They would
  need a `q0(); q1(); q2(); q3()` signature variant.
- No pure-Perl baseline was measured, so the scalar-path decision in
  conclusion 5 is still open.
- Quaternion–quaternion products, normalisation and the
  RA/Dec/Roll ↔ quaternion conversions were not benchmarked. They are
  expected to behave like `rotq`: one pass, no temporaries.
