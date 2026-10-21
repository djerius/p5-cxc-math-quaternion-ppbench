# Rotating vectors with PDL: PDL::PP vs. PDL::Transform vs. plain PDL

I (Claude) wrote the code and this README.  Guidance was provided by Diab Jerius.

## Contents

- [Building and running](#building-and-running)
  - [Requirements](#requirements)
  - [Build](#build)
  - [Run](#run)
- [Hypotheses](#hypotheses)
- [Setup](#setup)
  - [Code](#code)
  - [Method](#method)
- [Results](#results)
  - [One quaternion, N vectors (bench.pl)](#one-quaternion-n-vectors-benchpl)
  - [N quaternions, N vectors (one each)](#n-quaternions-n-vectors-one-each)
  - [Allocation and threading (alloc.pl)](#allocation-and-threading-allocpl)
  - [glibc mmap threshold (mmap.pl, new out only)](#glibc-mmap-threshold-mmappl-new-out-only)
  - [Page faults per call (faults.pl)](#page-faults-per-call-faultspl)
- [Analysis](#analysis)
  - [H1: the PP kernel wins (confirmed)](#h1-the-pp-kernel-wins-confirmed)
  - [H2: quaternion kernel vs. matrix kernel (confirmed)](#h2-quaternion-kernel-vs-matrix-kernel-confirmed)
  - [H3: the large-N slowdown is allocation, not the kernel (H3a confirmed, H3b rejected)](#h3-the-large-n-slowdown-is-allocation-not-the-kernel-h3a-confirmed-h3b-rejected)
  - [Threading](#threading)
- [Conclusions](#conclusions)
- [Column-layout input](#column-layout-input)
  - [Question](#question)
  - [Hypotheses](#hypotheses-1)
  - [Benchmarks (cols.pl)](#benchmarks-colspl)
  - [Results](#results-1)
  - [Analysis](#analysis-1)
  - [Conclusions for column-layout data](#conclusions-for-column-layout-data)
- [Allocating outputs without zero-filling](#allocating-outputs-without-zero-filling)
  - [Question and hypothesis](#question-and-hypothesis)
  - [Benchmark (nfs.pl)](#benchmark-nfspl)
  - [Results](#results-2)
  - [Analysis](#analysis-2)
  - [Conclusions for output allocation](#conclusions-for-output-allocation)
- [Pure-Perl baseline](#pure-perl-baseline)
  - [Question](#question-1)
  - [Hypotheses](#hypotheses-2)
  - [Benchmark (pureperl.pl)](#benchmark-pureperlpl)
  - [Results](#results-3)
  - [Analysis](#analysis-3)
  - [Conclusions for the pure-Perl implementation](#conclusions-for-the-pure-perl-implementation)
- [Shared C implementation: XS on Perl data](#shared-c-implementation-xs-on-perl-data)
  - [Question](#question-2)
  - [Hypotheses](#hypotheses-3)
  - [Benchmark (xs.pl)](#benchmark-xspl)
  - [Results](#results-4)
  - [Analysis](#analysis-4)
  - [Conclusions for the XS implementation](#conclusions-for-the-xs-implementation)
- [Writing into existing Perl arrays](#writing-into-existing-perl-arrays)
  - [Question](#question-3)
  - [Hypotheses](#hypotheses-4)
  - [Benchmark (xsinto.pl, svlayout.pl)](#benchmark-xsintopl-svlayoutpl)
  - [Results](#results-5)
  - [Analysis](#analysis-5)
  - [Conclusions for writing into existing arrays](#conclusions-for-writing-into-existing-arrays)
- [Comparison with PDL::LinearAlgebra](#comparison-with-pdllinearalgebra)
  - [Question](#question-4)
  - [Hypotheses](#hypotheses-5)
  - [Benchmark (lalg.pl, autopthread.pl)](#benchmark-lalgpl-autopthreadpl)
  - [Results](#results-6)
  - [Analysis](#analysis-6)
  - [Conclusions for PDL::LinearAlgebra](#conclusions-for-pdllinearalgebra)
- [Argument dispatch](#argument-dispatch)
  - [Question](#question-5)
  - [Hypotheses](#hypotheses-6)
  - [Benchmark (dispatch.pl, dispatch-parts.pl)](#benchmark-dispatchpl-dispatch-partspl)
  - [Results](#results-7)
  - [Analysis](#analysis-7)
  - [Conclusions for argument dispatch](#conclusions-for-argument-dispatch)
- [Caveats and open questions](#caveats-and-open-questions)

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
- For [lalg.pl](lalg.pl) only: PDL::LinearAlgebra, which needs BLAS
  and LAPACK libraries and a Fortran compiler to build. The results
  here used PDL::LinearAlgebra 0.437 with OpenBLAS 0.3.26 (pthreads
  build, which also supplies LAPACK 3.12.0), installed on Ubuntu with
  `apt-get install libopenblas-dev liblapack-dev`, then
  `cpanm PDL::LinearAlgebra`.
- Linux with glibc for the allocation experiments:
  - [faults.pl](faults.pl) reads `/proc/self/stat`.
  - The `MALLOC_*` environment variables below are glibc-specific.

  The other scripts run anywhere PDL does.

### Build

```sh
perl Makefile.PL
make
```

This compiles the module `QBench` from [qbench.pd](qbench.pd) into
`blib/`. It contains the PP kernels (`rotq`, `rotm`, `rotq3`,
`rotq_shared`) and plain XS functions on Perl data (`xs_rot1`,
`xs_rot1_list`, `xs_rot_aoa`, `xs_rot_cols`, and the `_into` variants
`xs_rot1_into`, `xs_rot_aoa_into`, `xs_rot_cols_into`), plus the
dispatch prototype `xsd_rotate_vec`. `rotq_shared`
and the XS functions share the C code in [qrot.h](qrot.h).
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
perl pureperl.pl  # pure-Perl rotation vs. PDL, conversions        (~2 min)
perl xs.pl        # XS on Perl data vs. pure Perl; shared PP kernel (~70 s)
perl xsinto.pl    # writing into existing arrays vs. new arrays    (~70 s)
perl svlayout.pl clean; perl svlayout.pl churn   # scalar layout   (~30 s)
perl lalg.pl      # PDL::LinearAlgebra (BLAS) vs. PP kernel        (~2 min)
OPENBLAS_NUM_THREADS=1 perl lalg.pl                # BLAS on one thread
perl autopthread.pl  # rotq with PDL's threading threshold lowered (~20 s)
perl dispatch-test.pl   # check the two dispatch prototypes            (~2 s)
perl dispatch.pl        # cost of argument dispatch, Perl vs. XS       (~2 min)
perl dispatch-parts.pl  # cost of each step of the Perl dispatcher     (~10 s)
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

The other `results/*-min.txt` files were made the same way, from three
runs of [nfs.pl](nfs.pl), [ovh.pl](ovh.pl), [pureperl.pl](pureperl.pl),
[xs.pl](xs.pl), [xsinto.pl](xsinto.pl), [lalg.pl](lalg.pl) (3 runs with
each OpenBLAS threading setting, alternating) and
[autopthread.pl](autopthread.pl). (`results/svlayout.txt` is
simply the six runs of [svlayout.pl](svlayout.pl), in order.) `min.pl`
matches lines by position and replaces
every plain decimal number with its minimum over the runs, so it works
on any of the scripts' output.

Timings on a laptop or a shared machine vary by 10–30% from run to
run, so compare ratios rather than absolute numbers.

## Hypotheses

1. <a id="h1"></a>**H1: a dedicated PP kernel beats the generic alternatives.** A PP
   kernel does the whole rotation in one pass over the data with no
   temporaries. It should beat:
   - plain PDL operations, which allocate a temporary for every
     intermediate;
   - `x` (matrix multiply);
   - `PDL::Transform::t_linear`, which adds Perl-level dispatch and
     generic machinery.
2. <a id="h2"></a>**H2: rotating with the quaternion directly is at least as fast as
   building a 3×3 matrix and using that.** The quaternion form,
   $`v' = v + 2w(u \times v) + 2u \times (u \times v)`$ with
   $`u = (x,y,z)`$, needs about the same number of flops as a
   matrix–vector product. It also avoids the matrix construction and
   reads only 4 numbers instead of 9.
3. <a id="h3"></a>**H3 (added after the first run): the slowdown at large $`N`$ comes
   from the kernel or from memory bandwidth.**
   - <a id="h3a"></a>**H3a:** it is from memory allocation, not from the kernel. glibc
     `malloc` serves blocks above the dynamic mmap threshold with fresh
     `mmap` calls. Its default ceiling is 32 MB on 64-bit systems. A
     fresh mapping is page-faulted in one 4 KiB page at a time when the
     kernel first writes to it, and it is unmapped again when the
     output piddle is freed. So every call pays that cost again.
   - <a id="h3b"></a>**H3b (alternative):** it is the memory bandwidth limit, reached
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
| BLAS/LAPACK  | OpenBLAS 0.3.26 (pthreads), LAPACK 3.12.0; PDL::LinearAlgebra 0.437 ([lalg.pl](lalg.pl) only) |

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
- [pureperl.pl](pureperl.pl) times pure-Perl rotation in several
  forms. It compares them with the PP kernels, both with data already
  in ndarrays and with Perl data converted to ndarrays and back. See
  [Pure-Perl baseline](#pure-perl-baseline).
  - The pure-Perl implementations live in [QRotPerl.pm](QRotPerl.pm),
    shared with [xs.pl](xs.pl).
- [qrot.h](qrot.h) holds the rotation as one `static inline` C
  function, `qrot_vec`, with no Perl or PDL API. Two front ends in
  [qbench.pd](qbench.pd) use it:
  - the PP kernel `rotq_shared`;
  - plain XS functions that read and write Perl arrays and scalars.
- [xs.pl](xs.pl) compares those XS functions with the pure-Perl code,
  and `rotq_shared` with the hand-written `rotq`. See
  [Shared C implementation: XS on Perl data](#shared-c-implementation-xs-on-perl-data).
- [xsinto.pl](xsinto.pl) compares XS and pure-Perl functions that
  write into arrays the caller supplies (the `_into` variants) with the
  versions that return new arrays. [svlayout.pl](svlayout.pl) checks
  how the memory layout of Perl scalars affects XS on large arrays. See
  [Writing into existing Perl arrays](#writing-into-existing-perl-arrays).
- [lalg.pl](lalg.pl) compares PDL::LinearAlgebra's BLAS matrix
  multiply (`mmult`, `gemm`) with the PP kernel and with PDL's own
  `matmult`. [autopthread.pl](autopthread.pl) times `rotq` with PDL's
  threading threshold lowered. See
  [Comparison with PDL::LinearAlgebra](#comparison-with-pdllinearalgebra).
- [QDispatch.pm](QDispatch.pm) and `xsd_rotate_vec` in
  [qbench.pd](qbench.pd) prototype the API's argument dispatch, in pure
  Perl and in XS. [dispatch-test.pl](dispatch-test.pl) checks them,
  [dispatch.pl](dispatch.pl) times them against the kernels, and
  [dispatch-parts.pl](dispatch-parts.pl) times the steps of the Perl
  one. See [Argument dispatch](#argument-dispatch).

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

### One quaternion, N vectors ([bench.pl](bench.pl))

|     N | rotq(PP) | rotm(PP) | matmult | t_linear | PDL ops |
|------:|---------:|---------:|--------:|---------:|--------:|
| 1 | 1633.7 ns | 1615.3 ns | 6220.0 ns | 27285.8 ns | 24136.8 ns |
| 1e2 | 23.7 ns | 29.5 ns | 85.4 ns | 423.9 ns | 297.2 ns |
| 1e3 | 7.7 ns | 12.3 ns | 25.2 ns | 128.6 ns | 51.5 ns |
| 1e4 | 6.2 ns | 10.6 ns | 19.0 ns | 119.9 ns | 70.2 ns |
| 1e5 | 6.2 ns | 10.5 ns | 18.2 ns | 150.2 ns | 77.6 ns |
| 1e6 | 4.6 ns | 5.1 ns | 18.5 ns | 153.4 ns | 75.9 ns |
| 1e7 | 16.2 ns | 17.1 ns | 29.1 ns | 154.0 ns | 105.8 ns |

### N quaternions, N vectors (one each)

|     N | rotq(PP) | PDL ops |
|------:|---------:|--------:|
| 1e3 | 7.8 ns | 52.6 ns |
| 1e5 | 6.8 ns | 31.6 ns |
| 1e7 | 20.5 ns | 118.8 ns |

### Allocation and threading ([alloc.pl](alloc.pl))

|     N | new out | prealloc | 1 thread, prealloc |
|------:|--------:|---------:|-------------------:|
| 1e5 | 6.2 ns | 6.1 ns | 6.1 ns |
| 1e6 | 5.3 ns | 3.8 ns | 6.6 ns |
| 3e6 | 17.7 ns | 4.0 ns | 7.5 ns |
| 1e7 | 22.9 ns | 4.2 ns | 7.9 ns |
| 3e7 | 22.6 ns | 3.8 ns | 7.0 ns |

### glibc mmap threshold ([mmap.pl](mmap.pl), `new out` only)

|     N | output size | default glibc | mmap/trim threshold raised |
|------:|------------:|--------------:|---------------------------:|
| 1e6 | 24 MB | 5.7 ns | 7.5 ns |
| 3e6 | 72 MB | 20.8 ns | 6.4 ns |
| 1e7 | 240 MB | 28.2 ns | 6.1 ns |

### Page faults per call ([faults.pl](faults.pl))

The "pages in output" column is $`3 \cdot 8 \cdot N / 4096`$.

|     N | new out, default | new out, threshold raised | pages in output | prealloc |
|------:|-----------------:|--------------------------:|----------------:|---------:|
| 1e6 | 8 faults | 8 faults | 5859 pages | 8 faults |
| 3e6 | 17 587 faults | 8 faults | 17 578 pages | 8 faults |
| 1e7 | 58 602 faults | 8 faults | 58 594 pages | 8 faults |

## Analysis

### [H1](#h1): the PP kernel wins (confirmed)

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

### [H2](#h2): quaternion kernel vs. matrix kernel (confirmed)

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

### [H3](#h3): the large-N slowdown is allocation, not the kernel ([H3a](#h3a) confirmed, [H3b](#h3b) rejected)

The first table shows a jump in `rotq` between $`10^6`$ (≈ 5 ns) and
$`10^7`$ (≈ 16–20 ns). Three measurements isolate the cause:

1. **Preallocating the output removes the jump entirely.** With
   preallocation, `rotq` stays at 3.8–4.2 ns/vector all the way to
   $`3\times10^7`$ vectors (720 MB of input plus output), well past L3.
   So the kernel and memory bandwidth are not the problem, and [H3b](#h3b) is
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

Conclusions are numbered [C1](#c1)–[C34](#c34) through the whole
document, like the hypotheses ([H1](#h1)–[H19](#h19)), so later
sections can refer to earlier ones.

- <a id="c1"></a>**C1: Write the `CXC::Math::Quaternion` vector core in PDL::PP.** It is
  2–4× faster than `x` and 8–30× faster than `PDL::Transform`. It also
  avoids every temporary that a composition of PDL operations
  allocates.
- <a id="c2"></a>**C2: Rotate directly with the quaternion formula; do not go through a
  rotation matrix.** It is 1.7× faster in cache and equal when
  memory-bound. Matrices remain useful as an output (`transform`), not
  as the rotation path.
- <a id="c3"></a>**C3: Use `PDL::Transform` only for what it uniquely provides**:
  resampling/`map`, projections and composition. A thin
  `to_transform` bridge is enough. Do not use it for bulk rotation.
- <a id="c4"></a>**C4: Let PP functions accept a preallocated `[o]` output.** PDL::PP
  already allows this. Document it as the way to avoid the mmap and
  page-fault penalty for large repeated calls (outputs over 32 MB,
  i.e. more than about 1.4 M 3-vectors or 1 M quaternions). Users who
  cannot preallocate can raise `MALLOC_MMAP_THRESHOLD_`. That is an
  environment-level setting, and the module should not set it.
- <a id="c5"></a>**C5: Single-element operations are dominated by about 1.6 µs of PDL
  call overhead.** That supports keeping a pure-Perl path for the
  scalar API (single quaternion, single vector). The
  [Pure-Perl baseline](#pure-perl-baseline) confirms this: pure Perl
  is faster for one vector, and much faster when the vector starts out
  as Perl data ([C16](#c16)–[C18](#c18)).
- <a id="c6"></a>**C6: Autopthreading helps less than 2× on this machine,** because the
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

- <a id="h4"></a>**H4: building `(3, N)` costs at least as much as the rotation.**
  Building means one full read and write of the data. That is the same
  memory traffic as `rotq` itself, so it should at least double the
  total time. The usual idioms (`cat`, `glue`) should cost more,
  because they create intermediate ndarrays.
- <a id="h5"></a>**H5: a column kernel is about as fast as `rotq` on `(3, N)`.** It
  does the same arithmetic over the same 48 bytes per vector, just
  spread over six streams instead of two.
- <a id="h6"></a>**H6: the stride inside a vector does not matter.** PP kernels read
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

I (Claude) ran each configuration 3 times, alternating between default glibc
malloc and `MALLOC_MMAP_THRESHOLD_` / `MALLOC_TRIM_THRESHOLD_` = 4 GiB
(which removes the page-fault penalty described under [H3](#h3)). Each run
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

|     N | cat+xchg+copy | glue | fill prealloc |
|------:|--------------:|------------:|--------------:|
| 1e2 | 278.7 ns (267.4 ns) | 338.2 ns (317.6 ns) | 138.3 ns (141.0 ns) |
| 1e3 | 62.9 ns (60.9 ns) | 41.7 ns (43.2 ns) | 16.9 ns (18.2 ns) |
| 1e4 | 60.9 ns (42.5 ns) | 27.5 ns (12.1 ns) | 4.4 ns (4.4 ns) |
| 1e5 | 67.1 ns (37.7 ns) | 9.1 ns (9.0 ns) | 3.1 ns (3.0 ns) |
| 1e6 | 69.9 ns (43.0 ns) | 9.4 ns (9.6 ns) | 6.3 ns (6.1 ns) |
| 1e7 | 72.0 ns (37.5 ns) | 30.0 ns (6.8 ns) | 7.1 ns (7.1 ns) |

For comparison, `rotq` on data that is already `(3, N)` costs 4–7 ns
per vector for $`N \ge 10^4`$.

**(2) Rotating, default malloc:**

|     N | rotq(3,N) | rotq(N,3)^T | cat+xchg | cat+xchg+copy | glue | fill prealloc | rotq3 | rotq3 prealloc |
|------:|----------:|------------:|---------:|--------------:|------:|--------------:|------:|---------------:|
| 1e2 | 27.9 ns | 40.4 ns | 232.1 ns | 300.4 ns | 349.6 ns | 161.9 ns | 47.1 ns | 23.6 ns |
| 1e3 | 8.5 ns | 9.8 ns | 31.2 ns | 70.9 ns | 51.8 ns | 25.8 ns | 11.4 ns | 7.8 ns |
| 1e4 | 6.3 ns | 6.6 ns | 11.5 ns | 47.2 ns | 17.8 ns | 10.6 ns | 7.3 ns | 6.3 ns |
| 1e5 | 6.8 ns | 6.1 ns | 10.3 ns | 49.4 ns | 16.7 ns | 9.1 ns | 7.1 ns | 6.1 ns |
| 1e6 | 4.2 ns | 4.4 ns | 37.6 ns | 83.6 ns | 49.4 ns | 10.1 ns | 7.7 ns | 6.3 ns |
| 1e7 | 16.8 ns | 17.4 ns | 34.2 ns | 82.6 ns | 49.6 ns | 10.6 ns | 19.4 ns | 3.7 ns |

**(2') Rotating, mmap threshold raised** (allocation page faults
removed):

|     N | rotq(3,N) | rotq(N,3)^T | cat+xchg | cat+xchg+copy | glue | fill prealloc | rotq3 | rotq3 prealloc |
|------:|----------:|------------:|---------:|--------------:|------:|--------------:|------:|---------------:|
| 1e2 | 25.4 ns | 37.6 ns | 215.6 ns | 281.6 ns | 326.7 ns | 154.2 ns | 43.8 ns | 21.6 ns |
| 1e3 | 7.6 ns | 8.9 ns | 28.6 ns | 65.8 ns | 45.2 ns | 22.7 ns | 10.2 ns | 7.3 ns |
| 1e4 | 5.8 ns | 6.0 ns | 11.9 ns | 50.6 ns | 20.1 ns | 10.8 ns | 7.6 ns | 6.9 ns |
| 1e5 | 6.5 ns | 6.1 ns | 9.4 ns | 45.3 ns | 16.2 ns | 9.6 ns | 7.1 ns | 6.1 ns |
| 1e6 | 5.2 ns | 4.9 ns | 10.1 ns | 49.4 ns | 14.7 ns | 9.5 ns | 7.6 ns | 6.6 ns |
| 1e7 | 4.6 ns | 4.5 ns | 9.1 ns | 47.6 ns | 12.3 ns | 10.5 ns | 4.8 ns | 3.7 ns |

### Analysis

**[H4](#h4): building `(3, N)` is expensive (confirmed).**
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

**[H5](#h5): the column kernel matches `rotq` (confirmed).**
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

**[H6](#h6): the stride inside a vector does not matter (confirmed).**
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

- <a id="c7"></a>**C7: Offer an API that takes the three columns directly, backed by
  column-signature PP kernels like `rotq3`.** It avoids a build step
  that costs at least as much as the rotation itself, and 2–15× more
  with the obvious `cat`/`glue` idioms. Example names:
  `rotate_xyz($q, $x, $y, $z)`, returning `($x', $y', $z')`.
  - It should also accept preallocated output columns.
  - Its outputs are columns, so callers who want columns back pay
    nothing on the way out either.
- <a id="c8"></a>**C8: Keep the `(3, N)` kernels as well.** They are the natural form
  for the rest of the API (vectors as one ndarray, broadcasting with
  `(4, N)` quaternions).
  - An `(N, 3)` ndarray can be passed as `$v->xchg(0,1)` at no extra
    cost, so 2-D column stores are already handled and need only a
    note in the documentation.
  - The API should not guess the layout from the shape, because
    $`N = 3`$ is ambiguous. Use separate entry points or an explicit
    option.
- <a id="c9"></a>**C9: If someone must build `(3, N)`, document the cheap options:**
  - `cat(@c)->xchg(0,1)`, as a view passed straight to the kernel;
  - slice-assigning into a preallocated ndarray.

  Avoid `->copy` on a transposed view. It was the slowest method
  measured.
- <a id="c10"></a>**C10: A column kernel triples the per-call ndarray overhead for outputs
  it allocates itself.** At small $`N`$ that is about 2 µs more per
  call than `rotq`. This only matters for many small calls, and
  passing preallocated outputs ([C4](#c4)) avoids it.
- <a id="c11"></a>**C11: Writing a column kernel costs little.** It is the same `Code`
  block with `$a()`/`$b()`/`$c()` in place of `$v(m=>0..2)`. If many
  operations need both forms, they could be generated from one
  template in the `.pd` file.

## Allocating outputs without zero-filling

### Question and hypothesis

`zeroes` writes zeros over the whole new ndarray. That write is wasted
when a PP kernel is about to overwrite every element.

- <a id="h7"></a>**H7:** `PDL->new_from_specification(double, @dims)` allocates
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

|     N | zeroes | nfs | faults, zeroes | faults, nfs |
|------:|------------:|------------:|---------------:|------------:|
| 1e3 | 2.2 ns (2.4 ns) | 2.4 ns (2.7 ns) | 0 faults (0 faults) | 0 faults (0 faults) |
| 1e5 | 0.6 ns (0.6 ns) | 0.6 ns (0.6 ns) | 0 faults (0 faults) | 0 faults (0 faults) |
| 1e6 | 1.1 ns (1.1 ns) | 1.1 ns (1.1 ns) | 0 faults (0 faults) | 0 faults (0 faults) |
| 3e6 | 15.2 ns (1.1 ns) | 14.9 ns (1.1 ns) | 17 579 faults (0 faults) | 17 579 faults (0 faults) |
| 1e7 | 15.6 ns (1.1 ns) | 15.5 ns (1.1 ns) | 58 594 faults (0 faults) | 58 594 faults (0 faults) |

**Rotating with a freshly allocated output on every call**, default
malloc:

|     N | rotq | rotq zeroes | rotq nfs | rotq pre | rotq3 | rotq3 zeroes | rotq3 nfs | rotq3 pre |
|------:|-----:|------------:|---------:|---------:|------:|-------------:|----------:|----------:|
| 1e3 | 8.6 ns | 10.0 ns | 10.4 ns | 7.4 ns | 11.4 ns | 13.7 ns | 14.5 ns | 8.2 ns |
| 1e5 | 6.4 ns | 6.5 ns | 6.4 ns | 5.8 ns | 7.1 ns | 6.8 ns | 7.0 ns | 6.6 ns |
| 1e6 | 5.1 ns | 5.2 ns | 5.2 ns | 3.8 ns | 8.0 ns | 8.2 ns | 8.1 ns | 6.7 ns |
| 3e6 | 18.8 ns | 19.4 ns | 18.9 ns | 3.8 ns | 15.1 ns | 14.7 ns | 14.9 ns | 3.8 ns |
| 1e7 | 19.6 ns | 19.6 ns | 19.4 ns | 3.6 ns | 20.7 ns | 19.3 ns | 19.5 ns | 3.8 ns |

**The same, mmap threshold raised:**

|     N | rotq | rotq zeroes | rotq nfs | rotq pre | rotq3 | rotq3 zeroes | rotq3 nfs | rotq3 pre |
|------:|-----:|------------:|---------:|---------:|------:|-------------:|----------:|----------:|
| 1e3 | 8.9 ns | 9.7 ns | 9.7 ns | 7.0 ns | 11.2 ns | 14.0 ns | 14.3 ns | 8.0 ns |
| 1e5 | 6.5 ns | 6.5 ns | 6.5 ns | 6.3 ns | 7.2 ns | 6.9 ns | 7.0 ns | 6.5 ns |
| 1e6 | 5.2 ns | 5.3 ns | 5.3 ns | 3.9 ns | 7.8 ns | 8.0 ns | 8.1 ns | 7.0 ns |
| 3e6 | 4.9 ns | 5.0 ns | 4.9 ns | 3.6 ns | 5.4 ns | 5.3 ns | 5.4 ns | 3.9 ns |
| 1e7 | 4.9 ns | 4.8 ns | 4.9 ns | 3.5 ns | 5.2 ns | 5.2 ns | 5.1 ns | 3.7 ns |

### Analysis

**[H7](#h7) is rejected for PDL 2.106: `new_from_specification` +
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

I (Claude) measured this per call with [ovh.pl](ovh.pl) and
[ovh2.pl](ovh2.pl), at $`N = 1`$ where the data work is negligible.
Results are the per-cell minimum of 3 runs,
[results/ovh-min.txt](results/ovh-min.txt):

| Operation (includes destroying the result) | N = 1 | N = 1000 |
|-----------------------------------------------------------|-------:|---------:|
| `zeroes(3,N)` (Core.pm wrapper) | 1255 ns | 1998 ns |
| `PDL->zeroes(3,N)` (XS directly) | 1168 ns | 1871 ns |
| `PDL->new_from_specification(double,3,N)`, no data | 1043 ns | 1042 ns |
| … then `->make_physical` | 1158 ns | 1865 ns |
| `PDL->null` | ≈ 850 ns | – |
| `rotq($q,$v)`, kernel creates output | 1580 ns | 7864 ns |
| `rotq($q,$v,$r)`, output reused | 1136 ns | 6490 ns |
| `rotq($q,$v,zeroes(3,N))` | 2413 ns | 9073 ns |
| `rotq($q,$v,PDL->new_from_specification(double,3,N))` | 2169 ns | 8895 ns |

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

### Conclusions for output allocation

- <a id="c12"></a>**C12: Replacing `zeroes` with `new_from_specification` gains nothing in
  PDL 2.106.** `zeroes` is `new_from_specification` plus an immediate
  `make_physical`. The clearing pass is the unconditional `memset` in
  `pdl_allocdata`, and every path that allocates data goes through
  it. Inside `CXC::Math::Quaternion`, use whichever reads better.
- <a id="c13"></a>**C13: The only way to avoid that pass today is to reuse an output
  ndarray across calls.** That is [C4](#c4) (accept a preallocated `[o]`
  output), which also removes the page faults. It
  is worth about 25–40% at large $`N`$ without page faults, and 4–5×
  with them.
- <a id="c14"></a>**C14: Skipping the zero fill would need a change in PDL itself,** at
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

## Pure-Perl baseline

### Question

`CXC::Math::Quaternion` will have a pure-Perl implementation, for use
when PDL is not installed. That raises three questions:

1. How fast is pure Perl, and in what form should it be written?
2. At what $`N`$ does PDL start to win?
3. When the caller's data is in Perl arrays, is it worth converting it
   to ndarrays, running the PP kernel and converting back?

### Hypotheses

- <a id="h8"></a>**H8: pure Perl wins for one vector and loses badly at large
  $`N`$.** A PDL call has about 1.1–1.6 µs of fixed overhead. A
  pure-Perl rotation is a few dozen interpreted operations, so it
  should be comparable for one vector and 100× or more slower per
  vector at large $`N`$.
- <a id="h9"></a>**H9: in pure Perl, the matrix form beats the quaternion form once
  there are a few vectors.** The interpreter's cost is per operation,
  not per byte. The matrix form needs 9 multiplies and 6 adds per
  vector. The quaternion form needs about 18 multiplies and 12 adds,
  plus temporaries. Building the matrix (about 30 operations) is paid
  once per call. This is the opposite of the PP result in [H2](#h2), where
  memory traffic dominated.
- <a id="h10"></a>**H10: converting Perl data to ndarrays and back costs about as much
  as rotating it in Perl.** Both walk every element in the
  interpreter, so going through PDL should gain little or nothing for
  Perl data.

### Benchmark ([pureperl.pl](pureperl.pl))

The quaternion is a Perl list `(x, y, z, w)`. The Perl vector data is
either an array of `[x,y,z]` refs, as produced by `unpdl`, or three
column arrays.

| Label       | Implementation |
|-------------|----------------|
| `rot1 map`  | `map { rot1(\@q, $_) } @aoa`: one sub call per vector, as a scalar API would be used. Returns a new `[x,y,z]`. |
| `aoa`       | One call, loop inlined, quaternion form. Array of `[x,y,z]` in and out. |
| `aoa mat`   | As `aoa`, but builds the 3×3 matrix once and uses it in the loop. |
| `cols`      | Three column arrays in and out, quaternion form. |
| `cols mat`  | Three column arrays, matrix form. |
| `rotq`, `rotq3` | The PP kernels, with data already in ndarrays (as in earlier sections). |
| `AoA->rotq` | Round trip: `unpdl(rotq($q, pdl($aoa)))`. |
| `cols->rotq3` | Round trip: `pdl()` each column, `rotq3`, then `[ ->list ]` each result. |
| `pdl(AoA)`, `unpdl`, `pdl(cols)`, `cols list` | The conversions alone. |

- All pure-Perl forms agree with `rotq` to ≤ 2.2e-16.
- Same protocol as before: ≥ 0.5 s per cell, 3 runs, per-cell minimum.
  - Raw output: `results/pureperl-{1,2,3}.txt`.
  - Minima: [results/pureperl-min.txt](results/pureperl-min.txt).
- The machine's load average was about 1.3 during these runs.

### Results

All figures are ns per vector.

**(1) Pure Perl:**

|     N | rot1 map | aoa | aoa mat | cols | cols mat |
|------:|---------:|-----:|--------:|-----:|---------:|
| 1 | 1404 ns | 1371 ns | 2203 ns | 1737 ns | 2420 ns |
| 3 | 1119 ns | 948 ns | 1114 ns | 1199 ns | 1132 ns |
| 10 | 1071 ns | 807 ns | 658 ns | 924 ns | 736 ns |
| 30 | 1061 ns | 774 ns | 545 ns | 863 ns | 634 ns |
| 1e2 | 1135 ns | 803 ns | 520 ns | 815 ns | 512 ns |
| 1e3 | 1074 ns | 744 ns | 506 ns | 762 ns | 486 ns |
| 1e4 | 1077 ns | 746 ns | 499 ns | 748 ns | 469 ns |
| 1e5 | 1095 ns | 770 ns | 525 ns | 860 ns | 624 ns |
| 1e6 | 1119 ns | 819 ns | 569 ns | 798 ns | 691 ns |

**(2) PDL with data already in ndarrays, and the round trip from Perl
data:**

|     N | rotq | rotq3 | AoA->rotq | cols->rotq3 |
|------:|-----:|------:|----------:|------------:|
| 1 | 1789 ns | 3472 ns | 7873 ns | 16689 ns |
| 3 | 613 ns | 1129 ns | 3489 ns | 5495 ns |
| 10 | 199 ns | 386 ns | 2011 ns | 2060 ns |
| 30 | 69 ns | 136 ns | 1564 ns | 908 ns |
| 1e2 | 25 ns | 45 ns | 1427 ns | 510 ns |
| 1e3 | 8 ns | 10 ns | 1383 ns | 346 ns |
| 1e4 | 6 ns | 7 ns | 1483 ns | 350 ns |
| 1e5 | 6 ns | 6 ns | 1721 ns | 620 ns |
| 1e6 | 5 ns | 8 ns | 1881 ns | 1106 ns |

**(3) The conversions alone:**

|     N | pdl(AoA) | unpdl | pdl(cols) | cols list |
|------:|---------:|------:|----------:|----------:|
| 1 | 4665 ns | 561 ns | 8724 ns | 2671 ns |
| 3 | 2386 ns | 306 ns | 3079 ns | 1009 ns |
| 10 | 1530 ns | 222 ns | 1030 ns | 455 ns |
| 30 | 1213 ns | 197 ns | 446 ns | 293 ns |
| 1e2 | 1131 ns | 190 ns | 218 ns | 231 ns |
| 1e3 | 1127 ns | 228 ns | 123 ns | 207 ns |
| 1e4 | 1176 ns | 248 ns | 116 ns | 233 ns |
| 1e5 | 1268 ns | 422 ns | 113 ns | 443 ns |
| 1e6 | 1431 ns | 589 ns | 121 ns | 1016 ns |

### Analysis

**[H8](#h8): PDL wins from 3 vectors up, if the data is already in ndarrays
(confirmed).**
- For one vector, pure Perl is a little faster: 1371–1404 ns against
  1789 ns for `rotq`.
- At $`N = 3`$, PDL is already faster: 613 against 948 ns per vector.
  So the crossover is at about 2 vectors.
- From there PDL's lead grows as its fixed cost is spread over more
  vectors. It is 20× at $`N = 100`$ against the fastest pure-Perl form, and
  60–160× for $`N \ge 10^3`$ (5–8 ns against 470–820 ns).
- Pure Perl's cost per vector is flat, at about 0.5–0.8 µs inlined and
  about 1.1 µs with a sub call per vector. It rises slightly at
  $`10^5`$–$`10^6`$, presumably from allocating millions of Perl
  scalars and arrays.

**[H9](#h9): in pure Perl, the matrix form wins from about 10 vectors
(confirmed).**
- The matrix form is about 1.5× faster per vector than the quaternion
  form: 499–569 against 744–819 ns for arrays of `[x,y,z]`, and similar
  for columns.
- Building the matrix costs about 0.8 µs once per call. That makes the
  matrix form slower at $`N = 1`$ (2203 against 1371 ns), about even at
  $`N = 3`$, and faster from $`N = 10`$ on.
- This is the reverse of the PP result ([H2](#h2)). In the interpreter, the
  number of operations is what counts. In C, memory traffic dominates
  and both forms are about equally cheap.
- A sub call per vector (`rot1 map`) adds about 300 ns per vector
  (1071–1135 against 744–819 ns) for argument unpacking and the call
  itself. A bulk API should loop inside one call.
- Whether the data is an array of `[x,y,z]` or three column arrays
  makes little difference in pure Perl.

**[H10](#h10): converting Perl data for PDL is not worth it (confirmed, and
worse than expected for arrays of `[x,y,z]`).**
- **Arrays of `[x,y,z]`:** the round trip is slower than pure Perl at
  *every* $`N`$: 1.4–1.9 µs per vector, against 0.5–0.8 µs in Perl.
  - `pdl($aoa)` alone costs 1.1–1.4 µs per vector, more than the whole
    pure-Perl rotation. That is about 400 ns per element, from walking
    the nested Perl structure generically.
  - `unpdl` back costs another 0.2–0.6 µs.
- **Column arrays:** these convert about 10× faster. `pdl(\@x)` costs
  about 40 ns per element. Even so:
  - The round trip ties the best pure-Perl form (`cols mat`) at
    $`N = 100`$ and $`10^5`$. It only beats it at $`10^3`$–$`10^4`$, by
    1.3–1.4× (346 against 486 ns at $`10^3`$).
  - It loses at $`10^6`$ (1106 against 691 ns), where turning the
    results back into Perl lists (`->list`) costs about 1 µs per
    vector.
- **For one vector,** the round trip costs 7.9 µs, against 1.4 µs in
  pure Perl.

In every case, converting the data costs more than the PP rotation
itself: about 3× at $`N = 1`$, and over 100× at $`10^6`$. PDL's speed is only
available to callers whose data is *already* in ndarrays.

### Conclusions for the pure-Perl implementation

- <a id="c15"></a>**C15: Write the pure-Perl bulk rotation with the matrix form.** Build
  the 3×3 matrix once per call, then do 9 multiplies and 6 adds per
  vector inside a single loop. It is about 1.5× faster than the
  quaternion form for $`N \ge 10`$. This is the reverse of [C2](#c2), which applies
  to the PP kernels.
- <a id="c16"></a>**C16: Use the quaternion form directly for a single vector.** It avoids
  about 0.8 µs of matrix construction.
- <a id="c17"></a>**C17: Choose the implementation by the type of the input, not by
  whether PDL is installed.**
  - Perl arrays in, Perl arrays out: never go through PDL.
    Converting to ndarrays and back is slower than pure Perl for arrays
    of `[x,y,z]`, and at best 1.4× faster for columns. Use XS when it is
    available and pure Perl otherwise ([C21](#c21)).
  - ndarrays in: use the PP kernels.
  - Do not quietly "upgrade" Perl data to PDL inside the module.
- <a id="c18"></a>**C18: The scalar API (one quaternion, one vector as Perl values) should
  not use PDL, even when PDL is available.** In pure Perl it costs
  about 1.4 µs, against about 1.8 µs with data already in ndarrays, and
  7.9 µs if the values have to be converted to ndarrays first. XS is
  faster still ([C23](#c23)).
- <a id="c19"></a>**C19: Tell users who want PDL speed to keep their data in ndarrays.**
  Pure Perl is about 100× slower than the PP kernels at large $`N`$.
  If users must convert, three column arrays (`pdl(\@x)`) convert
  about 10× faster than an array of `[x,y,z]` refs (`pdl(\@aoa)`).

## Shared C implementation: XS on Perl data

### Question

The plan for `CXC::Math::Quaternion` includes one C implementation
of the quaternion math, shared by two front ends:
- PDL::PP kernels for ndarrays;
- plain XS functions for Perl data, which do not need PDL.

The pure-Perl code stays as the fallback when no compiler is
available. This section asks:

1. Does routing the PP kernel through a shared C function cost
   anything?
2. How much faster is XS on Perl data than pure Perl?
3. Where does XS stand against PDL?

### Hypotheses

- <a id="h11"></a>**H11: the shared C function costs nothing in the PP kernel.** The
  function is `static inline` in a header, so the compiler should
  produce the same code as the hand-written `rotq`.
- <a id="h12"></a>**H12: XS on Perl data is several times faster than pure Perl, but
  far slower than PP on ndarrays.** The arithmetic is a few ns per
  vector in C. What remains is reading Perl scalars (`av_fetch`,
  `SvNV`) and creating new ones (`newSVnv`, `newAV`). That should cost
  tens of ns per vector, against about 500 ns for pure Perl and about
  5 ns for PP.
- <a id="h13"></a>**H13: for a single vector, XS beats both pure Perl and PDL.** An XS
  call has far less fixed cost than a PDL call (about 1 µs, see
  [H8](#h8)), and does no interpreted arithmetic.

### Benchmark ([xs.pl](xs.pl))

The shared code is [qrot.h](qrot.h): one `static inline` function,
`qrot_vec(q, v, r)`, written in plain C with no Perl or PDL API. It
uses the same quaternion formula as `rotq`. [qbench.pd](qbench.pd)
includes it with `pp_addhdr`, so it is visible both to the PP kernels
and to the XS code added with `pp_addxs`:

| Function | Front end | Input → output |
|----------|-----------|----------------|
| `rotq_shared($q, $v)` | PP | Same signature as `rotq`. The kernel body copies `$q()`/`$v()` into local arrays and calls `qrot_vec`. |
| `xs_rot1(\@q, \@v)` | XS | One vector as an array ref → new `[x,y,z]` |
| `xs_rot1_list(@q, @v)` | XS | Seven scalars → list of three scalars |
| `xs_rot_aoa(\@q, \@aoa)` | XS | Array of `[x,y,z]` → new array of `[x,y,z]` |
| `xs_rot_cols(\@q, \@x, \@y, \@z)` | XS | Three column arrays → three new column arrays |

- The XS functions check their arguments (array refs, lengths) and
  `croak` on bad input, as a real module would.
- They read elements with `av_fetch` and `SvNV`, so they handle any
  numeric scalar.
- Here they share one shared object with the PP code, for convenience.
  In the real module the XS part would be built without PDL, from the
  same header.

The pure-Perl counterparts are the subs from the
[Pure-Perl baseline](#pure-perl-baseline), now in
[QRotPerl.pm](QRotPerl.pm), plus `rot1_list` (seven scalars in, three
out).
- Moving them into a module did not change their speed. An A/B test of
  the module sub against an inline copy, in one process, gave 818–843
  against 821–834 ns per vector.
- The machine was busier during these runs (load average about 2.1,
  against 1.3 for the pure-Perl section). That made pure Perl up to
  about 15% slower than in the earlier tables. To keep the comparison fair,
  pure Perl and XS are timed side by side in the same runs, and only
  those numbers are compared here.

The tables:
- **(1) One vector per call.** Each timed closure makes 1000 calls, so
  the closure's own cost is spread out. An empty loop of 1000
  iterations is the baseline: 14 ns per iteration.
- **(2) $`N`$ vectors per call.**
- **(3) `rotq_shared` against `rotq`,** with fresh and preallocated
  outputs.

All results agree exactly ($`\lvert \Delta\rvert = 0`$): XS against
pure Perl, and `rotq_shared` against `rotq`.

Same protocol as before: ≥ 0.5 s per cell, 3 runs, per-cell minimum.
- Raw output: `results/xs-{1,2,3}.txt`.
- Minima: [results/xs-min.txt](results/xs-min.txt).

### Results

**(1) One vector per call**, ns per call:

| Call | Time |
|------|---:|
| empty loop (baseline) | 14 ns |
| Perl `rot1(\@q, \@v)` → `[..]` | 979 ns |
| XS `xs_rot1(\@q, \@v)` → `[..]` | 288 ns |
| Perl `rot1_list(@q, @v)` → list | 877 ns |
| XS `xs_rot1_list(@q, @v)` → list | 178 ns |
| PDL `rotq($q, $v)`, ndarrays | 1487 ns |
| PDL `rotq($q, $v, $r)`, preallocated | 948 ns |

**(2) $`N`$ vectors per call**, ns per vector. "Perl" is the
quaternion form and "Perl mat" the matrix form ([C15](#c15)):

|     N | Perl aoa | Perl mat | XS aoa | Perl cols | Perl mat | XS cols | PDL rotq (ndarrays) |
|------:|---------:|---------:|-------:|----------:|---------:|--------:|--------------------:|
| 1 | 1431 ns | 2277 ns | 485 ns | 1818 ns | 2560 ns | 612 ns | 1706.8 ns |
| 3 | 1011 ns | 1075 ns | 246 ns | 1151 ns | 1178 ns | 263 ns | 689.5 ns |
| 10 | 951 ns | 756 ns | 186 ns | 1037 ns | 815 ns | 129 ns | 215.4 ns |
| 30 | 913 ns | 650 ns | 168 ns | 956 ns | 677 ns | 84 ns | 70.1 ns |
| 1e2 | 863 ns | 552 ns | 155 ns | 833 ns | 576 ns | 72 ns | 25.8 ns |
| 1e3 | 814 ns | 557 ns | 168 ns | 773 ns | 488 ns | 64 ns | 8.5 ns |
| 1e4 | 828 ns | 526 ns | 168 ns | 757 ns | 488 ns | 69 ns | 5.7 ns |
| 1e5 | 770 ns | 503 ns | 182 ns | 780 ns | 633 ns | 204 ns | 6.4 ns |
| 1e6 | 786 ns | 521 ns | 209 ns | 773 ns | 674 ns | 340 ns | 4.5 ns |

**(3) PP kernels**, ns per vector:

|     N | rotq | rotq_shared | rotq, prealloc | rotq_shared, prealloc |
|------:|-----:|------------:|---------------:|----------------------:|
| 1e2 | 25.3 ns | 24.4 ns | 17.1 ns | 16.7 ns |
| 1e4 | 5.7 ns | 5.7 ns | 5.2 ns | 5.2 ns |
| 1e6 | 4.4 ns | 4.5 ns | 3.3 ns | 3.6 ns |
| 1e7 | 16.9 ns | 17.4 ns | 3.3 ns | 3.3 ns |

### Analysis

**[H11](#h11): the shared C function costs nothing (confirmed).**
`rotq_shared` and `rotq` agree within the run-to-run noise at every
$`N`$, with fresh or preallocated outputs (for example 5.7 against
5.7 ns, and 3.3 against 3.3 ns). The compiler inlines `qrot_vec`, so
the copies into local arrays disappear. Writing the math once, in a
plain-C header, has no cost on the PDL side.

**[H12](#h12): XS on Perl data is 2–8× faster than the best pure Perl,
and up to 75× slower than PP on ndarrays (confirmed).**
- **Arrays of `[x,y,z]`:** XS costs 155–209 ns per vector for
  $`N \ge 10`$. That is 2.5–4× faster than the pure-Perl matrix form
  (503–756 ns), and 4–5.5× faster than the quaternion form.
- **Column arrays:** XS costs 64–84 ns per vector at $`N = 30`$ to
  $`10^4`$, 7–8× faster than the pure-Perl matrix form (488–677 ns).
- **Where XS spends its time.** The arithmetic is a small part: the PP
  kernel does the whole rotation, memory traffic included, in about
  5 ns per vector. The rest is Perl data handling:
  - each input number costs an `av_fetch` and an `SvNV`;
  - each output number costs a `newSVnv` and an `av_store`;
  - the array-of-`[x,y,z]` form also creates a new array and a
    reference for every output vector.

  That is why columns are 2–2.6× cheaper than arrays of `[x,y,z]`
  in XS (64–84 against 155–168 ns), while in pure Perl the two layouts
  cost about the same.
- **Large $`N`$.** At $`10^5`$–$`10^6`$, XS columns slow to 204–340 ns
  per vector, and their lead over pure Perl falls to 2–3×. Each call
  creates millions of new scalars and frees the previous call's
  millions. That costs about 250 ns per vector at $`10^6`$, even in a
  fresh process ([H15](#h15)). Pure Perl was already slow per element,
  so it degrades much less.
- **Against PDL:** PP on data already in ndarrays is 3–6× faster than
  XS at $`N = 100`$, and 7.5–75× faster for $`N \ge 10^3`$ (4.5–8.5
  against 64–340 ns), because it never touches a Perl scalar per
  element.
  - XS columns beat PDL up to about $`N = 20`$ (129 against 215 ns at
    $`N = 10`$; 84 against 70 ns at $`N = 30`$).
  - XS arrays of `[x,y,z]` beat PDL up to about $`N = 10`$.

**[H13](#h13): for one vector, XS beats both (confirmed).**
- `xs_rot1_list` takes 178 ns per call (about 165 ns over the empty
  loop). That is 4.9× faster than pure Perl's `rot1_list` (877 ns),
  and 5–8× faster than PDL (948–1487 ns).
- Returning a new `[x,y,z]` (`xs_rot1`, 288 ns) costs about 110 ns
  more than returning three scalars, for the array and the reference.
  It is still 3.4× faster than pure Perl's `rot1` (979 ns).

### Conclusions for the XS implementation

- <a id="c20"></a>**C20: Keep the quaternion math in one plain-C header, shared by the
  PP kernels and the XS functions.** It costs nothing on the PP side
  ([H11](#h11)). The header uses no Perl or PDL API, so the XS module
  can be built without PDL.
- <a id="c21"></a>**C21: For Perl data, use XS when it is available and pure Perl
  otherwise.** XS is 2–8× faster than the best pure-Perl form for
  many vectors, and 3.4–4.9× faster for one. Pure Perl is still needed
  as the fallback when no compiler is available. Neither path should
  convert to PDL ([C17](#c17)).
- <a id="c22"></a>**C22: The bulk XS API on Perl data should prefer column arrays over
  arrays of `[x,y,z]`.** Columns are 2–2.6× cheaper in XS (64–84
  against 155–168 ns per vector), because an array of `[x,y,z]` needs
  a new array and reference for every output vector. This matches the
  ndarray-side column support ([C7](#c7)). It applies to functions
  that return new arrays. When writing into existing arrays, the two
  layouts cost about the same ([C26](#c26)).
- <a id="c23"></a>**C23: The scalar API should take and return plain scalars.**
  `xs_rot1_list` costs 178 ns per call, against 288 ns for the
  array-ref form, about 1 µs for pure Perl, and 1–1.5 µs for PDL.
  Writing into an existing `[x,y,z]` (`xs_rot1_into`) is nearly as
  fast as the list form ([C27](#c27)).
- <a id="c24"></a>**C24: XS does not replace PDL for large data.** Every element of
  Perl data costs tens of ns in XS just to read and write. PP on
  ndarrays is 7.5–75× faster for $`N \ge 10^3`$. That reinforces
  [C19](#c19): users with large data should keep it in ndarrays.

## Writing into existing Perl arrays

### Question

The XS functions in the previous section return new arrays. For each
output vector they create new scalars, and for arrays of `[x,y,z]` also
a new inner array and a reference. The analysis there
([H12](#h12)) found that this Perl data handling, not the arithmetic,
is where XS spends its time. So:

1. How much faster is XS when it writes into arrays the caller already
   has?
2. Does pure Perl gain the same way?
3. Does writing in place, with the input array as the output, help
   further?

### Hypotheses

- <a id="h14"></a>**H14: writing into existing arrays removes most of XS's per-vector
  cost for arrays of `[x,y,z]`, less for columns, and little for pure
  Perl.**
  - The array-of-`[x,y,z]` form creates five things per vector (an
    array, a reference and three scalars), so it should gain most.
  - The column form creates three scalars per vector, so it should
    gain less.
  - Pure Perl's cost is mostly interpreting the arithmetic, so it
    should gain little.
- <a id="h15"></a>**H15 (added after the first run): XS's slowdown on large column
  arrays comes from where the scalars sit in memory, not from the
  `_into` code.** In the first run, writing into separate column
  arrays was as slow as creating new ones at $`N \ge 10^5`$ (about
  260 ns per vector), while writing in place stayed at about 45 ns.
  - Perl reuses freed scalars in an order that has nothing to do with
    array index.
  - The "new arrays" measurement just before it, in the same process,
    creates and frees millions of scalars on every call.
  - Arrays filled after that point have their scalars scattered, so
    walking them by index misses the cache.

### Benchmark ([xsinto.pl](xsinto.pl), [svlayout.pl](svlayout.pl))

New functions, all in [qbench.pd](qbench.pd) and
[QRotPerl.pm](QRotPerl.pm):

| Function | Writes into |
|----------|-------------|
| `xs_rot1_into(\@q, \@v, \@r)` | `@r` (three elements) |
| `xs_rot_aoa_into(\@q, \@vs, \@rs)` | `@rs`, an array of `[x,y,z]` |
| `xs_rot_cols_into(\@q, \@x, \@y, \@z, \@rx, \@ry, \@rz)` | `@rx`, `@ry`, `@rz` |
| `rot1_into`, `rot_aoa_mat_into`, `rot_cols_mat_into` | the same, in pure Perl (matrix form for the bulk versions, [C15](#c15)) |

- **Reuse.** The XS versions set existing scalars with `sv_setnv` and
  write into existing inner `[x,y,z]` arrays.
  - Missing elements, and elements that are not array refs, are
    created. So the first call into an empty array allocates, and later
    calls reuse.
  - Output arrays are not truncated: an output longer than the input
    keeps its extra elements.
  - A read-only output element raises Perl's usual "Modification of a
    read-only value" error.
- **In place.** The output may be the input itself. Each vector is read
  before it is written.
- **Pure Perl** assigns to existing elements (`$r->[0] = …`), which
  also reuses the scalars.
- **Correctness.** Every `_into` result agrees exactly with the
  corresponding new-array function, including in-place use.

[xsinto.pl](xsinto.pl) times four kinds of call for each layout:
- the new-array version;
- the `_into` version, into a separate output that already holds the
  previous call's results ("into");
- the `_into` version with the input as the output ("in place");
- for columns, `rotq3` on ndarrays with preallocated outputs, for
  reference.

Same protocol as before (3 runs, per-cell minimum, ≥ 0.5 s per cell):
raw `results/xsinto-{1,2,3}.txt`, minima in
[results/xsinto-min.txt](results/xsinto-min.txt).

These runs include an `av_extend(v, 2)` on every inner `[x,y,z]` array
in `xs_rot_aoa_into`. It guarantees room for three elements through
the documented API, rather than relying on how much Perl happens to
allocate. An A/B test against a build without it (four alternating
processes per build, $`N = 10^3`$ and $`10^4`$) gave 57–65 against
54–61 ns per vector, so its cost is within the noise (about 2 ns).

After these runs, `av_extend` was also added to `xs_rot1_into` (for its
three-element output) so that every function that writes into an array
reserves the room first; `xs_rot_aoa_into` now calls it on the outer
array unconditionally, once per call. Table (1) was measured before
that change. A separate A/B test of `xs_rot1_into` (four alternating
processes per build) gave 193–214 ns per call with the call and 181–209
without, so it adds about 10 ns (about 5%), barely above the noise.

[svlayout.pl](svlayout.pl) tests [H15](#h15) at $`N = 10^6`$ on
columns. Each run is a fresh process that builds the input and output
arrays and then times "into", "in place" and "new". The two modes
differ in one step:
- `clean`: the arrays are built straight away;
- `churn`: the process first creates and frees 3 × 3 million scalars,
  as a new-array benchmark does.

It was run three times in each mode, alternating:
[results/svlayout.txt](results/svlayout.txt).

### Results

**(1) One vector per call**, ns per call:

| Call | Time |
|------|---:|
| empty loop (baseline) | 15 ns |
| Perl `rot1` → new | 1021 ns |
| Perl `rot1_into` → into | 1015 ns |
| XS `xs_rot1` → new | 291 ns |
| XS `xs_rot1_into` → into | 187 ns |
| XS `xs_rot1_list` → list | 176 ns |

**(2) Array of `[x,y,z]`**, ns per vector:

|     N | Perl new | Perl into | XS new | XS into | XS in place |
|------:|---------:|----------:|-------:|--------:|------------:|
| 1 | 2157 ns | 2126 ns | 472 ns | 361 ns | 317 ns |
| 10 | 726 ns | 642 ns | 178 ns | 91 ns | 89 ns |
| 1e2 | 552 ns | 515 ns | 150 ns | 69 ns | 66 ns |
| 1e3 | 560 ns | 509 ns | 174 ns | 68 ns | 67 ns |
| 1e4 | 551 ns | 509 ns | 182 ns | 70 ns | 62 ns |
| 1e5 | 553 ns | 474 ns | 194 ns | 82 ns | 67 ns |
| 1e6 | 566 ns | 484 ns | 214 ns | 81 ns | 64 ns |

**(3) Three column arrays**, ns per vector:

|     N | Perl new | Perl into | XS new | XS into | XS in place | PDL `rotq3`, prealloc |
|------:|---------:|----------:|-------:|--------:|------------:|----------------------:|
| 1 | 2469 ns | 2241 ns | 582 ns | 321 ns | 312 ns | 1598.0 ns |
| 10 | 763 ns | 619 ns | 114 ns | 71 ns | 70 ns | 163.9 ns |
| 1e2 | 551 ns | 452 ns | 65 ns | 45 ns | 43 ns | 21.7 ns |
| 1e3 | 463 ns | 429 ns | 63 ns | 48 ns | 44 ns | 7.3 ns |
| 1e4 | 454 ns | 446 ns | 76 ns | 55 ns | 42 ns | 5.9 ns |
| 1e5 | 624 ns | 676 ns | 268 ns | 258 ns | 44 ns | 5.8 ns |
| 1e6 | 566 ns | 540 ns | 252 ns | 260 ns | 45 ns | 6.7 ns |

**(4) Scalar layout, $`N = 10^6`$ columns** ([svlayout.pl](svlayout.pl)),
ns per vector, range over 3 fresh processes per mode:

| Mode | XS into | XS in place | XS new |
|------|--------:|------------:|-------:|
| `clean` | 43–48 ns | 43–45 ns | 254–256 ns |
| `churn` | 72–81 ns | 54–59 ns | 245–258 ns |

### Analysis

**[H14](#h14): writing into existing arrays makes XS 2.2–2.6× faster for
arrays of `[x,y,z]`, about 1.4× for columns, and pure Perl only a
little (confirmed).**
- **Arrays of `[x,y,z]`:** XS goes from 150–182 to 68–70 ns per vector
  for $`N = 10^2`$–$`10^4`$. That is 2.2–2.6× faster, and about 7.5×
  faster than pure Perl writing into existing arrays (509–515 ns).
- **Columns:** XS goes from 63–76 to 45–55 ns per vector in the same
  range, 1.3–1.45× faster, and 8–10× faster than pure Perl writing
  into existing arrays.
- **Layout no longer matters much in XS.** With new arrays, columns
  were 2.3–2.8× cheaper than arrays of `[x,y,z]` in these runs
  ([C22](#c22)). When writing into existing arrays the gap shrinks to
  1.3–1.5× (45–55 against 68–70 ns), because the per-vector inner array and reference
  are no longer created.
- **In place** is slightly faster again (62–67 ns per vector for arrays
  of `[x,y,z]`, 42–44 ns for columns), since each scalar is read and
  written while it is in cache.
- **Pure Perl** gains little from writing into existing arrays: between
  −8% and +19% across the cells, mostly 5–15%. As expected, its cost
  is interpreting the arithmetic.
- **One vector:** `xs_rot1_into` (187 ns) is 1.6× faster than
  `xs_rot1` (291 ns), and about as fast as `xs_rot1_list` (176 ns).
  Pure Perl's `rot1_into` gains nothing (1015 against 1021 ns).
- **Against PDL:** writing into existing arrays narrows the gap to
  ndarrays, but does not close it. For columns, PP with preallocated
  outputs is still 6–9.5× faster for $`N \ge 10^3`$ (5.8–7.3 against
  42–55 ns). XS beats PDL somewhere between $`N = 10`$ and 100 (71
  against 164 ns at $`N = 10`$; 45 against 22 ns at $`N = 100`$).

**[H15](#h15): the large-$`N`$ slowdown for separate output columns is a
layout effect (confirmed).**
- In a fresh process ([svlayout.pl](svlayout.pl) `clean`), writing
  into separate output columns at $`N = 10^6`$ costs 43–48 ns per
  vector. That is the same as in cache, and nowhere near the 258–260 ns
  measured in [xsinto.pl](xsinto.pl).
- After freeing 9 million scalars first (`churn`), the same call costs
  72–81 ns, and in place 54–59 ns. Arrays built from reused, scattered
  scalars are slower to walk.
- In [xsinto.pl](xsinto.pl), the output columns are filled right after
  the "XS new" cell has run for half a second, creating and freeing
  millions of scalars per call. That scatters them far more than
  `churn` does, which explains the 258–260 ns.
- The new-array functions cost about 250 ns per vector at $`10^6`$ in
  both modes. That cost comes from creating 3 million scalars and
  freeing the previous 3 million on every call, whatever the layout.
- So the "XS into" figures for columns at $`N \ge 10^5`$ in table (3)
  overstate the intrinsic cost. Measured in a fresh process, writing
  into existing arrays is **5–6× faster than creating new ones** at
  $`10^6`$ (43–48 against about 250 ns). Arrays of `[x,y,z]` are less
  affected (81–82 ns into, 64–67 in place), presumably because each
  inner array's three scalars were created together.

### Conclusions for writing into existing arrays

- <a id="c25"></a>**C25: Offer `_into` variants in the XS API for bulk Perl data.** The
  caller passes output arrays, which are filled on first use and
  reused afterwards. In-place use (output = input) should be allowed.
  - Arrays of `[x,y,z]`: 2.2–2.6× faster than returning new arrays.
  - Columns: 1.3–1.45× faster in cache, and 5–6× faster at $`10^6`$
    measured in a fresh process.
- <a id="c26"></a>**C26: With `_into`, both Perl layouts are acceptable.** Arrays of
  `[x,y,z]` and column arrays are within 1.3–1.5× of each other, so
  the preference for columns ([C22](#c22)) only matters for functions
  that return new arrays.
- <a id="c27"></a>**C27: For one vector, `xs_rot1_into` and the list form are about
  equally fast (187 and 176 ns).** Either suits the scalar API ([C23](#c23)).
  Choose for convenience: the list form for one-off calls, `_into`
  when the caller keeps a result array.
- <a id="c28"></a>**C28: Pure Perl does not need `_into` for speed (mostly a 5–15% gain).** It
  is still worth providing for API symmetry with XS, so that code
  works the same with either backend ([C21](#c21)).

## Comparison with PDL::LinearAlgebra

### Question

PDL::LinearAlgebra wraps BLAS and LAPACK. Its `mmult` and `gemm`
multiply matrices with an optimised BLAS library. Rotating vectors by
a rotation matrix is a matrix multiply, so:

1. Is BLAS faster than the PP quaternion kernel for one quaternion and
   many vectors?
2. Can it handle many quaternions, one per vector?
3. How does it compare with PDL's own matrix multiply?

### Hypotheses

- <a id="h16"></a>**H16: BLAS is no faster than the PP kernel for one quaternion.**
  Multiplying a 3×3 matrix by a 3×N matrix has an inner dimension of
  3, so it does very little arithmetic per number read. Like `rotq`,
  it is limited by memory traffic ([H3](#h3)), and BLAS is tuned for
  large matrices, where blocking and vectorising the inner loop pay
  off. At best it should match `rotq`.
- <a id="h17"></a>**H17: BLAS is much slower for many quaternions, one per vector.**
  That needs one 3×3 matrix per vector. `mmult` broadcasts over them
  by calling BLAS once per vector, so every 3×3 product pays the cost
  of a BLAS call. The matrices also have to be built from the
  quaternions first.

### Benchmark ([lalg.pl](lalg.pl), [autopthread.pl](autopthread.pl))

[lalg.pl](lalg.pl) times, in ns per vector:

| Label | Call |
|-------|------|
| `rotq`, `rotq pre` | the PP kernel, with a fresh or preallocated output |
| `matmult` | PDL's own matrix multiply, `PDL::matmult($v, $Rt)`, as `x` in [bench.pl](bench.pl) |
| `x (=mmult)` | `$v x $Rt` with PDL::LinearAlgebra loaded (see below) |
| `mmult` | `mmult($v, $Rt)`: BLAS `dgemm`, fresh output |
| `gemm pre` | `gemm($R, 1, 0, $v, 1, 0, $c)`: BLAS `dgemm` into a preallocated `$c` |
| `mmult (mats pre)` | many quaternions: `mmult($v->dummy(1), $Rs)`, broadcasting over (3,3,N) matrices built beforehand |
| `build + mmult` | the same, including building the N matrices from the quaternions with PDL operations |

- `gemm` follows BLAS's Fortran convention, in which a PDL `(m,n)`
  ndarray is an m-row matrix. So the rotation is written as
  $`C = R_F^{\mathsf T} V`$, with the transpose flag on `$R`.
- All methods agree with `rotq` to ≤ 2.8e-16.
- OpenBLAS starts its own thread pool, separate from PDL's automatic
  threading. [lalg.pl](lalg.pl) was therefore run 3 times with
  OpenBLAS's default threading and 3 times with
  `OPENBLAS_NUM_THREADS=1`, alternating. The tables give the per-cell
  minimum:
  - [results/lalg-default-min.txt](results/lalg-default-min.txt);
  - [results/lalg-1thread-min.txt](results/lalg-1thread-min.txt).
- The machine was busy during these runs (load average 2.5–3.7).

**Loading PDL::LinearAlgebra changes `x` for the whole program.**
Loading `PDL::LinearAlgebra::Real`, which `PDL::LinearAlgebra` loads,
replaces the overload of PDL's `x` operator: it goes from
`PDL::overload_matmult` to the module's BLAS-based
`PDL::overload_mmult`. In separate processes at $`N = 10^5`$, `x` took
18.9 ns per vector without the module and 4.9–5.4 ns with it. The
`PDL::matmult` function itself is not replaced (it is still defined
in PDL's `Primitive.pd`), so the `matmult` column measures PDL's own
code. The `x` column in [bench.pl](bench.pl) is likewise PDL's own,
since that script does not load PDL::LinearAlgebra.

[autopthread.pl](autopthread.pl) follows up a result from
[lalg.pl](lalg.pl). It times `rotq` with a preallocated output, with
PDL's default threading threshold (1 M elements) and with the
threshold set to 0 (always split across threads), back to back in the
same process. 3 runs, per-cell minimum:
[results/autopthread-min.txt](results/autopthread-min.txt).

### Results

**(1) One quaternion, $`N`$ vectors**, ns per vector. Each cell gives
OpenBLAS default threading, then `OPENBLAS_NUM_THREADS=1` in
parentheses. The PP and PDL columns don't use BLAS, so their two
values differ only by noise:

|     N | rotq | rotq pre | matmult | x (=mmult) | mmult | gemm pre |
|------:|-----:|---------:|--------:|-----------:|------:|---------:|
| 1 | 2150 ns (1607 ns) | 1183 ns (1086 ns) | 3375 ns (3303 ns) | 4380 ns (4161 ns) | 1893 ns (1680 ns) | 4873 ns (4636 ns) |
| 1e2 | 26.6 ns (25.0 ns) | 17.8 ns (16.0 ns) | 56.8 ns (55.1 ns) | 54.4 ns (53.5 ns) | 28.7 ns (27.6 ns) | 55.5 ns (50.6 ns) |
| 1e3 | 7.9 ns (7.7 ns) | 6.7 ns (6.4 ns) | 22.3 ns (21.4 ns) | 10.8 ns (10.6 ns) | 8.1 ns (7.7 ns) | 10.9 ns (10.2 ns) |
| 1e4 | 6.2 ns (6.1 ns) | 5.4 ns (5.5 ns) | 17.8 ns (17.9 ns) | 6.6 ns (6.4 ns) | 6.4 ns (6.1 ns) | 6.1 ns (6.0 ns) |
| 1e5 | 6.0 ns (6.0 ns) | 5.4 ns (5.4 ns) | 17.6 ns (18.4 ns) | 4.9 ns (6.2 ns) | 4.8 ns (6.2 ns) | 2.8 ns (5.8 ns) |
| 1e6 | 4.6 ns (4.5 ns) | 3.3 ns (3.3 ns) | 19.0 ns (19.1 ns) | 8.8 ns (7.8 ns) | 10.3 ns (7.8 ns) | 7.0 ns (7.0 ns) |
| 1e7 | 17.1 ns (16.9 ns) | 3.3 ns (3.3 ns) | 32.0 ns (32.0 ns) | 28.7 ns (21.4 ns) | 29.3 ns (21.3 ns) | 8.0 ns (7.1 ns) |

**(2) $`N`$ quaternions, one vector each**, ns per vector (OpenBLAS
default threading, then one thread):

|     N | rotq | mmult (mats pre) | build + mmult |
|------:|-----:|-----------------:|--------------:|
| 1e2 | 26.6 ns (25.1 ns) | 272 ns (268 ns) | 2233 ns (2235 ns) |
| 1e4 | 6.2 ns (6.2 ns) | 206 ns (196 ns) | 408 ns (407 ns) |
| 1e6 | 6.0 ns (6.0 ns) | 827 ns (740 ns) | 1176 ns (1091 ns) |

**(3) `rotq` with a preallocated output, PDL threading threshold**
([autopthread.pl](autopthread.pl)), ns per vector:

|     N | default threshold (1 M elements) | threshold 0 |
|------:|---------------------------------:|------------:|
| 1e3 | 7.5 ns | 228.5 ns |
| 1e4 | 6.0 ns | 25.8 ns |
| 1e5 | 6.5 ns | 4.1 ns |
| 3e5 | 6.5 ns | 3.5 ns |

### Analysis

**[H16](#h16): BLAS matches the PP kernel at best, and loses at large
$`N`$ (confirmed, with one exception that comes from threading).**
- **Up to $`10^4`$ vectors,** `mmult` and `rotq` cost about the same
  with fresh outputs (for example 6.4 against 6.2 ns at $`10^4`$).
  `gemm` into a preallocated output is slower than `rotq` into one at
  small $`N`$ (55.5 against 17.8 ns at $`N = 100`$; 4.9 against 1.2 µs
  for one vector). It has more per-call cost, presumably in its wrapper;
  that was not investigated.
- **At $`N = 10^5`$, OpenBLAS with its default threading wins:** `gemm`
  2.8 against 5.4 ns for `rotq`, both with preallocated outputs. With
  `OPENBLAS_NUM_THREADS=1` it is 5.8 ns, no better than `rotq`. So the
  win comes from threads, not from a faster kernel. OpenBLAS splits
  this size across threads; PDL does not, because $`10^5`$ vectors are
  $`3 \times 10^5`$ elements, below PDL's 1 M-element threshold
  ([C6](#c6)).
- **At $`10^6`$–$`10^7`$, where PDL threads too,** `rotq` with a
  preallocated output is 2.1–2.4× faster than `gemm` with one (3.3
  against 7.0–8.0 ns). `mmult` with fresh outputs is slower still
  (7.8–29 ns), because it also pays the fresh-output page faults
  ([H3](#h3)). It is worse with OpenBLAS threading than without at
  $`10^7`$ (29 against 21 ns).
- **BLAS against PDL's own `matmult`:** `mmult` is 2.8–3.7× faster than
  `PDL::matmult` for $`N = 10^4`$–$`10^5`$ (4.8–6.4 against
  17.6–18.4 ns), and about 2× faster at $`10^6`$. So for code that does
  use matrices, loading PDL::LinearAlgebra speeds up `x`.
- **The `x` operator after loading the module** costs the same as
  `mmult` for $`N \ge 10^4`$. At small $`N`$ it adds Perl-level
  overhead (54.4 against 28.7 ns at $`N = 100`$; 4.4 against 1.9 µs for
  one vector).

**[H17](#h17): one matrix per vector is far slower with BLAS (confirmed).**
- Broadcasting `mmult` over $`N`$ precomputed 3×3 matrices costs
  196–827 ns per vector, 10–140× slower than `rotq` with $`N`$
  quaternions (6.0–26.6 ns). Each vector is a separate BLAS call, and
  the call costs far more than a 3×3 product.
- Building the matrices from the quaternions with ordinary PDL
  operations adds about 200–2000 ns per vector more.
- `rotq` broadcasts over quaternions with no extra cost per vector.

**PDL's threading threshold** ([autopthread.pl](autopthread.pl)):
- With the threshold at 0, `rotq` is 1.6× faster at $`10^5`$ vectors
  (4.1 against 6.5 ns) and 1.9× faster at $`3 \times 10^5`$ (3.5
  against 6.5 ns). That is close to multi-threaded OpenBLAS at
  $`10^5`$ (2.8 ns).
- Below that, starting threads costs more than the work: 228.5 against
  7.5 ns per vector at $`10^3`$, and 25.8 against 6.0 at $`10^4`$.
- So on this machine, for this kernel, a threshold somewhere between
  $`10^4`$ and $`10^5`$ vectors ($`3 \times 10^4`$–$`3 \times 10^5`$
  elements) would be better than PDL's default 1 M elements. It
  should not be 0.

### Conclusions for PDL::LinearAlgebra

- <a id="c29"></a>**C29: Don't use PDL::LinearAlgebra for rotations.** For one
  quaternion it is never faster than the PP kernel on the same number
  of threads, and it is 2.1–2.4× slower at $`10^6`$–$`10^7`$ vectors.
  For one quaternion per vector it is 10–140× slower. It would also
  add a dependency on BLAS, LAPACK and a Fortran toolchain. This
  supports [C1](#c1) and [C2](#c2).
- <a id="c30"></a>**C30: Don't load PDL::LinearAlgebra from `CXC::Math::Quaternion`,
  and don't rely on what `x` does.** Loading it replaces PDL's `x`
  operator for the whole program: about 3× faster for large products, but
  2× slower than `mmult` for small ones. A library should not change
  that for its users. Any matrix products the module needs should
  call a specific function, not `x`.
- <a id="c31"></a>**C31: Document PDL's threading threshold for mid-sized inputs.** At
  $`10^5`$–$`3 \times 10^5`$ vectors, PDL runs `rotq` on one thread.
  Lowering the threshold (`PDL::set_autopthread_size`) makes it
  1.6–1.9× faster, but setting it to 0 makes small calls 4–30× slower.
  Like the malloc setting in [C4](#c4), this is a setting for the
  application, not something the module should change.

## Argument dispatch

### Question

The draft API for `CXC::Math::Quaternion` has one function per
operation. The function works out from its arguments which data form
it has been given (a single value, records, columns or an ndarray) and
which output mode (allocate, into a caller's container, or in place).
It then checks shapes and calls the matching kernel. With small inputs,
that dispatch could cost more than the work itself. How much does it
cost, and where should it be written?

### Hypotheses

- <a id="h18"></a>**H18: dispatch written in Perl costs more than the work for a single
  vector.** A Perl sub call alone is a few hundred ns. Classifying
  three arguments, checking them and detecting in place adds several
  more operations at 0.1–0.3 µs each, so it should total at least
  1 µs. That is several times the 0.2–0.3 µs that XS takes to rotate
  one vector ([H13](#h13)).
- <a id="h19"></a>**H19: dispatch written in XS costs tens of ns.** The same checks in C
  (`SvROK`, `SvTYPE`, looking at the first element) are a few ns each.

### Benchmark ([dispatch.pl](dispatch.pl), [dispatch-parts.pl](dispatch-parts.pl))

Two prototypes implement `rotate_vec($q, $v, [$out])` for every data
form and output mode, with output truncation:

- **Perl dispatcher, `QDispatch::rotate_vec`** ([QDispatch.pm](QDispatch.pm)).
  It classifies arguments with `reftype`, `blessed` and a look at the
  first element, checks the output's kind, and calls the existing
  kernels: the XS functions for Perl data, and `rotq` for ndarrays.
  - For ndarrays it checks that the output is `double`.
  - It converts a quaternion given as Perl data with `PDL->pdl`.
- **XS dispatcher, `xsd_rotate_vec`** (in [qbench.pd](qbench.pd)). It
  does the same in C for Perl data, with its own loops (which also
  shorten each output record to three elements). Blessed arguments,
  ndarrays or objects, are passed to the Perl dispatcher with `call_pv`.

Both look only at the top level and the first element. Checking every
record is left to the kernel's loop, which reads every element anyway.
That keeps dispatch a fixed cost per call. In both prototypes the
quaternion must be a single value.

[dispatch-test.pl](dispatch-test.pl) checks both dispatchers against
the kernels; all 33 checks pass:
- every form and output mode, with truncation;
- in place, including an ndarray passed as both input and output to
  `rotq`;
- the errors: a `float` ndarray output, the wrong kind of output
  container, a short record, and an argument that is not a container.

[dispatch.pl](dispatch.pl) times each case three ways: the matching
kernel called directly, through the Perl dispatcher, and through the
XS dispatcher. Small cases repeat the call 100–1000 times inside each
timed closure. [dispatch-parts.pl](dispatch-parts.pl) times the
individual steps of the Perl dispatcher. Both were run 3 times;
minima are in [results/dispatch-min.txt](results/dispatch-min.txt) and
[results/dispatch-parts-min.txt](results/dispatch-parts-min.txt).
- The load average was about 3 during these runs. In the first round,
  the same direct `rotq` call measured 1.1 µs in one row and 1.9 µs in
  another. These are the second round's figures.
- The overheads below are the differences between the minimum times,
  in ns per call.

### Results

**(1) Perl data.** The first column is the time of the kernel called
directly, with no dispatcher. The other two are the *extra* time each
dispatcher adds to that:

| Case | Kernel called directly | Added by Perl dispatch | Added by XS dispatch |
|------|--------------:|--------------:|------------:|
| single vector, new | 324 ns | +1124 ns | +12 ns |
| single vector, into | 215 ns | +1276 ns | +58 ns |
| single vector, in place | 227 ns | +1355 ns | +68 ns |
| 10 records, new | 1624 ns | +1379 ns | +177 ns |
| 10 records, into | 804 ns | +1460 ns | +153 ns |
| 10 columns, into | 623 ns | +3626 ns | +469 ns |

**(2) ndarrays**, laid out as table (1). The XS dispatcher passes
ndarrays to the Perl one, so both include the Perl dispatcher's cost:

| Case | `rotq` called directly | Added by Perl dispatch | Added by XS dispatch |
|------|--------------:|--------------:|------------:|
| N=1, quaternion as an ndarray | 1.10 µs | +3.18 µs | +3.68 µs |
| N=1, quaternion as Perl data | 1.08 µs | +8.07 µs | +8.82 µs |
| N=1000, quaternion as an ndarray | 7.23 µs | +3.04 µs | +4.14 µs |

**(3) Steps of the Perl dispatcher**, ns each, including about 25 ns
of loop overhead:

| Step | Time |
|------|---:|
| call a Perl sub with three arguments | 290 ns |
| `reftype` and `blessed` on an array ref | 170 ns |
| records or single value? (`ref $x->[0]`) | 107 ns |
| `refaddr($out) == refaddr($in)` | 133 ns |
| hash slice `@{$h}{qw(x y z)}` | 300 ns |
| `blessed` and `->isa('PDL')` | 323 ns |
| `$out->isnull` | 169 ns |
| `$out->type == PDL::double()` | 1620 ns |
| `$out->get_datatype == $double` (integer) | 204 ns |
| `PDL->pdl` of four numbers | 4432 ns |

With 1000 records or columns, the differences come mostly from the
dispatchers' kernel loops and from memory layout ([H15](#h15)), not
from dispatch. For example, the XS dispatcher's loop also shortens
each record. So they are left out of table (1).

### Analysis

**[H18](#h18): Perl dispatch costs several times the work for a single vector
(confirmed).**
- The Perl dispatcher adds 1.1–1.5 µs per call on single values and
  records, and 3.6 µs on columns. That is 4–6× the time XS takes to
  rotate one vector (215–324 ns).
- Table (3) accounts for it: about 0.3 µs for the sub call, plus 0.1–0.3
  µs for each classification, `refaddr` or hash access. Columns cost
  more because of the hash slices and the per-column checks.
- The overhead stays under 10% of the work only from about 200
  vectors: XS writing into existing records costs about 60–70 ns per
  vector ([C25](#c25)).

**[H19](#h19): XS dispatch costs tens of ns for single values (confirmed),
and a few hundred for columns.**
- Single values: 12–68 ns, i.e. 5–30% of the kernel.
- 10 records: 153–177 ns, about 10–20%. This includes shortening each
  output record, which the direct kernel does not do.
- 10 columns: 469 ns. Most of that is six hash lookups (three keys in,
  three out). Precomputed key hashes might reduce it; that was not
  measured.

**ndarrays.** The Perl-side checks cost about 3 µs, nearly three times
PDL's own fixed cost per call (about 1.1 µs). Two steps dominate:
- **The type check:** comparing PDL Type objects
  (`$out->type == PDL::double()`) takes 1.6 µs. Comparing the integer
  type code (`get_datatype`) takes 0.2 µs.
- **Converting a quaternion given as Perl data** with `PDL->pdl` takes
  4.4 µs on every call.

With the integer check, the ndarray path should cost about 1.5 µs. I
estimated that from table (3); it was not measured. Without the
`pdl` conversion (for example, passing a single Perl quaternion to the
PP kernel as four scalar parameters) it should cost about the same
whether the quaternion is Perl data or an ndarray. From about 1000
vectors (PP costs about 5–7 ns each), even 3 µs is under half the
call, and it shrinks from there.

### Conclusions for argument dispatch

- <a id="c32"></a>**C32: Write the dispatch for Perl data in XS.** In Perl it adds
  1.1–1.5 µs per call: 4–6× the work for a single vector, and more than
  10% until about 200 vectors. In XS it adds 12–68 ns for a single
  vector and about 150–180 ns for 10 records. So one entry point per
  operation, as in the draft `CXC::Math::Quaternion` API, is
  affordable, provided that entry point is an XSUB.
- <a id="c33"></a>**C33: The pure-Perl fallback can dispatch in Perl.** Its own work is
  about 1 µs per vector ([H8](#h8)). Perl dispatch therefore roughly
  doubles a single-vector call, but is under 10% from about 25–30
  vectors. The fallback only runs where no compiler is available.
- <a id="c34"></a>**C34: Keep the ndarray path's checks cheap.**
  - Compare integer type codes, not PDL Type objects (0.2 against
    1.6 µs).
  - Don't convert a Perl quaternion with `pdl` on every call (4.4 µs).
    Pass it to the kernel as scalar parameters instead, or have the
    XSUB recognise ndarrays and call the PP kernel's C function
    directly.

  With those changes the estimated overhead, about 1.5 µs, is close to
  PDL's own per-call cost, and negligible from about 1000 vectors.

## Caveats and open questions

- These results come from a single machine: a laptop CPU inside a
  container, with frequency scaling. Absolute numbers will differ on
  servers. The ratios and the mmap behaviour, which is glibc-specific,
  should carry over to Linux/glibc. Other platforms (macOS, musl,
  Windows) have different allocators.
- Only `double` was benchmarked.
- BLAS results depend on the library and its build. Only OpenBLAS
  0.3.26 (pthreads) was measured. Reference BLAS, MKL, BLIS or another
  OpenBLAS build could differ, especially in how they thread small
  products.
- The column-layout runs used only the one-quaternion case. Quaternions
  supplied as four columns, `(N)` each, were not measured. They would
  need a `q0(); q1(); q2(); q3()` signature variant.
- The pure-Perl baseline used only the one-quaternion case, and only
  arrays of `[x,y,z]` and three column arrays as Perl data.
  - A flat Perl array `(x0,y0,z0,x1,…)` was not measured.
  - Faster ways to build an ndarray from Perl data, such as `pack 'd*'`
    into the data of a `new_from_specification` ndarray, were not
    measured.
- Flat Perl arrays, and `pack`ed strings of doubles as XS input, were
  not measured.
- Large-$`N`$ timings on Perl data depend on where Perl's allocator put
  the scalars, which depends on what the process allocated and freed
  before ([H15](#h15)). The XS figures for Perl data at
  $`N \ge 10^5`$ should be read with that in mind, and future
  benchmarks of Perl data at that size should run each case in a fresh
  process.
- Here the XS functions are linked into the same shared object as the
  PP code. A separate XS module without PDL should perform the same,
  but that was not measured.
- Quaternion–quaternion products, normalisation and the
  RA/Dec/Roll ↔ quaternion conversions were not benchmarked. They are
  expected to behave like `rotq`: one pass, no temporaries.
