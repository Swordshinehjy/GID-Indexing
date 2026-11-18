# GIDIndexing — GIWAXS peak indexing tool (C++ / wxWidgets)

Starting from a grazing-incidence diffraction peak position table `(q_xy, q_z)`, it derives the unit cell parameters and completes indexing (hkl assignment).
This project is a next-generation port of the GIDInd family, offering three indexing modes, with the numerical kernel fully separated from the wxWidgets UI.

The artifact is a **single statically linked exe** (`/MT` + wxWidgets 3.2 static libraries, no DLL dependencies).

## Three Indexing Modes

| Mode | Prerequisite | Algorithm summary |
| --- | --- | --- |
| **With q_spec (two-step method)** | Specular rod peaks exist in the data (q_xy≈0) | The first step derives the in-plane metric from q_xy only; the second step fixes the in-plane parameters and solves for the out-of-plane component, then applies Buerger/Niggli reduction and numerical refinement. This mode has five implementations: improved / **default PeriodVote+DDL** / RANSAC / integer least squares / **SCRM (w-eliminated two-step method)** |
| **Without q_spec (single-step method)** | No specular peaks, q_spec unknown | The unknowns are the 6 components of the reciprocal metric `G*` + the out-of-plane coefficient `ν`; seed `(h,k,l)` 3×3 solve → beam search extending peak by peak (largest `\|ν\|` component as solve axis) → linear metric solve ((A) residual layer-by-layer pruning) → reassignment/refinement |
| **Fixed contact plane (uvw)** | Contact plane known (or assumed) | Searches the unit cell and each peak's `(h,k,l)` simultaneously. Main path **FP-DL**: three-peak integer Cramer solve for `ν` → linear solve for `G*` (`G tᵀ=(t·ν)ν` + `q²=xGxᵀ`) → alternating "assignment + linear re-solve" → a small amount of LM polishing; there is also an initial-value fast lane and a rod fundamental two-step method. |

## Build

Dependencies:

- MSVC (VS 2022, x64)
- wxWidgets 3.2 static build (`D:\wxWidgets\install`, located via CMake CONFIG mode)
- **xlsxio** (vcpkg static libraries, `D:\vcpkg\installed\x64-windows-static`): reads `.xlsx / .xlsm`,
  its dependencies minizip + expat must be installed as well, and `zlib.lib` / `bz2.lib` must be provided

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

Artifacts:

- `build\Release\GIDIndexing.exe` — main program (single-file static link)
- `build\Release\GIDSelfTest.exe` — core algorithm self-test (optional, `GID_BUILD_SELFTEST=ON`)
- `build\Release\GIDDipProbe.exe` — forward-model/diagnostics utility
- `build\Release\GIDNoSpec.exe` — no-q_spec algorithm probe (DIP 11 peaks + manual solution comparison)
- `build\Release\GIDSpecBench.exe` — timing benchmark of the two improved with-q_spec implementations (improved / PeriodVote+DDL)
  (`--l/--w/--R/--uv/--hk1/--max1/--refine/--no-early` override parameters, `--micro` prints the fine-grained cost of the refinement stage)
- `build\Release\GIDScrmProbe.exe` — SCRM self-check probe (synthetic DIP / synthetic (001) / parameter sweep, including the second-step funnel)
- `build\Release\GIDScrmFile.exe <peak-table.csv> [q_spec]` — runs SCRM on a real peak table and compares it with the improved version;
  when `a b ga` (the 4th–6th parameters) are given, it prints the per-`seedL` solution counts of the "in-plane candidate closest to the ground truth"

The third-party paths are adjusted in `VCPKG_STATIC_DIR` / `CMAKE_PREFIX_PATH` at the top of `CMakeLists.txt`.


## Directory Structure

```
src/
  core/                  numerical kernel (pure computation, no UI dependency)
    types.h              unified input/output: Peak / Solution / Progress
    cell.{h,cpp}         unit cell <-> reciprocal lattice, forward model, crystal-system constraints, Niggli / Buerger reduction
    linalg.{h,cpp}       small dense linear algebra + MATLAB-style array utilities (solve/QR/sortrows/unique…)
    refine.{h,cpp}       post-indexing processing (hkl assignment, residual statistics, LM refinement)
    parallel.h           parallelFor (multi-threaded parallelism)
    specindex.{h,cpp}    with-q_spec two-step method orchestration and dispatch: Scrm → scrm::runSpecIndex,
                         all other Part1Algorithm values → improved::runSpecIndex
                         (the original GIDInd part1/part2 porting paths have been removed)
    specindex_improved.cpp improved two-step method (from index2_other: Part ONE + Part TWO + refinement,
                         plus the PeriodVote/DDL, RANSAC and integer least squares Part TWO variants)
    scrm.{h,cpp}         SCRM w-eliminated two-step method: the first step enumerates (u,v) only to solve the in-plane metric,
                         the second step solves a 3×3 system for (κ,λ,ρ) and analytically recovers w and the per-peak l
    nospec.{h,cpp}       no-q_spec single-step method
    fixedplane.{h,cpp}   fixed contact plane (uvw)
    engine.{h,cpp}       unified entry runIndexing, automatic q_spec detection, hkl alignment of solutions
    dataio.{h,cpp}       .csv / .xlsx reading, export of indexed solutions to CSV
  ui/
    MainFrame.{h,cpp}    main window: peak table + scatter plot + three parameter pages + candidate solution table
    PeakPlot.{h,cpp}     overlay of measured (blue dots) and predicted (red crosses) peak positions
  App.cpp                wxApp entry point
app.rc / app.ico         window icon and resources (toolbar icons are drawn in code)
tools/
  SelfTest.cpp           synthetic-data inversion self-test (forward -> index -> verify)
  DipProbe.cpp           forward-model/residual probe
  NoSpecProbe.cpp        no-q_spec probe (DIP 11 peaks, specular peaks removed, compared with the manual solution)
  ScrmProbe.cpp          SCRM probe (synthetic DIP / synthetic (001) degenerate case / parameter sweep of seedL and R)
  ScrmFile.cpp           runs SCRM on a real peak table CSV and compares it with the improved version (ideal in-plane candidate → second step)
  verify_core_sync.ps1   verifies that the shared kernel files of this project and ../analyzer are byte-for-byte identical
```

## Usage Workflow

1. **Import a peak table** (`Ctrl+O` / toolbar "Import"):
   supports `.csv / .txt / .dat / .xlsx / .xlsm`; at least two columns are required, and if the column count > 2 a dialog pops up to select
   the columns holding `q_xy` and `q_z`. Units are Å⁻¹.
   **The sign of q_xy** only reflects the mirror handedness of the detection geometry; on import `|q_xy|` is always taken (computation, display,
   the delta column and export all use this normalized data; the core entry point is `normalizeQxySigns`).
2. **q_spec**: detected automatically on import (the smallest non-zero q_z among the specular rod peaks with q_xy≈0),
   it can also be entered manually, or re-obtained by clicking the "Auto-detect" button on the parameter page.
3. **Part 1 algorithm** (the "Algorithm" subsection of the with-q_spec parameter page):
   - `PeriodVote/DDL (this project)`: **default**. Part ONE is the same as the improved version; Part TWO instead uses **difference vectors + period voting**
     to determine `ρ` and `(κ,λ)` directly, no longer enumerating `l` triples, so the runtime is essentially independent of the `l`/`w` ranges (see below).
     If a whole round yields no solution, it automatically falls back to the improved version's enumeration method and runs again.
   - `Improved (this project)`: the starting three reflections are taken simultaneously from both seed sets "the first n peaks with the smallest \|q\|" (same as the original)
     and "the first n peaks with the smallest q_xy"; candidates are processed in ascending order of in-plane parallelogram area;
     `(h,k)=(0,0)` is allowed ((00l) reflections on the specular rod); Part TWO uses a 4×3 least squares
     (specular constraint + three starting peaks) to solve `(κ,λ,ρ)`, and performs a q_spec self-consistency check on the refined solution.
     It has the strongest noise resistance and fallback capability, at about 3–4× the runtime of PeriodVote/DDL.
   - `RANSAC (this project)`: Part ONE is unchanged, while Part TWO uses a **2-peak minimal sample + specular row** (3 equations) to solve
     `(κ,λ,ρ)` directly, with `[(2L+1)R]²` hypotheses (52× fewer than the 3-peak enumeration), and early exit on per-peak `g_z` residuals; deterministic, no random numbers.
     On DIP: 0.73 s / 32 solutions (the first solution matches the manual solution), about 10 ms per candidate.
   - `Integer least squares (this project)`: difference vectors `|Δg_z|/n` are clustered by **relative tolerance** to determine `ρ` (no histogram bins, no harmonic heuristics),
     then an **interval-endpoint sweep** over `t` picks the solution covering the most peaks (exact global optimum for a given `(ρ,w)`, no top-K truncation).
     On DIP: 0.40 s, the first solution matches the manual solution.
   - `SCRM (w-eliminated two-step method)`: the first step enumerates `(u,v)` only (`w` does not enter the search at all), using `t = q_z/q_s` to solve the in-plane metric
     `(A,B,C) = (z_a², z_b², z_a z_b cos γ)`; the second step solves a 3×3 system from three seed peaks to obtain `(κ,λ,ρ)`,
     and `w` and the per-peak `l` are **all recovered analytically**, with no more `l` triple enumeration. Two pitfalls are handled: both sign branches of
     `S3 = κu+λv+ρw` must be tried (`q_s²−g_uv = S3²` only gives `|S3|`), and the sign must follow `sgn` when solving `l` analytically;
     the seed peaks need `(h,k)` pairwise non-collinear and at least one `l≠0`, otherwise the 3×3 system degenerates near a `(001)` contact plane and `ρ` cannot be solved
     (the implementation falls back to multiple seed sets).
4. **Select a mode** (toolbar dropdown / "Compute" menu): after switching among the three modes, the available parameter pages change accordingly.
5. **Start indexing** (`F5`): runs in a background thread, with the status bar showing progress and stage information;
   `Esc` requests a stop (responded to between candidates).
6. **View the results**: the "Candidate Solutions" table on the right lists `(uvw)`, unit cell, volume, q_spec and the three RMSDs;
   after selecting a row, the peak table on the left is filled with that solution's per-peak `(h,k,l)` and `Δq_xy / Δq_z`,
   and the scatter plot overlays the predicted peak positions (red crosses).
7. **Export** (`Ctrl+S`): the current solution is exported as CSV, with the same layout as the original:
   an 11-row parameter header + one row `[q_xy q_z h k l]` per data point.

Rows can be appended/deleted in the peak table in the UI (`Ctrl+N` / delete selected rows), or the table can be cleared and recomputed.

## Parameters (Default Values)

- **With q_spec**: Part 1 algorithm=PeriodVote/DDL (improved / RANSAC / integer least squares / SCRM optional),
  per-peak (h,k) candidates R=4, (u,v) limit=2,
  Part1 (h,k) range=3, indexing (h,k) range=6, l range=6, w limit=2, starting three reflections take the first 4 peaks,
  SCRM seed l radius=1 (automatically widened to 3 when there is no "clearly good" solution),
  Part1 candidate limit=100, refinement candidate count=200, maximum number of output solutions=50;
  Niggli reduction, optional unit-cell range limits (a/b/c 3–30 Å, volume 50–1e5 Å³, angles 60–120°),
  dqxy / dqz / dqxyz cutoff=0.05, dqspec=0.01.
- **Without q_spec**: method (default **beam search**, optionally LLL+sphere-decoding assignment, see
  the file header of `src/core/nospec.cpp`), crystal system (triclinic by default), index enumeration box h/k/l=2/2/6 (a long c axis needs
  a sufficiently large `l`), seed index shell `|hkl|²≤6`, first 6 seed peaks, minimum 8 rows for the metric solve, beam width 8000,
  reassignment grid 6/6/6, refinement 300, at most 50 solutions.
- **Fixed contact plane**: (u,v,w) defaults to (0,0,1), h/k/l range 4, initial values a=10 b=10 c=15,
  seed peaks 5, 2000 trials, 6 alternating refinements, keep 10 solutions; q_xy is always computed as `|q_xy|`
  (the original "auto-compensate q_xy sign" switch is still present in the UI, but no longer affects the results).

## Self-Test

```powershell
build\Release\GIDSelfTest.exe
```

Coverage: forward-model consistency, Frame round-trip, agreement of the two ν algorithms, automatic q_spec detection,
the with-q_spec two-step method, the no-q_spec single-step method (both the beam search / B2 LLL+sphere-decoding methods run on
the monoclinic and cubic cases), fixed contact plane (exact convergence from ground-truth hkl / convergence from near-truth initial values / unit-cell recovery).

Current status:

- With-q_spec two-step method (`Improved` algorithm): PASS, a synthetic triclinic cell reproduces the unit cell and hkl exactly;
- With-q_spec two-step method (`Improved` algorithm, DIP measured-data regression): PASS, with `linesToPermute=4 /
  uvLimit=2 / hkLimit1=2` it yields `uvw=(-1 2 1)`,
  `a≈7.16 b≈8.49 c≈16.70 Å`, `rmsQxyz≈0.004`, consistent with the manual;
- With-q_spec two-step method (`PeriodVote/DDL` algorithm): PASS, both the same synthetic data and the DIP data give the
  same unit cell and hkl as the improved version;
- SCRM identities: PASS, all four residuals `|g|²=pla+gp²`, `q_z·q_s=s3+gp·ζ`, `q_s²=guv+ζ²`,
  `q_xy²=⟨(H,K)⟩+(Hκ+Kλ+Lρ)²` are < 5e-15;
- SCRM two-step method: PASS, both the `uvw=(1,2,3)` synthetic triclinic cell and the DIP measured data recover the correct unit cell
  (DIP: `a=7.157 b=8.491 c=16.696`, `uvw=(-1 2 1)`, `rmsQxyz=0.0039`, 2 solutions);
  the synthetic `(001)` degenerate case (triclinic cell with an in-plane candidate rms=0) was fixed and recovered from 0 solutions to 2 solutions;
- No-q_spec single-step method: PASS; on the DIP measured data (11 peaks after removing the specular rod peak `{0.000, 1.776}`) the 1st solution
  is the manual cell `a=7.12 b=8.52 c=16.65 Å / α=90.6 β=92.4 γ=90.0`, `V=1009 Å³`,
  `rmsQxyz≈0.003`, maximum per-peak deviation 0.012, about 6 s;
- Fixed contact plane (FP-DL): PASS, zero residual for the ground-truth hkl assignment, exact convergence from near-truth initial values, exact
  recovery of the orthorhombic `10×12×17` cell (`rmsQxyz=0.00000`); the 3 legacy failures (convergence from near-truth initial values,
  unit-cell parameter recovery, RMSD) were fixed along with the progressive assignment in `refineAssign` and the FP-DL refactor.

Setting the environment variable `GID_DEBUG=1` makes Part ONE print seed and candidate statistics to stderr.

## Sources

The algorithms are consistent with the MATLAB original *GIDInd* (GIDInd_2021_04_05, E. Kerscher et al.),
the original functionality and algorithms belong to the original authors; this project is for localized use and performance improvement of a scientific tool.

### Original Implementation and Citations

- Reference implementation (MATLAB, original author's repository): <https://github.com/m-kainz/GIDInd>
- Paper: Kainz, M. P., Legenstein, L., Holzer, V., Hofer, S., Kaltenegger, M., Resel, R. &
  Simbrunner, J. (2021). *GIDInd: an automated indexing software for grazing-incidence
  X-ray diffraction data.* J. Appl. Cryst. **54**, 1256–1267.
  DOI: [10.1107/S1600576721006609](https://doi.org/10.1107/S1600576721006609)
- Simbrunner, J. et.al.;. Acta Crystallogr A Found Adv 2018, 74 (Pt 4), 373-387. 
  DOI: [10.1107/S2053273318006629](https://doi.org/10.1107/S2053273318006629)


**If you use the "with q_spec two-step method" (`Improved` / `PeriodVote/DDL` / `RANSAC` /
`Integer least squares` / `SCRM`) of this program, please cite the above references in your publications.**

`SCRM` (the w-eliminated two-step method) is a two-step method independently implemented in this project; it shares the `(A,B,C)+(κ,λ,ρ)` parameterization and the
refinement/reduction kernel with GIDInd, but the search method is different; when citing, it is likewise recommended to list the above references as well.
