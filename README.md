# Polyjam

Polyjam is a small research playground for **packing identical convex polyhedra inside a convex polyhedral shell**. The piece and shell may be different shapes. Every piece gets an independent translation and rotation; the shell is centered and uniformly scaled. The optimizer minimizes that shell scale.

It is intentionally two things at once:

- a dependency-light C++20 search backend built around exact convex geometry tests, and
- a local HTML inspector that can launch searches and orbit/slice/explode saved results.

## What is exact, and what is heuristic?

Geometry checks are deterministic for a given pose (up to floating-point tolerance):

- **Containment:** every transformed piece vertex is tested against every support plane of the convex shell.
- **Piece overlap:** the 3D separating-axis theorem tests piece face normals plus all edge × edge axes.
- **Convex hull:** OBJ input only needs `v x y z` records. Polyjam infers support planes/faces, recenters the hull, and normalizes its circumradius to 1.

The **search is heuristic**. It uses parallel stochastic local searches with translation, quaternion rotation, shell-pressure moves, annealing, constraint-penalty tightening, and occasional shakes. Searches are seeded with feasible lattice constructions fitted to both the piece's oriented bounding box and the shell's support planes, rather than an overlapping random cloud. A run finding scale `2.71` is an upper bound/construction, not a proof that `2.71` is optimal.

That distinction is reflected in every result file: `scale`, `feasible`, `maxViolation`, and the search budget are all stored separately.

## Build and run the browser lab

Requirements: a C++20 compiler, CMake 3.20+, and Python 3.10+ (the Python server uses only the standard library).

```bash
./run.sh
```

Then open <http://127.0.0.1:8080/>.

`run.sh` builds `build/polyjam` and starts a localhost-only HTTP server. The browser can launch builtin-shape searches and polls their progress. Results are written to `web/results/` and can also be opened directly from JSON files.

## CLI

```bash
./build/polyjam search \
  --piece cube \
  --shell cube \
  --count 10 \
  --seconds 30 \
  --threads 8 \
  --seed 123 \
  --output web/results/cube10.json
```

Builtin convex shapes:

```text
tetra  cube  octa  icosa  dodeca
```

The piece and shell need not match:

```bash
./build/polyjam search --piece tetra --shell dodeca --count 12 --seconds 60 --output web/results/tetra-in-dodeca.json
```

For custom convex polyhedra, pass an OBJ file as either shape. Faces in the OBJ are ignored; the convex hull of all `v` records is used:

```bash
./build/polyjam search --piece examples/wedge.obj --shell icosa --count 7 --seconds 60 --output web/results/wedges.json
```

All loaded shapes are recentered at their vertex centroid and scaled to circumradius 1. Consequently, when piece and shell are the same shape, `scale` is simply the outer/inner linear size ratio.

## Search knobs

- `--seconds`: wall-clock budget shared by all workers.
- `--threads`: independent stochastic workers sharing a global best result.
- `--seed`: reproducible master seed.
- `--start-scale`: override the loose automatic starting scale.
- `--tolerance`: maximum permitted geometric penetration/excess for a result to be marked feasible.
- `--clearance`: request positive separation between pieces (useful for fabrication/visual breathing room).
- `--progress-json`: newline-delimited machine-readable progress events; used by `tools/server.py`.

The search currently targets small/medium instances. SAT cost grows roughly quadratically with piece count and with the number of edges in the piece.

## Frontend controls

The browser renderer is plain Canvas 2D with a tiny software 3D projection—no Three.js/WebGL dependency. Drag to orbit, wheel to zoom, click a piece for its pose, adjust shell opacity, use an exploded view, or sweep a Z slice through the packing.

## Result format

`polyjam-result-v1` embeds the normalized piece and shell geometry plus all poses, so a result is portable and inspectable without access to the original OBJ:

```json
{
  "format": "polyjam-result-v1",
  "piece": { "vertices": [], "faces": [] },
  "shell": { "vertices": [], "faces": [] },
  "count": 10,
  "scale": 2.72,
  "feasible": true,
  "metrics": { "maxViolation": 0.0 },
  "poses": [
    { "p": [0, 0, 0], "q": [1, 0, 0, 0] }
  ]
}
```

Quaternions are `[w,x,y,z]`.

## Where I would take it next

The current stochastic engine is deliberately compact and hackable. Natural upgrades are:

1. contact-graph-aware local moves and per-piece violation gradients;
2. a deterministic polishing stage after a good contact topology is found;
3. branch-and-bound / interval arithmetic for **proofs** on small instances;
4. symmetry reduction for identical pieces and highly symmetric shells;
5. checkpoint/resume plus a persistent run database;
6. SIMD/GPU batching of SAT tests for larger populations.

The clean boundary is `evaluate(...)`: replacement optimizers can reuse the same hull, containment, SAT, result, server, and viewer machinery.

## License

MIT. See `LICENSE`.
