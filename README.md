# Designed interfaces vs natural interfaces

Do de novo designed protein interfaces occupy the same structural space as
natural ones?

This repository adds one experiment to the Foldseek-Interface analysis code:
47,178 BindCraft-designed interfaces and 47,178 natural PDB interfaces are
clustered **together in a single run**, and we count how many clusters end up
containing both.

The notebook is [`Iface-HumanBindome.ipynb`](Iface-HumanBindome.ipynb).
Everything else in this repository is the original Foldseek-Interface analysis
code from [steineggerlab/foldseek-interface-analysis](https://github.com/steineggerlab/foldseek-interface-analysis).

---

## Result

```
interfaces   94,286   (47,178 designed + 47,108 natural)
clusters     58,614

  designs only    31,930   54.5%
  natural only    26,149   44.6%
  both               535    0.9%

designs sharing a cluster with a natural interface:  3.57%
```

**96.4% of designed interfaces land in a cluster containing no natural
interface at all.**

Restricting to clusters with at least two members and comparing against a
null where source labels are shuffled across members:

| cluster size | observed mixed | expected if random |
|---|---|---|
| 2 | 1.6% | 50.0% |
| 3 | 4.3% | 75.0% |
| 6 | 9.6% | 96.9% |
| 10 | 22.2% | 99.8% |
| **all ≥2** | **4.3%** | **66.7%** |

The pool is 50.04% designed, so chance alone would put both sources in
two-thirds of multi-member clusters. Observed is 15× below that. Designed and
natural interfaces do not merely fail to mix — they segregate.

---

## Data

| source | what it provides |
|---|---|
| [HumanBindome](https://huggingface.co/datasets/wjulius/HumanBindome) | 306,137 BindCraft target–binder complexes, `metadata.csv` + 8,296 tar archives |
| [Zenodo 10.5281/zenodo.22040892](https://doi.org/10.5281/zenodo.22040892) | `pdb_clusters_annotated.tsv`, 3,121,961 natural chain pairs with cluster labels |
| RCSB | the assemblies themselves, fetched on demand |

The Zenodo release ships cluster labels but no coordinates, so every natural
interface structure in this experiment is downloaded from RCSB and cut here.

**foldseek** must be the interface build (commit `463739e0` or later). Released
versions, including bioconda's, do not have `createinterfacedb`.

---

## Procedure

### 1. Extract the designed complexes

Each archive in `structures/` holds 102 members for one target: the apo target,
the apo domain, and 50 complex / 50 apo-binder pairs. Only the complexes are
needed. `metadata.csv` names them in `target_binder_complex_pdb_path`.

```
306,137 two-chain complexes
```

### 2. Build chain and dimer databases

```bash
foldseek createdb      complexes/ chains
foldseek createdimerdb chains dimers tmp \
         --distance-threshold 10 --min-interface-residues-perchain 4
```

```
612,274 chains → 306,137 dimers
```

### 3. Deduplicate the designs

Sister ProteinMPNN variants from the same seed share a backbone and are
redundant. TM 0.9 is strict enough that genuinely distinct interfaces never
merge. No exhaustive search at this stage — that is what keeps 306,137 dimers
tractable.

```bash
foldseek multimercluster dimers dedup tmp \
         --multimer-tm-threshold 0.9 --chain-tm-threshold 0.5 \
         --interface-lddt-threshold 0 -c 0 --cov-mode 0 --cluster-mode 0
foldseek createtsv dimers dimers dedup dedup_cluster.tsv
```

```
306,137 designs → 188,710 representatives   (62% retention, ~2.7 h)
```

### 4. Build representative interfaces

Rebuild from the representatives' source PDB files, then sever.

```bash
foldseek createdb          reps/ rep
foldseek createdimerdb     rep rep_dimer tmp  <interface flags>
foldseek createinterfacedb rep_dimer rep_int  <interface flags>
```

```
188,710 interfaces, covering 8,296 targets
```

> **Note on `createsubdb`.** The published `scripts/cluster_pdb.sh` selects
> representatives with `createsubdb clu PDB_dimer PDB_dimer_reps`. On foldseek
> `463739e0` that selects at cluster-key granularity — 188,710 entries where a
> dimer DB needs 377,420 — and carries the parent `.lookup` over unchanged.
> `createinterfacedb` then yields **35,712 interfaces instead of 188,710**, with
> exit code 0 and no warning. Rebuilding from the source structures avoids this.
> Worth checking the count explicitly whichever route is taken.

### 5. Sample both sides

**Designed:** 47,178 of the 188,710 representatives, uniform.

**Natural:** uniform over the 189,830 rows with `dimerrep == 1` — the authors'
own TM 0.9 deduplicated set, and the same point in the pipeline the designed
side has reached. Sampling raw rows would pit deduplicated designs against
undeduplicated PDB.

No centroid enrichment. A uniform draw contains 40.7% interface-cluster
centroids, which is the PDB's composition at this level; forcing it higher
would inflate mixing for reasons unrelated to the designs.

```
47,178 per side   (25.0% and 24.9% of the respective pools)
spanning 25,884 distinct natural interface clusters
```

### 6. Fetch and cut the natural interfaces

Each sampled row names an assembly and two chains. Assemblies are downloaded,
then the two named chains are written to their own file.

```
23,274 assemblies needed → 23,262 fetched (12 withdrawn from RCSB)
47,129 two-chain structures
```

mmCIF column order is not fixed, so `auth_asym_id` is located per file rather
than by index.

### 7. Stage one pool

```
pool/
  W<row>C<intcluster>.cif    natural
  <design_name>.pdb          designed
```

The `W` prefix is the only thing distinguishing the sources downstream; design
names always contain `no_hotspot`, so the two cannot collide.

### 8. Build interfaces and cluster

```bash
foldseek createdb          pool/ pooldb
foldseek createdimerdb     pooldb pooldimer tmp  <interface flags>
foldseek createinterfacedb pooldimer int        <interface flags>

foldseek easy-multimercluster int clu tmpc \
         --exhaustive-search 1 -e 10000000 --lddt-threshold 0.20 \
         --chain-tm-threshold 0 -c 0 --interface-lddt-threshold 0 \
         --cov-mode 0 --cluster-mode 0 --max-iterations 1 \
         --multimer-tm-threshold 0.4
```

`createinterfacedb` keeps only contact residues, concatenates the two chains
and recomputes 3Di — roughly 35 residues from a 136-residue complex. What gets
clustered is the interface, not the parent fold.

```
94,286 interfaces → 58,614 clusters
```

### 9. Count

Per cluster: designed only, natural only, or both.

---

## Cluster size distribution

```
size      designed   natural   both     total
1           25,060    21,024      0    46,084
2            3,702     2,683    101     6,486
3-4          2,392     1,525    193     4,110
5-9            710       641    162     1,513
10-49           65       247     76       388
50+              1        29      3        33
```

78.6% singletons. Largest cluster 2,056 members, all natural.

Two features worth noting.

**Mixed clusters are concentrated at larger sizes** — none at size 1, 1.6% at
size 2, 20% at size 10–49. A cluster must be reasonably large before it can
contain both, so the 3.57% is partly a statement about cluster size.

**The two arms have different size profiles.** The designed side has one
cluster above 50 members; the natural side has 29, including a 2,056-member
one. The largest natural clusters are near-duplicate families — amyloid
fibrils (`IC33341`, 10,996 members in the catalogue), gas vesicle protein,
viral capsids — where a single deposited structure contains hundreds of copies
of the same subunit contact. Designs form no equivalent.

---

## Interpreting the number

3.57% is uninterpretable on its own, because the maximum achievable rate is not
100%. The control is to run the identical experiment with natural interfaces on
**both** sides, split so the two halves share no protein. That split matters:
without it, most matches are two depositions of the same complex — a match a de
novo design can never produce.

At the same thresholds, on a smaller sample:

| interface TM | designs vs natural | natural vs natural | ratio |
|---|---|---|---|
| 0.4 | 0.51% | 5.62% | 11.0× |
| 0.3 | 3.44% | 11.02% | 3.2× |
| 0.25 | 7.38% | 17.23% | 2.3× |
| 0.2 | 13.74% | 28.81% | 2.1× |

Loosening the threshold raises both columns, so a bare cross-source number at a
looser cutoff means nothing. The ratio narrows from 11× to 2.1×, but designs
never reach parity.

Both numbers are low for the same structural reason: in a 47k sample, 81.3% of
the natural interface clusters present contribute exactly one structure, so
most structures have no possible partner in the pool. This is a property of
sampling density in interface space, not of either source.

The authors' own benchmark found TM 0.3 optimal; 0.4 was used for the full PDB
clustering. 0.25 and 0.2 are below anything they tested.

---

## Running it

Edit five lines in the notebook's **Setup** cell:

```python
FOLDSEEK  = ".../foldseek"                      # interface build
BINDOME   = Path(".../HumanBindome")            # the HF download
COMPLEXES = Path(".../complexes")               # output of step 1
PDB_TABLE = Path(".../pdb_clusters_annotated.tsv")
WORK      = Path(".../repro")                   # everything is written here
```

Everything else derives from those. `--threads 48` is set in the foldseek cells.

Approximate cost on 48 threads:

| step | time | peak |
|---|---|---|
| extract complexes | ~1 h | 98 GB disk |
| chain + dimer DBs | 90 s | |
| dedup TM 0.9 | 2.7 h | 107 GB RAM |
| representative interfaces | 35 s | |
| sample | 1 min | |
| fetch assemblies | 17 min | 15 GB disk |
| stage pool | 35 min | 30 GB disk |
| interface DB | 31 s | |
| cluster TM 0.4 | several hours | ~580 GB scratch |

The clustering scratch directory is removed only on success; delete it by hand
if a run is interrupted.
