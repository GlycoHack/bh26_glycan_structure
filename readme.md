# WURCS → SMILES Batch Conversion and ResCode Extraction

Bulk conversion of GlyTouCan WURCS strings to SMILES using `molwurcs.jar`, followed by
extraction of the residue codes (ResCode), SkeletonCodes and MAP codes contained in the
WURCS strings.

All paths in this document are relative to the root of this repository.

`molwurcs.jar` is **not tracked** here: it is a 45 MB build artefact of another project and,
being a zip already, git can neither compress nor delta it. Place the MolWURCS 0.12.2 jar in
this directory before running anything — build it from
[glycoinfo/MolWURCS](https://gitlab.com/glycoinfo/molwurcs) with `mvn package`, which writes it
to `target/molwurcs.jar`. `WURCSFilter.jar` (section 9) is untracked for the same reason. Every
other file below is committed as-is.

- Date of run: 2026-09-18
- Environment: macOS 26.5.2, 12 CPU cores, 64 GB RAM, OpenJDK 17.0.7 (Temurin)
- Converter: `molwurcs.jar` (MolWURCS 0.12.2), WURCS parsing per `wurcsframework` 1.3.1

---

## 1. Input data

| File | Description | Lines | Size |
|---|---|---|---|
| `input.txt` | Source data. Tab-separated, two columns, every field double-quoted: `"WURCS=..."<TAB>"G12345XX"` | 256,453 | 61 MB |
| `gtc.txt` | Byte-identical to `input.txt` (verified with `cmp`); the original download | 256,453 | 61 MB |
| `output.txt` | Intermediate from an earlier attempt: column 1 with quotes stripped, i.e. WURCS only. Byte-identical to `cut -f1 input.txt \| tr -d '"'` | 256,453 | 58 MB |

Properties verified on `input.txt`: every line has exactly 2 tab-separated fields; there are
no blank lines, no `#` comment lines, and no empty ID fields. Row alignment with the input is
therefore exact throughout the pipeline.

---

## 2. The conversion script

### Usage

```bash
./convert_wurcs.sh [input_file] [output_file]
```

Defaults are `input.txt` → `wurcs-smi-id.txt`.

| Environment variable | Default | Meaning |
|---|---|---|
| `JOBS` | `4` | Number of JVMs run in parallel |
| `CHUNK_LINES` | `8000` | Lines processed per JVM invocation |
| `ONLY_SUCCESS` | `0` | Set to `1` to emit only rows that converted successfully |

```bash
JOBS=6 ./convert_wurcs.sh            # more parallelism
ONLY_SUCCESS=1 ./convert_wurcs.sh    # drop failed rows
```

### Output format

Tab-separated, three columns, in the same row order as the input:

```
ID <TAB> WURCS <TAB> SMILES
```

Rows whose conversion failed keep their ID and WURCS and have an **empty SMILES column**,
so the output always has the same number of rows as the (filtered) input.

### How it works

1. **Extract and clean.** `awk` takes column 1 (WURCS) and column 2 (ID) as a pair, strips
   double quotes and surrounding whitespace, and drops blank and `#` lines. Because the ID is
   carried alongside its WURCS through this step, filtering cannot desynchronise the columns.
2. **Convert in parallel.** The WURCS column is split into chunks of `CHUNK_LINES` lines
   (`split -a 4`), and `xargs -P "$JOBS"` runs one JVM per chunk, each reading its chunk from
   stdin. The jar is invoked once per chunk, not once per line.
3. **Reassemble.** Chunk outputs are concatenated in filename order, which is the original
   row order. The script aborts without writing the output file if the line count does not
   match the input.
4. **Join.** `paste` combines ID, WURCS and SMILES into the three-column result.

Two details make the row order safe:

- `-n` (`--output-no-result`) makes the jar emit an empty line for every failed conversion,
  so each chunk's output line count equals its input line count.
- `split -a 4` produces lexicographically ordered chunk names, so plain concatenation
  restores the original order.

### Correct jar arguments

The format names accepted by `molwurcs.jar` are lowercase short names:

```bash
java -jar molwurcs.jar -i wurcs -o smi -n < chunk.txt
```

`-o SMILES` and `-i WURCS` are **not** valid and fail with
`Output format must be specified.` Valid values for both `-i` and `-o` are
`wurcs`, `mdlv2000`, `sdf`, `sdfv3000`, `smi`.

### Previous version

`convert_wurcs.sh.bak` is the original script, kept for reference. It produced a 0-byte
output file for three reasons:

1. It passed `-o SMILES` / `-i WURCS`, which are invalid, so every single conversion failed.
   The output file had already been truncated by `> "$OUTPUT_FILE"` at the top of the script,
   so it stayed empty.
2. It fed the raw TSV line (quotes and ID column included) to the jar instead of the bare
   WURCS string.
3. Its `echo "$wurcs$smi"` had no tab between the two values, so successful rows would not
   have been split into columns.

It also launched one JVM per line — 256,453 JVM startups.

---

## 3. Results

Full run over all 256,453 records:

```
real  6m1.845s    user  26m49.178s
converted: 85,114 / 256,453  (33.2%)
failed:    171,339
```

Failure breakdown, from the error log of the full run:

| Reason | Count |
|---|---|
| `This WURCS cannot be converted to molecule` / `Incompatible` | 165,771 |
| Exception | 24 |
| `Unknown CarbonDescriptor is contained.` | 1 |

These failures come from MolWURCS itself, not from the script. The overwhelming majority are
WURCS strings that cannot be expanded to an atomic-level structure because they describe the
glycan incompletely — for example `a2122h-1x_1-5` (unknown anomeric configuration) or
`axxxxh-1x_1-?` (unknown carbon descriptors and unknown linkage position). GlyTouCan contains
many such partially-determined entries, so a failure rate of this size is expected.

The rate is not uniform across the file: the first 20,000 records converted at 55%, while the
file as a whole converted at 33%, i.e. indeterminate structures are more frequent in the
later part of the data.

---

## 4. Output files

### Conversion results

| File | Columns | Rows | Size |
|---|---|---|---|
| `wurcs-smi-id.txt` | ID, WURCS, SMILES — all records, empty SMILES where conversion failed | 256,453 | 78 MB |
| `wurcs-smi-id-ok.txt` | ID, WURCS, SMILES — successful conversions only | 85,114 | 28 MB |
| `wurcs-smi.txt` | WURCS, SMILES — earlier two-column form, all records | 256,453 | 75 MB |
| `wurcs-smi-ok.txt` | WURCS, SMILES — earlier two-column form, successes only | 85,114 | 27 MB |

The two-column files predate the three-column ones and are kept only for reference; the
three-column files carry the same SMILES data plus the GlyTouCan ID.

### ResCode extraction

Derived from `wurcs-smi-id-ok.txt`, i.e. only from glycans that converted to SMILES.

| File | Columns | Rows | Size |
|---|---|---|---|
| `rescode-by-id.txt` | ID, ResCode — long format, one row per residue code per glycan | 225,757 | 7.5 MB |
| `rescode-unique.txt` | ResCode, number of glycans containing it — descending | 35,711 | 1.7 MB |
| `skeletoncode-unique.txt` | SkeletonCode, number of glycans, number of distinct ResCodes it appears in | 2,585 | 32 KB |
| `mapcode-unique.txt` | MAP code, number of glycans, number of distinct ResCodes it appears in | 2,615 | 111 KB |
| `mapcode-features.txt` | MAP code, counts, and parser-derived structural features (see section 8) | 2,615 | 198 KB |
| `mapcode-class.txt` | MAP code, feature flags, `primary_class` (see section 8) | 2,615 | 177 KB |

ResCodes are written without the enclosing square brackets. MAP codes **retain their leading
`*`**, matching the framework's own representation (`t_strMAP = substring(indexOf("*"))`).

### Logs

| File | Description |
|---|---|
| `convert_wurcs.log` | stdout of the full run (progress, final counts) |
| `convert_wurcs.err` | stderr of the jar, one entry per failed WURCS with the reason |

**Note:** `convert_wurcs.err` is truncated at the start of every run. The file currently
present is from a later 3,000-line test run, **not** from the full run whose numbers are
quoted in section 3. Re-run the full conversion to regenerate the complete error log.

---

## 5. Performance

One JVM already uses about 3.3 CPU cores (measured 69.8 s user / 21.1 s real), so
parallelism saturates the 12-core machine at around `JOBS=4`. Going higher oversubscribes it
and is slower.

| `JOBS` | Throughput |
|---|---|
| 1 (serial) | 377 rows/s |
| 4 | 963 rows/s |
| 6 | 937 rows/s |
| 12 | 827 rows/s |

Measured on 8,000–24,000 row samples. The chosen default (`JOBS=4`, `CHUNK_LINES=8000`) gives
about 2.5× the serial throughput, converting 26.8 minutes of CPU time in 6 minutes of wall
clock.

---

## 6. Caveats

### SMILES output is not canonical

`molwurcs.jar` emits non-canonical SMILES: **the same input produces different atom ordering
on different runs.** Two consecutive serial runs over the same 2,000-row chunk differed in 848
and 613 rows respectively. This is jar behaviour and is unrelated to parallelisation.

The molecules are nevertheless identical — converting two differing SMILES for the same row
back with `-i smi -o wurcs` yields the same original WURCS string. Do not use these files for
byte-level diffing between runs, and do not treat the SMILES as a structural key.

### Non-pipeline files

`a.sh`, `convert_wurcs_fixed.sh` (0 bytes) and `input_fixed.txt` (0 bytes) are earlier
scratch files and are not used by this pipeline. Note that `a.sh` is Shift_JIS-encoded and
displays as mojibake in a UTF-8 terminal; `convert_wurcs.sh` and this readme are UTF-8.

---

## 7. ResCode / SkeletonCode / MAP code extraction

These steps were run ad hoc after the conversion; the commands are recorded here for
reproducibility.

### ResCode

A WURCS string lists its unique residue codes in square brackets in the third section:

```
WURCS=2.0/6,9,8/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5].../1-1-2-3-1-4-5-3-6/a4-b1_a6-i1_...
```

Because a ResCode itself contains `/` and `*`, the string cannot be split on `/`. Brackets do
not nest (verified: zero occurrences of `[` before a closing `]`), so a non-greedy bracket
match is safe:

```bash
# ID <TAB> ResCode, long format
awk -F'\t' 'BEGIN{OFS="\t"} {
  s = $2
  while (match(s, /\[[^]]*\]/)) {
    print $1, substr(s, RSTART+1, RLENGTH-2)
    s = substr(s, RSTART+RLENGTH)
  }
}' wurcs-smi-id-ok.txt > rescode-by-id.txt

# ResCode <TAB> glycan count, descending
cut -f2 rescode-by-id.txt | sort | uniq -c | sort -rn \
  | awk '{c=$1; $1=""; sub(/^ /,""); print $0 "\t" c}' > rescode-unique.txt
```

### SkeletonCode and MAP code

The decomposition follows `WURCSImporter.extractMS()` and `extractMOD()` in
`wurcsframework` 1.3.1:

1. Split the ResCode on `_`.
2. **Field 1** is `SkeletonCode[-AnomericPositionSymbol]`; the SkeletonCode is everything
   before the first `-`.
3. **Fields 2..n** are MODs. Only those containing `*` carry a MAP; the MAP code is the
   substring from the **first** `*` onward (inclusive). Fields without `*`, such as the ring
   closure `1-5`, have no MAP code.

Splitting on every `*` would be wrong. A `*` is a `MAPStar`: an attachment point to a backbone
carbon. A modification that bridges two or more backbone carbons carries one numbered star per
attachment, so a single MAP code can hold several — of the 2,615 distinct MAP codes, 1,997 have
one star, 599 have two, 18 have three and 1 has four. Combined across the MOD fields of one
ResCode, up to 12 `*` occur. Only the **first** `*` in each `_`-field is the MOD/MAP delimiter;
the framework uses `indexOf("*")`.

Two further MAP characters matter when reading these codes, and neither means what a
SMILES-trained eye expects:

- `$n` is a `MAPAtomCyclic`: a **ring closure** back to atom *n*.
- `(` and `)` do **not** delimit a branch — they switch **aromaticity** on and off
  (`MAPGraphImporter` sets the aromatic flag for atoms inside them). Branches are written with
  the `/<atom><bond><atoms>` connection syntax instead.

So `*OC(CCCCCC$4)` is not an alkyl chain: it is an O-CH2 linker carrying six aromatic carbons
closed into a ring — a benzyl ether. Adding `/3=O` puts a carbonyl on atom 3, giving a benzoyl
ester.

```bash
awk -F'\t' 'BEGIN{OFS="\t"}
{
  res=$1; g=$2
  nf = split(res, f, "_")
  p = index(f[1], "-")
  sc = (p ? substr(f[1], 1, p-1) : f[1])
  skel[sc] += g; skel_n[sc]++
  for (i = 2; i <= nf; i++) {
    q = index(f[i], "*")
    if (q == 0) continue
    map = substr(f[i], q)
    mp[map] += g; mp_n[map]++
  }
}
END{
  for (k in skel) print k, skel[k], skel_n[k] > "/tmp/sc.raw"
  for (k in mp)   print k, mp[k],   mp_n[k]   > "/tmp/mp.raw"
}' rescode-unique.txt
sort -t$'\t' -k2,2nr -k1,1 /tmp/sc.raw > skeletoncode-unique.txt
sort -t$'\t' -k2,2nr -k1,1 /tmp/mp.raw > mapcode-unique.txt
```

### Distribution

ResCodes per glycan: mean 2.65, maximum 16 distinct residue codes.

Of the 35,711 distinct ResCodes, 29,941 (84%) occur in only one glycan, 4,898 in 2–10, and
872 in 11 or more. The long tail consists mostly of heavily modified residues carrying lipid
chains or cyclic substituents. 2,751 ResCodes carry no MAP code at all.

Most frequent SkeletonCodes (by number of glycans):

| SkeletonCode | Glycans | Typical residue |
|---|---|---|
| `a2122h` | 55,186 | Glc / GlcNAc series |
| `a2112h` | 39,473 | Gal / GalNAc series |
| `a1122h` | 33,434 | Man series |
| `a1221m` | 12,440 | Fuc |
| `Aad21122h` | 9,361 | Sialic acid (Neu5Ac) |
| `a2211m` | 4,253 | Rha |
| `a2122A` | 3,146 | GlcA |
| `axxxxh` | 2,898 | Unknown configuration |

SkeletonCode length equals the carbon count; hexoses (6 characters) are the largest group
with 1,067 distinct codes.

Most frequent MAP codes (by number of glycans):

| MAP code | Glycans | Modification |
|---|---|---|
| `*NCC/3=O` | 56,026 | N-acetyl |
| `*OCC/3=O` | 23,318 | O-acetyl |
| `*OC(CCCCCC$4)` | 18,535 | Benzyl ether (aromatic ring) |
| `*OC(CCCCCC$4)/3=O` | 10,627 | Benzoyl ester (aromatic ring) |
| `*OC` | 10,093 | O-methyl |
| `*OSO/3=O/3=O` | 8,159 | Sulfate |
| `*N` | 7,356 | Amino |
| `*OPO/3O/3=O` | 2,679 | Phosphate |
| `*F` | 1,673 | Fluoro |

---

## 8. MAP code feature classification

MAP codes are classified by structural features — lipid chains, aromatic rings, macrocycles and
so on. Because a MAP code's grammar is deceptive (see the notes in section 7), the features are
**not** derived from string patterns. `MapFeatures.java` parses each code with the framework's own
`MAPFactory`, resolves `MAPAtomCyclic` references into ring-closure edges, and measures the
resulting graph.

```bash
javac -cp molwurcs.jar MapFeatures.java
cut -f1 mapcode-unique.txt | java -cp molwurcs.jar:. MapFeatures > features.raw
```

All 2,615 MAP codes parsed without error.

### Files

| File | Description | Rows |
|---|---|---|
| `mapcode-features.txt` | Objective per-code measurements: atom counts by element, ring count and sizes, aromatic atom count, carbonyl and ester counts, longest acyclic carbon chain, attachment atom | 2,615 |
| `mapcode-class.txt` | Feature flags (multi-label) plus one `primary_class` per code | 2,615 |
| `MapFeatures.java` | The feature extractor | — |

`mapcode-class.txt` carries independent 0/1 flags for `aromatic`, `lipid`, `multi_acyl`,
`alicyclic`, `macrocycle`, `sulfate`, `phosphate`, `halogen` and `acyl`, because these overlap —
a tosyl group is both aromatic and a sulfonate. `primary_class` assigns each code to exactly one
bucket by the priority order below, for tallying.

### Rules

| Class | Rule |
|---|---|
| `macrocycle` | largest ring ≥ 10 atoms |
| `lipid` | longest acyclic carbon chain ≥ 8 |
| `aromatic` | ≥ 1 ring whose atoms all carry the parser's aromatic flag |
| `alicyclic` | has a ring, none aromatic |
| `phosphate` | contains P |
| `sulfate` | contains S |
| `multi_acyl` | ≥ 2 C=O on a chain of ≥ 4 carbons |
| `acyl_short` | ≥ 1 C=O, no longer feature |
| `halogen` | contains F, Cl, Br or I |
| `simple` | none of the above |

`primary_class` applies these top to bottom, first match wins.

### Distribution

| primary_class | Codes | Glycan occurrences |
|---|---|---|
| `acyl_short` | 159 | 85,718 |
| `aromatic` | 1,419 | 43,268 |
| `simple` | 143 | 24,434 |
| `sulfate` | 61 | 11,197 |
| `lipid` | 333 | 6,057 |
| `phosphate` | 108 | 4,827 |
| `halogen` | 17 | 3,360 |
| `alicyclic` | 217 | 1,822 |
| `multi_acyl` | 46 | 241 |
| `macrocycle` | 112 | 113 |

Totals: 2,615 codes, 181,037 occurrences — matching the sum in `mapcode-unique.txt`.

Aromatic groups dominate the distinct codes (1,419 of 2,615) while N-/O-acetyl dominates the
occurrences. Both are expected: GlyTouCan holds many synthetic oligosaccharides carrying benzyl
and benzoyl protecting groups, each variant a distinct code, whereas acetamido appears in nearly
every natural glycan.

### Confidence per category

**Lipid — reliable.** The chain-length distribution is bimodal, and above C12 it is strongly
biased to even lengths (206 even vs 46 odd), which is the signature of biological fatty acyl
chains. The most frequent members read correctly as such: `*OCCCCCCCCCCCCCCCC/3=O` is palmitoyl
(C16), `*OCCCCCCCCCCCCCCCCCC/3=O` stearoyl (C18), and `*OCCCCCCCCC=^ZCCCCCCCCC/3=O` oleoyl, whose
`^Z` marks the cis double bond in the right position.

**Aromatic — reliable.** Aromaticity is not inferred; the parser flags it directly from the MAP
notation. The top members are benzyl (18,535 glycans), benzoyl (10,627) and tosyl (829).

**Polyketide — not determinable.** Polyketide is a classification by biosynthetic origin, and a
substituent's connectivity does not carry that information. A first attempt using "≥ 2 carbonyls
on a chain" was discarded: inspection showed it collects dicarboxylic acyl groups such as
`*OCCCCC/6=O/3=O` (glutaryl-like, 158 glycans) and diacyl glycerolipid chains, none of them
polyketides. That rule is retained under the accurate name `multi_acyl`.

The closest defensible proxy is `macrocycle`, and genuine polyketides do appear in it — for
example a 16-membered macrolactone with methyl branches, conjugated E-configured dienes and
hydroxyls,
`*C^SC^SC^ROCC=^ECC=^ECC^RC^SC^SCC=^ECC=^ECC^SOC/19$4/15C/13C/12O/11C/9C/7C/6=O/3C/2O`,
which is a classic macrolide skeleton. But the class is small (112 codes, 113 glycans) and mixed:
it also contains cyclic peptides (`*OCCOCNCCCNCOCC/12=O/6=O/2$15`) and polyamine macrocycles
(`*NCCNCCNCCNCC$2/...`). Treat `macrocycle` as a shortlist of candidates to inspect by hand, not
as a polyketide assignment.

### Adjusting the thresholds

The C8 lipid cutoff and the 10-atom macrocycle cutoff are conventions, not results. The raw
measurements in `mapcode-features.txt` let you re-cut without re-running the parser — for example
a C12 lipid threshold:

```bash
awk -F'\t' 'NR>1 && $17>=12' mapcode-features.txt
```

The full chain-length distribution is:
C0 312, C1 837, C2 317, C3 435, C4 228, C5 77, C6 44, C7 32, C8 37, C9 17, C10 21, C11 6, C12 17,
C13 9, C14 29, C15 11, C16 41, C17 7, C18 47, C19 5, C20 24, C22 23, C24 9, C26 5, and a thin
tail out to C79.

Ring sizes present: 5 (563), 6 (2,672), 7 (25), 8 (126), 10 (4), 11 (3), 12 (103), 14 (1), 16 (1),
plus one 2-atom artefact.

---

## 9. WURCSFilter (STORM) subset

The same analysis, restricted to WURCS that pass
[WURCSFilter](https://gitlab.com/glycoinfo/wurcsfilter) 0.6.1 with the **STORM** predefined
pattern set.

```bash
git clone https://gitlab.com/glycoinfo/wurcsfilter.git
cd wurcsfilter && mvn clean compile assembly:single   # -> target/WURCSFilter.jar
```

WURCSFilter takes `Title<TAB>WURCS` on stdin and writes the passing lines in the same format.
`wurcs-smi-id.txt` already holds `ID<TAB>WURCS` in its first two columns:

```bash
cut -f1,2 wurcs-smi-id.txt \
  | java -jar WURCSFilter.jar -t storm > wurcsfilter-storm-pass.txt   # 65 s
```

STORM keeps structures built from a microbial-oriented dictionary (hexoses, pentoses, heptoses,
ulosonic acids, branched and small monosaccharides) and enables fuzzy matching. The default
`snfg` set is much stricter; both were run, for comparison.

| Predefined set | Passed | of 256,453 |
|---|---|---|
| `storm` | 169,603 | 66.1% |
| `snfg` (default) | 70,136 | 27.3% |

### How the subset was built

The passing IDs were joined back onto the existing conversion results rather than re-running
MolWURCS. This is deliberate: filtering selects a subset of the same WURCS strings, and reusing
the conversion keeps the two analyses **directly comparable**, whereas a re-run would emit
different SMILES strings for the same molecules (see section 6). GlyTouCan IDs are unique across
all 256,453 records, so the join is exact; it matched all 169,603 IDs with zero WURCS mismatches.

```bash
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {pass[$1]=$2; next} ($1 in pass) {print $1,$2,$3}' \
  wurcsfilter-storm-pass.txt wurcs-smi-id.txt > wurcs-smi-id-storm.txt
awk -F'\t' '$3 != ""' wurcs-smi-id-storm.txt > wurcs-smi-id-storm-ok.txt
./analyze_rescode.sh wurcs-smi-id-storm-ok.txt storm-
```

`analyze_rescode.sh` packages the extraction and classification steps of sections 7 and 8 so
they can be re-run on any `ID/WURCS/SMILES` file:

```bash
./analyze_rescode.sh <input file> [output prefix]
```

### Files

| File | Description | Rows |
|---|---|---|
| `wurcsfilter-storm-pass.txt` | ID, WURCS — passed the STORM filter | 169,603 |
| `wurcsfilter-snfg-pass.txt` | ID, WURCS — passed the default SNFG filter | 70,136 |
| `wurcs-smi-id-storm.txt` | ID, WURCS, SMILES — STORM subset, empty SMILES where conversion failed | 169,603 |
| `wurcs-smi-id-storm-ok.txt` | ID, WURCS, SMILES — STORM subset, converted only | 54,796 |
| `storm-rescode-by-id.txt` | ID, ResCode | 178,889 |
| `storm-rescode-unique.txt` | ResCode, glycan count | 11,683 |
| `storm-skeletoncode-unique.txt` | SkeletonCode, glycan count, distinct ResCodes | 987 |
| `storm-mapcode-unique.txt` | MAP code, glycan count, distinct ResCodes | 39 |
| `storm-mapcode-features.txt` | Structural features per MAP code | 39 |
| `storm-mapcode-class.txt` | Feature flags and `primary_class` per MAP code | 39 |
| `analyze_rescode.sh` | The extraction and classification pipeline | — |

### Comparison

| | Full dataset | STORM subset |
|---|---|---|
| Records | 256,453 | 169,603 |
| Converted to SMILES | 85,114 (33.2%) | 54,796 (32.3%) |
| ResCode occurrences | 225,757 | 178,889 |
| Distinct ResCodes | 35,711 | 11,683 |
| Distinct SkeletonCodes | 2,585 | 987 |
| **Distinct MAP codes** | **2,615** | **39** |

Two results stand out.

**The filter barely changes the SMILES conversion rate** — 32.3% against 33.2%. Passing the
STORM filter is not a predictor of being convertible to an atomic structure, because STORM
enables fuzzy matching and so admits the very patterns MolWURCS cannot expand: unknown anomeric
configuration, unknown carbon descriptors, unresolved linkage positions.

**The substituent variety collapses, from 2,615 MAP codes to 39.** The distinct ResCode count
falls by 67% and SkeletonCodes by 62%, but MAP codes fall by 98.5%. The surviving 39 are all
ordinary biological modifications:

| MAP code | Glycans | Modification |
|---|---|---|
| `*NCC/3=O` | 51,979 | N-acetyl |
| `*OCC/3=O` | 14,755 | O-acetyl |
| `*OC` | 7,980 | O-methyl |
| `*OSO/3=O/3=O` | 7,452 | Sulfate |
| `*N` | 5,722 | Amino |
| `*OPO/3O/3=O` | 2,084 | Phosphate |
| `*NCCO/3=O` | 1,446 | N-glycolyl |
| `*NSO/3=O/3=O` | 1,073 | N-sulfate |
| `*F` | 1,057 | Fluoro |
| `*OC^XO*/3CO/6=O/3C` | 88 | Pyruvate acetal (bridges two positions) |

The feature classification makes the same point from the other side. Every STORM MAP code falls
into just five classes:

| primary_class | Codes | Glycan occurrences |
|---|---|---|
| `acyl_short` | 13 | 68,623 |
| `simple` | 13 | 14,770 |
| `sulfate` | 6 | 8,954 |
| `phosphate` | 4 | 2,157 |
| `halogen` | 3 | 1,597 |

**`aromatic`, `lipid`, `macrocycle`, `alicyclic` and `multi_acyl` are all empty.** In the full
dataset those accounted for 1,419, 333, 112, 217 and 46 codes respectively — the benzyl and
benzoyl protecting groups, fatty acyl chains and macrolactones of section 8. The STORM
dictionary excludes them, which is the clearest single statement of what the filter does: it
separates naturally-occurring glycan chemistry from synthetic and exotic chemistry.

### Verification

| Check | Result |
|---|---|
| ID join between the filter output and `wurcs-smi-id.txt` | all 169,603 matched, 0 WURCS mismatches |
| GlyTouCan ID uniqueness (join safety) | 256,453 records, 256,453 distinct IDs, 0 duplicates |
| Extracted ResCodes vs WURCS header unique-counts | 178,889 = 178,889, 0 mismatched rows |
| SkeletonCode totals vs ResCode occurrences | 178,889 = 178,889 |
| MAP codes parsed by `MAPFactory` | all 39, zero errors |
| `analyze_rescode.sh` against the section 7–8 results | reproduces all six unfiltered files byte-identically |
| All 987 STORM SkeletonCodes resolved by `CarbonDescriptor` | zero unknown characters |
| `SkeletonFeatures.java` against known sugars — Glc, Gal, Man, Fuc, Neu5Ac, GlcA | carbon counts, terminal types, deoxy and acid counts all correct |
| `analyze_rescode.sh` with the SkeletonCode step added | reproduces all eight STORM files byte-identically |

---

## 10. SkeletonCode structural analysis

Section 9 showed that STORM cleans up the substituents (MAP codes) very effectively. This
section does the equivalent for the monosaccharide backbones, to find SkeletonCodes that remain
after filtering but are questionable as glycan residues.

A SkeletonCode is one character per backbone carbon, each a `CarbonDescriptor`. The same
character means different things at a terminal and a non-terminal position — `d` is `-CH2-` in
the middle of a chain while `m` is a terminal `-CH3` — so `SkeletonFeatures.java` resolves every
character through the framework's own `CarbonDescriptor.forCharacter(c, isTerminal)` rather than
matching on the text.

```bash
javac -cp molwurcs.jar SkeletonFeatures.java
cut -f1 storm-skeletoncode-unique.txt | java -cp molwurcs.jar:. SkeletonFeatures
```

Both are produced by `analyze_rescode.sh`:

| File | Description | Rows |
|---|---|---|
| `storm-skeletoncode-features.txt` | Per code: carbon count, first/last carbon type, counts of anomeric, defined/unknown stereocentres, deoxy, undefined, carbonyl, acid, double and triple bonds | 987 |
| `storm-skeletoncode-class.txt` | Feature flags and `primary_class` | 987 |
| `skeletoncode-features.txt`, `skeletoncode-class.txt` | The same for the unfiltered dataset | 2,585 |
| `SkeletonFeatures.java` | The extractor | — |

All 987 STORM SkeletonCodes resolved without an unknown character.

### Classes

| Class | Rule | Chemistry |
|---|---|---|
| `oversized` | more than 9 carbons | beyond the C3–C9 range of monosaccharides |
| `unsaturated` | a C=C in the backbone | e.g. Δ4,5-unsaturated uronic acid |
| `open_chain` | first carbon is an aldehyde | open-form aldose, no ring |
| `polyol` | no anomeric carbon, no carbonyl, terminal CH2OH | alditol / sugar alcohol |
| `undersized` | fewer than 5 carbons | triose, tetrose |
| `no_anomeric_other` | no anomeric carbon, not covered above | open-chain aldonic acid, open-form ketose |
| `standard` | none of the above | ring monosaccharide with an anomeric carbon |

### Distribution

| primary_class | Codes | Glycan occurrences |
|---|---|---|
| `standard` | 738 | 173,737 |
| `polyol` | 36 | 4,324 |
| `open_chain` | 109 | 294 |
| `unsaturated` | 4 | 203 |
| `no_anomeric_other` | 72 | 181 |
| `undersized` | 27 | 149 |
| `oversized` | 1 | 1 |

`standard` covers 97.1% of ResCode occurrences and spans C5–C9 (457 of its codes are hexoses,
103 are nonoses — the sialic acid family). Everything questionable is in the remaining 249
codes, which together account for 5,152 occurrences, under 3%.

### STORM already removed most of the anomalies

| Class | Full dataset | STORM | Removed |
|---|---|---|---|
| `standard` | 1,206 | 738 | 468 |
| `no_anomeric_other` | 408 | 72 | 336 |
| `unsaturated` | 266 | 4 | 262 |
| `undersized` | 209 | 27 | 182 |
| **`oversized`** | **165** | **1** | **164** |
| `polyol` | 154 | 36 | 118 |
| `open_chain` | 177 | 109 | 68 |

The unfiltered dataset also contains backbones with triple bonds (`11zz` and similar) and
undefined carbons (`u`, `U`, `Q`); **STORM leaves none of either**. The one over-long backbone it
does keep is `a212221122h`, an 11-carbon chain appearing in a single glycan.

### What is actually questionable

Not every class here is an error, and the largest ones are legitimate. Judgement is needed
before excluding anything:

**Legitimate — do not exclude.** The four `unsaturated` codes are `a21eEA`, `a11eEA`, `a12eEA`
and `a22eEA`: C6, anomeric, terminal COOH, with a double bond between C4 and C5. That is
Δ4,5-unsaturated uronic acid, the standard product of a glycosaminoglycan lyase digest, and its
presence is evidence of how the sample was prepared rather than a defect.

`open_chain` and `no_anomeric_other` are open-form representations of ordinary sugars —
`o2122h` is open-chain glucose, `A2122h` gluconic acid, `hO122h` an open-form ketose,
`AOd21122h` an open-chain ulosonic acid. They are correct descriptions of real molecules; they
simply are not ring forms.

**Judgement call.** The `polyol` class is the largest questionable group: `h2122h` (2,125
glycans) is glucitol, `h2112h` (1,349) galactitol, `h1122h` (216) mannitol. These are reduced
ends produced deliberately during glycan analysis. They are sugar alcohols, not monosaccharides,
so whether they belong depends on whether reduced glycans are in scope.

**Weakest as residues.** Within `polyol`, `h2h` and `hxh` (C3, 319 glycans between them) are
glycerol, and in `undersized` `hOh` is dihydroxyacetone — small polyols rather than sugars,
though glycerol is a genuine component of GPI anchors and glycoglycerolipids. The single
`oversized` code `a212221122h` is the clearest candidate for exclusion on structure alone.

### Suggested filter

If the goal is to keep only ring monosaccharides, the `standard` class is the criterion, and it
is one `awk` away:

```bash
awk -F'\t' 'NR>1 && $3=="standard" {print $1}' storm-skeletoncode-class.txt > keep-skeletons.txt
```

That retains 738 codes and 173,737 of 178,889 occurrences. Relaxing it to also keep
`unsaturated` — recommended, since those are real GAG residues — adds 4 codes and 203
occurrences. The raw measurements are in `storm-skeletoncode-features.txt`, so any other cut
(for example C5–C9 only, or requiring a defined anomeric configuration) can be made without
re-running the extractor.

---

## 11. Verification performed

| Check | Result |
|---|---|
| Row order of `wurcs-smi.txt` — column 1 vs cleaned `input.txt` column 1, all rows | `cmp` identical |
| ID/WURCS columns of `wurcs-smi-id.txt` vs `input.txt`, all rows | `cmp` identical |
| Row correspondence — SMILES converted back with `-i smi -o wurcs` must equal its own row's WURCS | 10 sampled rows, all match |
| ID correspondence — ID looked up in `input.txt` must give the same WURCS | 4 sampled rows (`G08770RO`, `G16155TD`, `G69397AC`, `G98994GL`), all match |
| Bracket parsing — number of `[...]` per WURCS vs the unique-residue count declared in the WURCS header | all 85,114 rows match |
| Total ResCodes extracted vs sum of declared unique-residue counts | 225,757 = 225,757 |
| ResCode decomposition — parts rejoined must reproduce the original string | all 35,711, zero failures |
| SkeletonCode totals vs total ResCode occurrences | 225,757 = 225,757 |
| Script output format, both modes, on a 3,000-row sample | 3 columns; `ONLY_SUCCESS=1` gave 1,613 rows with zero empty SMILES |
| MAP codes parsed by the framework's `MAPFactory` | all 2,615, zero parse errors |
| Feature extractor against known chemistry — N-acetyl, O-acetyl, sulfate, phosphate, benzyl, benzoyl, fluoro, amino, methyl, butoxy, N-glycolyl | 11 reference codes, all features correct |
| MAP column alignment between `mapcode-features.txt` and `mapcode-unique.txt` | `cmp` identical |
| Glycan-occurrence total in `mapcode-class.txt` vs `mapcode-unique.txt` | 181,037 = 181,037 |
| Lipid class plausibility — chain-length parity above C12 | 206 even vs 46 odd, consistent with fatty acyl chains |

---

## 12. Reproducing from scratch

```bash
./convert_wurcs.sh                                  # input.txt -> wurcs-smi-id.txt (~6 min)
awk -F'\t' '$3 != ""' wurcs-smi-id.txt > wurcs-smi-id-ok.txt
```

Then run the two `awk` blocks in section 7 to regenerate the ResCode, SkeletonCode and MAP
code files.
