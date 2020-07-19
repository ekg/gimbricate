## why this fork?

In this fork of [gimbricate](https://github.com/ekg/gimbricate/) we commented out the code for realigning overlaps (see this [discussion](https://github.com/ekg/gimbricate/issues/2)). Hence, this version of gimbricate should only be used when dealing with perfect overlaps, such as the ones produced by [Bwise](https://github.com/Malfoy/BWISE). Here the idea is to run first the present fork gimbricate on the GFA outputted by Bwise in order to generate a FASTA and a PAF, then run [seqwish](https://github.com/ekg/seqwish) on the FASTA and PAF in order to generate a GFA without any overlap. The original version of gimbricate seqfaulted on the large overlaps produced by Bwise, which is why we created this fork.
## installation

To build `gimbricate`, use git to download its source and cmake to build.
You'll need a recent gcc >= v7.4.

```
git clone --recursive https://github.com/ekg/gimbricate.git
cd gimbricate
cmake -H. -Bbuild && cmake --build build -- -j 4
```

By default, it builds into `bin/` in the directory root.

## correcting the overlaps with `gimbricate`

`gimbricate` reads a GFA and uses the CIGAR strings attached to `L` records to guide realignment of the ends of the sequences and the subsequent correction of the CIGAR strings on the overlaps.
It takes the sum of CIGAR operation lengths in the cigar to determine how much sequence to realign.
Then it realigns the appropriate ends of the `S`equence records referred to by the `L`ink record with [GSSW](https://github.com/vgteam/gssw).

For instance, these link records might be approximately accurate:

```
L       256     +       1072    +       11M
L       1072    +       1074    +       20M
L       110     +       1072    -       20M
L       1072    -       1074    -       20M
```

But realignment suggests that only 11bp of the last one match:

```
L       256     +       1072    +       11M
L       1072    +       1074    +       20M
L       110     +       1072    -       20M
L       1072    -       1074    -       9D11M9I
```

## pangenome construction

`gimbricate` can also convert its GFA output to a combination of FASTA and PAF.
These may then be fed into [`seqwish`](https://github.com/ekg/seqwish), which induces a variation graph including the original contigs as paths.
By adding additional alignments (in PAF format) to its input, `seqwish` can further compress redundant regions of the graph, but without loss of the original graph structure (which is given by the FASTA and PAF outputs of `gimbricate`).
By example (using test/h.gfa), this would bluntify the graph in h.gfa:

```
gimbricate -g h.gfa -n -p h.paf -f h.fasta >h.gimbry.gfa
seqwish -s h.fasta -p h.paf -g h.seqwish.gfa
```

While this would additionally compress it by adding alignments between the nodes:

```
gimbricate -g h.gfa -n -p h.paf -f h.fasta >h.gimbry.gfa
minimap2 -cx asm20 -k 11 -w 1 -X h.fasta h.fasta >h.self.paf
cat h.paf h.self.paf >h.self+.paf
seqwish -s h.fasta -p h.self+.paf -g h.seqwish+.gfa
```

Note that this techinque can be extended to provide graph to graph alignment and pangenome assembly.
All that we need to do is to concatenate together the gimbricate PAFs, FASTAs, and `minimap2` PAFs for the full collection of graphs, and then feed these to `seqwish`.

## segment name management

By default, `gimbricate` renames all input segments to numerical ids from 1 to N, where N = |Sequences|, but this can be disabled with `-n`.
This can simplify preprocessing when importing the graphs into `vg`, `odgi`, or other tools that expect numerical GFA names.

When building pangenomes from many inputs, it's often helpful to name input sequences uniquely.
Add a prefix to the names in the file with `-x` or `--name-prefix`.
This will propagate to the PAF, FASTA, and GFA outputs.
