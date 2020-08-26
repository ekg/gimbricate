# gimbricate

Corrects the overlaps in DNA sequence graphs by realigning sequence ends.

Optionally outputs FASTA and PAF representing the GFAv1-format graph.

## bad overlaps!

Almost no overlap-based assembly methods produce correct [GFAv1](https://github.com/GFA-spec/GFA-spec) output.
Invariably, some overlaps are incorrectly defined.
(The major exception to this are De Bruijn assemblers, which have fixed length overlaps that are correct by definition.)
In some methods (like [shasta](https://github.com/chanzuckerberg/shasta)), the overlaps are systematically slightly wrong, due to their derivation from run length encoded sequences.
It can help to correct these, because it lets us "bluntify" the graph.
This produces a graph in which each base in the graph exists on only one node, which is a desirable property for variation graph and other pangenomic reference models.

## installation

To build this fork of `gimbricate`, use git to download its source and cmake to build.
You'll need a recent gcc >= v7.4 as well as cmake >= v3.1

```
git clone --recursive https://github.com/eeg-ebe/gimbricate.git
cd gimbricate
cmake -H. -Bbuild && cmake --build build -- -j 4
```

By default, it builds into `bin/` in the directory root.

## Correcting the overlaps with `gimbricate`

gimbricate reads a GFA and uses the CIGAR strings attached to L records to guide realignment of the ends of the sequences and the subsequent correction of the CIGAR strings on the overlaps. It takes the sum of CIGAR operation lengths in the cigar to determine how much sequence to realign. Then it realigns the appropriate ends of the Sequence records referred to by the Link record with GSSW.

For instance, these link records might be approximately accurate:

L       256     +       1072    +       11M
L       1072    +       1074    +       20M
L       110     +       1072    -       20M
L       1072    -       1074    -       20M

But realignment suggests that only 11bp of the last one match:

L       256     +       1072    +       11M
L       1072    +       1074    +       20M
L       110     +       1072    -       20M
L       1072    -       1074    -       9D11M9I

## Pangenome construction

`gimbricate` can also convert its GFA output to a combination of FASTA and PAF.
These may then be fed into [`seqwish`](https://github.com/ekg/seqwish), which induces a variation graph including the original contigs as paths.
By adding additional alignments (in PAF format) to its input, `seqwish` can further compress redundant regions of the graph, but without loss of the original graph structure (which is given by the FASTA and PAF outputs of `gimbricate`).
By example (using test/h.gfa), this would bluntify the graph in h.gfa:

```
gimbricate -g h.gfa -n -p h.paf -f h.fasta >h.gimbry.gfa
seqwish -s h.fasta -p h.paf -g h.seqwish.gfa
```

Optional flag -P can be used when trying to bluntify a graph with already correct overlaps (such as those outputted by De Bruijn assemblers), skipping the time- and memory-consuming process of realignment.

This would additionally compress the graph by adding alignments between the nodes:

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
