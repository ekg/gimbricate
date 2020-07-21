## why this fork?

In this fork of [gimbricate](https://github.com/ekg/gimbricate/) we commented out the code for realigning overlaps (see this [discussion](https://github.com/ekg/gimbricate/issues/2)). Hence, this version of gimbricate should only be used when dealing with perfect overlaps, such as the ones produced by [Bwise](https://github.com/Malfoy/BWISE). Here the idea is to run first the present fork of gimbricate on the GFA outputted by Bwise in order to generate a FASTA and a PAF, then run [seqwish](https://github.com/ekg/seqwish) on the FASTA and PAF in order to generate a GFA without any overlap. The original version of gimbricate seqfaulted on the large overlaps produced by Bwise, which is why we created this fork.

## installation

To build this fork of `gimbricate`, use git to download its source and cmake to build.
You'll need a recent gcc >= v7.4 as well as cmake >= v3.1

```
git clone --recursive https://github.com/eeg-ebe/gimbricate.git
cd gimbricate
cmake -H. -Bbuild && cmake --build build -- -j 4
```

By default, it builds into `bin/` in the directory root.


## pangenome construction

`gimbricate` can also convert its GFA output to a combination of FASTA and PAF.
These may then be fed into [`seqwish`](https://github.com/ekg/seqwish), which induces a variation graph including the original contigs as paths.
By adding additional alignments (in PAF format) to its input, `seqwish` can further compress redundant regions of the graph, but without loss of the original graph structure (which is given by the FASTA and PAF outputs of `gimbricate`).
By example (using test/h.gfa), this would bluntify the graph in h.gfa:

```
gimbricate -g h.gfa -n -p h.paf -f h.fasta >h.gimbry.gfa
seqwish -s h.fasta -p h.paf -g h.seqwish.gfa
```


## segment name management

By default, `gimbricate` renames all input segments to numerical ids from 1 to N, where N = |Sequences|, but this can be disabled with `-n`.
This can simplify preprocessing when importing the graphs into `vg`, `odgi`, or other tools that expect numerical GFA names.

When building pangenomes from many inputs, it's often helpful to name input sequences uniquely.
Add a prefix to the names in the file with `-x` or `--name-prefix`.
This will propagate to the PAF, FASTA, and GFA outputs.
