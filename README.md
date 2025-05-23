## Table of Contents

   * [Introduction](#Introduction)
   * [Installation](#Installation)
   * [Subgenome phasing with WGDI](#Subgenome-phasing-with-WGDI)
   * [Subgenome phasing with SubPhaser](#Subgenome-phasing-with-SubPhaser)
   * [Citation](#Citation)

### Introduction ###
This is an example to phase subgenomes of an allopolyploid complex using 
[WGDI](https://github.com/SunPengChuan/wgdi) and [SubPhaser](https://github.com/zhangrengang/SubPhaser). 
Here we use the data of 
wheat complex (tetraploid–hexaploid reticulate allopolyploidization) as the example. 
The complex include an allotetraploid (AABB, *Triticum turgidum*, 2n = 4x = 28) and 
an allohexaploid (AABBDD, *T. aestivum*, 2n = 6x = 42).
Less than 0.8 million years ago (mya), a hybridization event between AA (T. urartu) and BB 
(a close relative of *Aegilops speltoides*) genomes
gave rise to the allopolyploid *T. turgidum* genome (AABB).
Subsequently, less than 0.4 mya, emmer wheat (AABB) hybridized with another wild wheat species
carrying the D genome (*A. tauschii*), resulting in the allohexaploid *T. aestivum* genome (AABBDD).
We assume that the diploid progenitors of allopolyploid wheats were either extinct or not sampled during the subgenome phasing process.

### Installation ###
Firstly, we need to install the required software ([WGDI](https://github.com/SunPengChuan/wgdi) and [SubPhaser](https://github.com/zhangrengang/SubPhaser)) for this example.
Here, we just install the software and dependencies via [conda](https://www.anaconda.com/).
```
git clone https://github.com/zhangrengang/SubPhaser
cd SubPhaser
conda env create -f SubPhaser.yaml -n SGphasing
conda activate SGphasing
python setup.py install

conda install -c bioconda wgdi diamond aster phytop newick_utils
```

### Subgenome phasing with WGDI ###
The rationales to assign subgenomes based on phylogenetic postions can be found in [the paper](https://doi.org/10.1093/bib/bbad513):
- The wheat complex (tetraploid–hexaploid reticulate allopolyploidization)
- The oat complex (tetraploid–hexaploid reticulate allopolyploidization)
- The poppy complex (tetraploid–octoploid reticulate allopolyploidization)
- Allotetraploids in the U’s triangle (tetraploid–tetraploid–tetraploid parallel allopolyploidization)
- The allooctoploid strawberry (single allooctoploid)

#### Prepare input data ####

1. Genomic data (protein sequences in fasta format and gene coordinates in custom gff format) of the allopolyploid complex are required.
2. Genomic data of potential diploid progenitors as far as possible are recommended (here is omitted).
3. Genomic data of outgroup, or [ancestral karyotype](https://github.com/SunPengChuan/wgdi-example/blob/main/Karyotype_Evolution.md) are required for phylogeny. 
One arbitrarily assigned subgenome (e.g. seven non-homoeologous chromosomes of allohexaploid wheat) of the ingroups can also be used as the reference for subgenome assignments (the outgroup is still required for rooting the phylogeny).
4. Configure files for WGDI.

Here, we just use the example data [*T. aestivum* (AABBDD) and *T. turgidum* (AABB), and the outgroup *Hordeum vulgare*] prepared in this repo:
```
git clone https://github.com/zhangrengang/subgenome_phasing_example
cd subgenome_phasing_example
cd wgdi
gunzip *gz
cat *.pep > pep.faa
cat *.cds > cds.fa
```
Now, all the required input data are present:
```
$ tree
├── ak.txt    # karyotype of the reference Hordeum_vulgare
├── Hordeum_vulgare.fasta
├── Hordeum_vulgare.gff
├── Hordeum_vulgare.lens
├── Triticum_aestivum.fasta
├── Triticum_aestivum.gff
├── Triticum_aestivum-Hordeum_vulgare.conf
├── Triticum_aestivum.lens
├── Triticum_aestivum-Triticum_aestivum.conf
├── Triticum_turgidum.fasta
├── Triticum_turgidum.gff
├── Triticum_turgidum-Hordeum_vulgare.conf
├── Triticum_turgidum.lens
├── Triticum_turgidum-Triticum_aestivum.conf
```
#### Run BLAST search ####
Blast results are also required for WGDI. Here, we run the pairwise BLAST search using DIAMOND:
```
diamond blastp -q Triticum_turgidum.pep -d Triticum_aestivum.pep -o Triticum_turgidum-Triticum_aestivum.blast --more-sensitive -p 40 --quiet -e 0.001
diamond blastp -q Triticum_turgidum.pep -d Hordeum_vulgare.pep -o Triticum_turgidum-Hordeum_vulgare.blast --more-sensitive -p 40 --quiet -e 0.001
diamond blastp -q Triticum_aestivum.pep -d Triticum_aestivum.pep -o Triticum_aestivum-Triticum_aestivum.blast --more-sensitive -p 40 --quiet -e 0.001
diamond blastp -q Triticum_aestivum.pep -d Hordeum_vulgare.pep -o Triticum_aestivum-Hordeum_vulgare.blast --more-sensitive -p 40 --quiet -e 0.001
```
These processes can be speed up by increasing `-p` or using parallel computation.

#### Detect synteny and calculate Ks ####
These are basic steps:
```
wgdi -icl Triticum_turgidum-Triticum_aestivum.conf
wgdi -ks Triticum_turgidum-Triticum_aestivum.conf
wgdi -bi Triticum_turgidum-Triticum_aestivum.conf

wgdi -icl Triticum_turgidum-Hordeum_vulgare.conf
wgdi -ks Triticum_turgidum-Hordeum_vulgare.conf
wgdi -bi Triticum_turgidum-Hordeum_vulgare.conf

wgdi -icl Triticum_aestivum-Triticum_aestivum.conf
wgdi -ks Triticum_aestivum-Triticum_aestivum.conf
wgdi -bi Triticum_aestivum-Triticum_aestivum.conf

wgdi -icl Triticum_aestivum-Hordeum_vulgare.conf
wgdi -ks Triticum_aestivum-Hordeum_vulgare.conf
wgdi -bi Triticum_aestivum-Hordeum_vulgare.conf
```

#### [Optional] Seek evidence from Ks-colored dot plots ####
To show Ks-colored dot plots:
```
wgdi -bk Triticum_turgidum-Triticum_aestivum.conf
wgdi -bk Triticum_aestivum-Triticum_aestivum.conf
```
From the resulted dot plots, 
we can find that the D subgenome of `Triticum_aestivum` shows 
higher Ks to `Triticum_turgidum`, while A or B subgenomes of `Triticum_aestivum` show lower Ks to `Triticum_turgidum` (Fig. 1 left). 
Thus, the D subgenome as a singleton can be phased out.

![Triticum_turgidum-Triticum_aestivum.blockks](wgdi/Triticum_turgidum-Triticum_aestivum.blockks.png) | ![Triticum_aestivum.blockks](wgdi/Triticum_aestivum-Triticum_aestivum.blockks.png)
---|---

**Fig. 1. Ks-colored dot plots between `Triticum_turgidum` and `Triticum_aestivum` (left) and within `Triticum_aestivum` (right).** The lower Ks indicates the higher similarity, and the lowest inter-genomic Ks (e.g. the red colored in the left panel) indicates orthology.

We hypothesize the A or B subgenome may be closer to the D subgenome. However, there is no such a pattern that Ks(A-D) is higher or lower than Ks(B-D) to distinguish A and B subgenomes (Fig. 1 right).

Alternatively, Ks can be replaced by [*Orthology index*](https://github.com/zhangrengang/orthoindex) 
which shows much more clear orthology relationshipes (Fig. 1c). 
It is highly recommended when Ks patterns are not clear:
```
soi dotplot -s Triticum_turgidum-Triticum_aestivum.collinearity -g *.gff@(.gz|) -c Triticum_turgidum-Triticum_aestivum.ctl \
     --ks-hist --max-ks 1 -o Triticum_turgidum-Triticum_aestivum.io \
     --plot-ploidy --gene-axis --xlabel '$Triticum~turgidum$' --ylabel '$Triticum~aestivum$' \
     --number-plots --ofdir ../OrthoFinder/OrthoFinder/Results_* --of-color

soi dotplot -s Triticum_aestivum-Triticum_aestivum.collinearity -g *.gff@(.gz|) -c Triticum_aestivum-Triticum_aestivum.ctl \
     --ks-hist --max-ks 1 -o Triticum_aestivum-Triticum_aestivum.io \
     --plot-ploidy --gene-axis --xlabel '$Triticum~aestivum$' --ylabel '$Triticum~aestivum$' \
     --number-plots --ofdir ../OrthoFinder/OrthoFinder/Results_* --of-color

```
[//]: ![Triticum_turgidum-Triticum_aestivum.orthoindex](wgdi/Triticum_turgidum-Triticum_aestivum.io.png)
<img src="wgdi/Triticum_turgidum-Triticum_aestivum.io.png" alt="Triticum_turgidum-Triticum_aestivum.orthoindex" width="400" >

**Fig. 1c. *Orthology index*-colored dot plots between `Triticum_turgidum` and `Triticum_aestivum`.** 
The highest inter-genomic *Orthology index* (i.e. the red colored dots) indicates orthology.


#### Assign subgenome preliminarily ####
Fisrt, we need to identify orthologous synteny between the outgroup reference and the polyploids, 
and to visually validate with the Ks-colored dot plots:
```
wgdi -c Triticum_turgidum-Hordeum_vulgare.conf
wgdi -bk Triticum_turgidum-Hordeum_vulgare.conf

wgdi -c Triticum_aestivum-Hordeum_vulgare.conf
wgdi -bk Triticum_aestivum-Hordeum_vulgare.conf
```
If there are non-orthologous syntenic blocks in the dot plots (Fig. 2), 
we need to adjust the parameters (including `homo`, `multiple` and `pvalue`) 
for `wgdi -c`, and re-run the above commands. Sometimes, we need to delete the 
out-paralogous blocks from the `blockinfo` file manually.

![Triticum_aestivum-Hordeum_vulgare.blockks](wgdi/Triticum_aestivum-Hordeum_vulgare.blockks.png) | ![Triticum_turgidum-Hordeum_vulgare.blockks](wgdi/Triticum_turgidum-Hordeum_vulgare.blockks.png)
---|---

**Fig. 2. Orthologous synteny.**

Alternatively, orthologous syntenic blocks can also be robustly identified via 
[*Orthology index*](https://github.com/zhangrengang/orthoindex)
by synthezising synteny with pre-inferred orthology 
(see [an example](https://github.com/zhangrengang/evolution_example)). 
It is highly recommended when `wgdi -c` do not work well:
```
mv Triticum_turgidum-Hordeum_vulgare.collinearity Triticum_turgidum-Hordeum_vulgare.collinearity.raw
soi filter -s Triticum_turgidum-Hordeum_vulgare.collinearity.raw -o ../OrthoFinder/OrthoFinder/Results_* -c 0.5 > Triticum_turgidum-Hordeum_vulgare.collinearity
wgdi -bi Triticum_turgidum-Hordeum_vulgare.conf
ln -f Triticum_turgidum-Hordeum_vulgare.blockinfo.csv Triticum_turgidum-Hordeum_vulgare.blockinfo.new.csv

mv Triticum_aestivum-Hordeum_vulgare.collinearity Triticum_aestivum-Hordeum_vulgare.collinearity.raw
soi filter -s Triticum_aestivum-Hordeum_vulgare.collinearity.raw -o ../OrthoFinder/OrthoFinder/Results_* -c 0.5 > Triticum_aestivum-Hordeum_vulgare.collinearity
wgdi -bi Triticum_aestivum-Hordeum_vulgare.conf
ln -f Triticum_aestivum-Hordeum_vulgare.blockinfo.csv Triticum_aestivum-Hordeum_vulgare.blockinfo.new.csv
```

Then, we map the karyotype of polyploids to the reference:
```
wgdi -km Triticum_turgidum-Hordeum_vulgare.conf
wgdi -km Triticum_aestivum-Hordeum_vulgare.conf
```
This step generates the karyotype mapping files of the two wheats: `Triticum_turgidum.ancestor.txt` and `Triticum_aestivum.ancestor.txt`.

At this stage, we need to manually edit the two files to assign subgenomes. 
We have phased the D subgenome based on Ks-based evidence, 
so we number the blocks of D subgenome as `3`, but have to randomly number those of A or B subgenomes as `1` or `2`. Meanwhile,
the assignments of A/B subgenomes of the two wheats are consistant based on the inter-genomic similarity/orthology. 
For eaxmple, if we assign `1A` as `2`, we need to assign `1B` as `1`, for `Triticum_aestivum`, based on their synteny; 
accordingly, we need to assign `1A` as `2` and `1B` as `1`, for `Triticum_turgidum`, based on their orthology to `Triticum_aestivum`: 
```
$ cat Triticum_aestivum.ancestor.txt
1A      1       4359    RoyalBlue       2
1B      1       4736    RoyalBlue       1
1D      1       4487    RoyalBlue       3
2A      1       5840    red     1
2B      1       6152    red     2
2D      1       5885    red     3
3A      1       5237    #99CC00 1
3B      1       5941    #99CC00 2
3D      1       5306    #99CC00 3
4A      1       3027    deepskyblue     1
4A      3028    3593    #339966 1
4A      3594    3907    fuchsia 2
4A      3908    4056    deepskyblue     1
4A      4057    4870    fuchsia 2
4B      1       3878    deepskyblue     2
4D      1       3582    deepskyblue     3
5A      1       4772    #339966 1
5A      4773    5450    deepskyblue     1
5B      1       5574    #339966 2
5D      1       5574    #339966 3
6A      1       4141    #FFCC00 2
6B      1       4627    #FFCC00 1
6D      1       4012    #FFCC00 3
7A      1       5573    fuchsia 1
7B      1       4892    fuchsia 2
7D      1       5419    fuchsia 3

$ cat Triticum_turgidum.ancestor.txt
1A      1       3906    RoyalBlue       2
1B      1       4136    RoyalBlue       1
2A      1       5192    red     1
2B      1       5463    red     2
3A      1       4956    #99CC00 1
3B      1       5832    #99CC00 2
4A      1       2869    deepskyblue     1
4A      2870    3347    #339966 1
4A      3348    3447    deepskyblue     1
4A      3448    4385    fuchsia 2
4B      1       3487    deepskyblue     2
5A      1       4057    #339966 1
5A      4058    4677    deepskyblue     1
5B      1       5037    #339966 2
6A      1       3651    #FFCC00 2
6B      1       4017    #FFCC00 1
7A      1       4880    fuchsia 1
7B      1       4244    fuchsia 2
```
The fragmented segments are assigned according to the complementarity of segments. For example, the large segment of chr4A-3' is assigned together with chr7B because they are complementary.

Now, we can apply the assignments (`wgdi -pc`) and generate alignments (`wgdi -a`):
```
wgdi -pc Triticum_turgidum-Hordeum_vulgare.conf
wgdi -a Triticum_turgidum-Hordeum_vulgare.conf

wgdi -pc Triticum_aestivum-Hordeum_vulgare.conf
wgdi -a Triticum_aestivum-Hordeum_vulgare.conf
```
![iteration1-Triticum_aestivum](wgdi/iteration1/Triticum_aestivum-Hordeum_vulgare.alignment.png) | ![iteration1/Triticum_turgidum](wgdi/iteration1/Triticum_turgidum-Hordeum_vulgare.alignment.png)
---|---

**Fig. 3. Subgenome assignments based on Ks evidence.** The same colored dot plots indicate the same subgenome assignments.

#### Reconstruct phylogeny by chromosomes and refine the assignments with the phylogeny-based evidence ####

We merge the alignments and build chromosome phylogeny to seek the phylogeny-based evidence (`wgdi -at`):
```
paste Triticum_turgidum-Hordeum_vulgare.alignment.csv Triticum_aestivum-Hordeum_vulgare.alignment.csv | perl -pe 's/\t[^,]+//g' > merged.alignment.csv

for chr in $(cut -f1 Hordeum_vulgare.lens)
do
    awk -v chr=$chr '$1==chr' Hordeum_vulgare.lens > Hordeum_vulgare.$chr.lens
    echo "[alignmenttrees]
alignment = merged.alignment.csv
gff = Hordeum_vulgare.gff
lens = Hordeum_vulgare.$chr.lens
dir = Hordeum_vulgare.$chr.tree
sequence_file = pep.faa
trees_file =  Hordeum_vulgare.$chr.trees.nwk
align_software = mafft
tree_software =  iqtree
model = MFP
trimming =  trimal
minimum = 4
delete_detail = true" > Hordeum_vulgare.$chr.conf
    wgdi -at Hordeum_vulgare.$chr.conf
    astral-pro --root 1 -i Hordeum_vulgare.$chr.trees.nwk -u 2 -t 8 -o Hordeum_vulgare.$chr.trees.nwk.astral
    phytop -pie -cp Hordeum_vulgare.$chr.trees.nwk.astral
    nw_topology -I Hordeum_vulgare.$chr.trees.nwk.astral | nw_order - | nw_display - -w 30
done
```
The following is the topology of these chromosome phylogenies:
```
chr1-3H:
 /-------------------------+ 1   /-------------------------+ 1   /-------------------------+ 1
 |                               |                               |
 |            /------------+ 2   |                  /------+ 2   |                  /------+ 2
=+      /-----+                 =+            /-----+           =+            /-----+
 |      |     \------------+ 4   |      /-----+     \------+ 4   |      /-----+     \------+ 4
 |      |                        |      |     |                  |      |     |
 \------+           /------+ 3   |      |     \------------+ 6   |      |     \------------+ 6
        |     /-----+            \------+                        \------+
        \-----+     \------+ 5          |     /------------+ 3          |     /------------+ 3
              |                         \-----+                         \-----+
              \------------+ 6                \------------+ 5                \------------+ 5

chr4-6H:
 /-------------------------+ 1   /-------------------------+ 1   /-------------------------+ 1
 |                               |                               |
 |                  /------+ 2   |                  /------+ 2   |            /------------+ 2
=+            /-----+           =+            /-----+           =+      /-----+
 |      /-----+     \------+ 4   |      /-----+     \------+ 4   |      |     \------------+ 4
 |      |     |                  |      |     |                  |      |
 |      |     \------------+ 6   |      |     \------------+ 6   \------+           /------+ 3
 \------+                        \------+                               |     /-----+
        |     /------------+ 3          |     /------------+ 3          \-----+     \------+ 5
        \-----+                         \-----+                               |
              \------------+ 5                \------------+ 5                \------------+ 6

chr7H:
 /-------------------------+ 1
 |
 |                  /------+ 2
=+            /-----+
 |      /-----+     \------+ 4
 |      |     |
 |      |     \------------+ 6
 \------+
        |     /------------+ 3
        \-----+
              \------------+ 5
```
The tip numbers in the trees are corresponding with the columns of `merged.alignment.csv`, so `1` = `Hordeum_vulgare`, 
`2-3` = `Triticum_turgidum 1-2`, `4-6` = `Triticum_aestivum 1-3`. 
The tip names can also be changed by:
```
for chr in $(cut -f1 Hordeum_vulgare.lens)
do
	nw_rename Hordeum_vulgare.$chr.trees.nwk.astral sg.idmap > Hordeum_vulgare.$chr.trees.nwk.astral.rename
done
```


We can find that all the topologies are identical (i.e. `((A, D), B)`) in fact, thus
we manually adjust the assignments according to the phylogenetic positions, by assigning chromosomes that is sister to `Triticum_aestivum 3` as `1` 
and assigning the sisters of `1+3` as `2`:
```
$ cat Triticum_aestivum.ancestor.txt
1A      1       4359    RoyalBlue       1
1B      1       4736    RoyalBlue       2
1D      1       4487    RoyalBlue       3
2A      1       5840    red     1
2B      1       6152    red     2
2D      1       5885    red     3
3A      1       5237    #99CC00 1
3B      1       5941    #99CC00 2
3D      1       5306    #99CC00 3
4A      1       3027    deepskyblue     1
4A      3028    3593    #339966 1
4A      3594    3907    fuchsia 2
4A      3908    4056    deepskyblue     1
4A      4057    4870    fuchsia 2
4B      1       3878    deepskyblue     2
4D      1       3582    deepskyblue     3
5A      1       4772    #339966 1
5A      4773    5450    deepskyblue     1
5B      1       5574    #339966 2
5D      1       5574    #339966 3
6A      1       4141    #FFCC00 1
6B      1       4627    #FFCC00 2
6D      1       4012    #FFCC00 3
7A      1       5573    fuchsia 1
7B      1       4892    fuchsia 2
7D      1       5419    fuchsia 3

$ cat Triticum_turgidum.ancestor.txt
1A      1       3906    RoyalBlue       1
1B      1       4136    RoyalBlue       2
2A      1       5192    red     1
2B      1       5463    red     2
3A      1       4956    #99CC00 1
3B      1       5832    #99CC00 2
4A      1       2869    deepskyblue     1
4A      2870    3347    #339966 1
4A      3348    3447    deepskyblue     1
4A      3448    4385    fuchsia 2
4B      1       3487    deepskyblue     2
5A      1       4057    #339966 1
5A      4058    4677    deepskyblue     1
5B      1       5037    #339966 2
6A      1       3651    #FFCC00 1
6B      1       4017    #FFCC00 2
7A      1       4880    fuchsia 1
7B      1       4244    fuchsia 2
```
We re-run the above commands (`wgdi -pc`, `-a`, `-at`):
![iteration2-Triticum_aestivum](wgdi/iteration2/Triticum_aestivum-Hordeum_vulgare.alignment.png) | ![iteration2/Triticum_turgidum](wgdi/iteration2/Triticum_turgidum-Hordeum_vulgare.alignment.png)
---|---

**Fig. 4. Refined subgenome assignments based on phylogeny evidence.**

Now, all the topologies are identical:
```
chr1-3H:
 /-------------------------+ 1   /-------------------------+ 1   /-------------------------+ 1
 |                               |                               |
 |                  /------+ 2   |                  /------+ 2   |                  /------+ 2
=+            /-----+           =+            /-----+           =+            /-----+
 |      /-----+     \------+ 4   |      /-----+     \------+ 4   |      /-----+     \------+ 4
 |      |     |                  |      |     |                  |      |     |
 |      |     \------------+ 6   |      |     \------------+ 6   |      |     \------------+ 6
 \------+                        \------+                        \------+
        |     /------------+ 3          |     /------------+ 3          |     /------------+ 3
        \-----+                         \-----+                         \-----+
              \------------+ 5                \------------+ 5                \------------+ 5

chr4-6H:
 /-------------------------+ 1   /-------------------------+ 1   /-------------------------+ 1
 |                               |                               |
 |                  /------+ 2   |                  /------+ 2   |                  /------+ 2
=+            /-----+           =+            /-----+           =+            /-----+
 |      /-----+     \------+ 4   |      /-----+     \------+ 4   |      /-----+     \------+ 4
 |      |     |                  |      |     |                  |      |     |
 |      |     \------------+ 6   |      |     \------------+ 6   |      |     \------------+ 6
 \------+                        \------+                        \------+
        |     /------------+ 3          |     /------------+ 3          |     /------------+ 3
        \-----+                         \-----+                         \-----+
              \------------+ 5                \------------+ 5                \------------+ 5

chr7H:
 /-------------------------+ 1
 |
 |                  /------+ 2
=+            /-----+
 |      /-----+     \------+ 4
 |      |     |
 |      |     \------------+ 6
 \------+
        |     /------------+ 3
        \-----+
              \------------+ 5
```
The above processes (`wgdi -pc`, `-a`, `-at`) may be iterated more than twice to generate such a consistant phylogeny.

#### [Optional] Seek evidence from biased fractionation ####
We analyze the gene retain (`wgdi -r`) patterns:
```
wgdi -r Triticum_aestivum-Hordeum_vulgare.conf
wgdi -r Triticum_turgidum-Hordeum_vulgare.conf
```
![Triticum_aestivum.retain](wgdi/Triticum_aestivum-Hordeum_vulgare.alignment.retain.png) | ![Triticum_turgidum.retain](wgdi/Triticum_turgidum-Hordeum_vulgare.alignment.retain.png)
---|---

**Fig. 5. Gene retain of subgenomes.**

However, biased fractionation patterns to distinguish subgenomes are not observed.

#### [Optional] Build subgenome phylogeny ####
Then we can build a final subgenome phylogeny:
```
cat Hordeum_vulgare.*.trees.nwk > Hordeum_vulgare.trees.nwk
astral-pro --root 1 -i Hordeum_vulgare.trees.nwk -u 2 -t 8 -o Hordeum_vulgare.trees.nwk.astral
phytop -pie -cp Hordeum_vulgare.trees.nwk.astral
nw_display Hordeum_vulgare.trees.nwk.astral
```
<img src="wgdi/Hordeum_vulgare.trees.nwk.astral.png" alt="subgenome phylogeny" width="500" >

**Fig. 6. Subgenome phylogeny from phased results of WGDI.** `1` = `Hordeum_vulgare`,
`2` = `Triticum_turgidum A`, `3` = `Triticum_turgidum B`, `4` = `Triticum_aestivum A`, 
`5` = `Triticum_aestivum B`, `6` = `Triticum_aestivum D`.

### Subgenome phasing with SubPhaser ###
#### Prepare input data ####

1. Genomic data (genome sequences in fasta format) of the allopolyploid complex are required.
2. Homoeologous relationships of chromosomes are required. These can be obtained from above synteny analyses or whole genome alignments..

Here, we just use the example data [*T. aestivum* (AABBDD) and *T. turgidum* (AABB)] prepared in this repo:
```
git clone https://github.com/zhangrengang/subgenome_phasing_example
cd subgenome_phasing_example
cd subphaser

# download genome sequences of Triticum_aestivum
wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v2.1/iwgsc_refseqv2.1_assembly.fa.zip -c && \
    unzip iwgsc_refseqv2.1_assembly.fa.zip && \
    mv iwgsc_refseqv2.1_assembly.fa Triticum_aestivum-genome.fasta && \
    gzip Triticum_aestivum-genome.fasta -f && \
    rm iwgsc_refseqv2.1_assembly.fa.zip

# download genome sequences of Triticum_turgidum
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/231/445/GCA_900231445.1_Svevo.v1/GCA_900231445.1_Svevo.v1_genomic.fna.gz -O Triticum_turgidum-genome.fasta.gz
```
Now, all the required input data are present:
```
$ tree
├── Triticum_aestivum-genome.fasta.gz
├── Triticum_aestivum-sg.config
├── Triticum_turgidum-genome.fasta.gz
└── Triticum_turgidum-sg.config
```

#### Run SubPhaser ####
```
subphaser -i Triticum_aestivum-genome.fasta.gz -c Triticum_aestivum-sg.config -pre Triticum_aestivum_
subphaser -i Triticum_turgidum-genome.fasta.gz -c Triticum_turgidum-sg.config -pre Triticum_turgidum_
```
Then we need to check whether the genomes have been well phased and whether the identified potential exchanges are confident: 

On the clustering heatmap (Fig. 6B) and PCA plot (Fig. 6C), a subgenome is defined as well-phased 
if it has clearly distinguishable patterns of both differential k-mers and homeologous chromosomes, 
indicating that each subgenome shares subgenome-specific features as expected.

SubPhaser outputs all windows whose enrichments do not match the subgenome assignments of 
their chromosome as potential exchanges, but further manual checks are still needed to 
identify them as bona fide exchanges. For example, at the 3’ end of wheat chr4A (Fig. 6D), 
the significant enrichments of subgenome B-specific k-mers are continuous (2nd from outer to inner circles), 
and the subgenome B-specific k-mers are as abundant as those on the chromosomes of subgenome B (5th circle) 
which contrasts to other subgenomes (4th and 6th circles). 
Combining all this information, we can confidently conclude that there has been an exchange, 
since the possibility of assembly errors has been ruled out with Hi-C data previously, etc. 
Since the distributions of subgenome-specific k-mers are usually uneven across the genome, 
inference should be careful and cautious to avoid type II errors.

![Triticum_aestivum](subphaser/Triticum_aestivum-merge_figures.png) | ![Triticum_turgidum](subphaser/Triticum_turgidum-merge_figures.png)
---|---

**Fig. 6. Subgenome assignments based on subgenome-specific kmers.**


#### [Optional] Convert to WGDI format and build subgenome phylogeny ####
Here for comparison purpose, we convert the output of SubPhaser to the format of WGDI, to build the subgenome phylogeny using the same method (`wgdi + astral-pro`).

Link files for WGDI:
```
ln ../wgdi/*gff ../wgdi/*lens ../wgdi/*Hordeum_vulgare.conf ../wgdi/ak.txt ../wgdi/*Hordeum_vulgare.blockinfo.new.csv ../wgdi/pep.faa .
```

Convert the output of SubPhaser, with discarding segments < 5 Mb:
```
cat Triticum_aestivum_phase-results/Triticum_aestivum_k15_q200_f2.bin.group | grep -v "#" | awk '$6>=5{print $1"\t"$2"\t"$3"\t"$4}' > Triticum_aestivum.sg.bed
python ../script/subphaser2wgdi.py Triticum_aestivum.sg.bed Triticum_aestivum.gff Triticum_aestivum.ancestor.txt

cat Triticum_turgidum_phase-results/Triticum_turgidum_k15_q200_f2.bin.group | grep -v "#" | awk '$6>=5{print $1"\t"$2"\t"$3"\t"$4}' > Triticum_turgidum.sg.bed
python ../script/subphaser2wgdi.py Triticum_turgidum.sg.bed Triticum_turgidum.gff Triticum_turgidum.ancestor.txt
```

Build the subgenome phylogeny using WGDI:
```
wgdi -pc Triticum_turgidum-Hordeum_vulgare.conf
wgdi -a Triticum_turgidum-Hordeum_vulgare.conf

wgdi -pc Triticum_aestivum-Hordeum_vulgare.conf
wgdi -a Triticum_aestivum-Hordeum_vulgare.conf

paste Triticum_turgidum-Hordeum_vulgare.alignment.csv Triticum_aestivum-Hordeum_vulgare.alignment.csv | perl -pe 's/\t[^,]+//g' > merged.alignment.csv

echo "[alignmenttrees]
alignment = merged.alignment.csv
gff = Hordeum_vulgare.gff
lens = Hordeum_vulgare.lens
dir = Hordeum_vulgare.tree
sequence_file = pep.faa
trees_file =  Hordeum_vulgare.trees.nwk
align_software = mafft
tree_software =  iqtree
model = MFP
trimming =  trimal
minimum = 4
delete_detail = true" > Hordeum_vulgare.conf

wgdi -at Hordeum_vulgare.conf
astral-pro --root 1 -i Hordeum_vulgare.trees.nwk -u 2 -t 8 -o Hordeum_vulgare.trees.nwk.astral
phytop -pie -cp Hordeum_vulgare.trees.nwk.astral
nw_display Hordeum_vulgare.trees.nwk.astral
```

<img src="subphaser/Hordeum_vulgare.trees.nwk.astral.png" alt="subgenome phylogeny" width="500" >

**Fig. 6. Subgenome phylogeny from phased results of SubPhaser.** `1` = `Hordeum_vulgare`,
`2` = `Triticum_turgidum A`, `3` = `Triticum_turgidum B`, `4` = `Triticum_aestivum A`,
`5` = `Triticum_aestivum B`, `6` = `Triticum_aestivum D`.

### Citation ###
If you use this pipeline, please cite:
> Zhang RG, Shang HY, Jia KH, Ma YP. Subgenome phasing for complex allopolyploidy: case-based benchmarking and recommendations [J]. *Brief. Bioinform.*, 2024, 25 (1): bbad513 [DOI:10.1093/bib/bbad513](http://doi.org/10.1093/bib/bbad513)

If you use the following tools, please cite:
* [WGDI](https://github.com/SunPengChuan/wgdi): Sun P, Jiao B, Yang Y et. al. WGDI: A user-friendly toolkit for evolutionary analyses of whole-genome duplications and ancestral karyotypes [J]. Mol. Plant., 2022, 15 (12): 1841–1851 [10.1016/j.molp.2022.10.018](http://doi.org/10.1016/j.molp.2022.10.018)
* [SubPhaser](https://github.com/zhangrengang/SubPhaser): Jia K, Wang Z, Wang L et. al. SubPhaser: a robust allopolyploid subgenome phasing method based on subgenome-specific k-mers [J]. New Phytol., 2022, 235 (2): 801–809 [10.1111/nph.18173](http://doi.org/10.1111/nph.18173)
* [ASTRAL](https://github.com/chaoszhang/ASTER): Zhang C, Mirarab S. ASTRAL-Pro 2: ultrafast species tree reconstruction from multi-copy gene family trees [J]. Bioinformatics, 2022, 38 (21): 4949–4950 [10.1093/bioinformatics/btac620](http://doi.org/10.1093/bioinformatics/btac620)
* [SOI](https://github.com/zhangrengang/SOI): Zhang RG, Shang HY, Milne RI et. al. SOI: robust identification of orthologous synteny with the Orthology Index and broad applications in evolutionary genomics [J]. Nucleic. Acids. Res., 2025, 53 (7): gkaf320 [10.1093/nar/gkaf320](https://doi.org/10.1093/nar/gkaf320)
* [phytop](https://github.com/zhangrengang/phytop): Shang H, Jia K, Li N et. al. Phytop: A tool for visualizing and recognizing signals of incomplete lineage sorting and hybridization using species trees output from ASTRAL [J]. Hortic. Res., 2024: uhae330 [10.1093/hr/uhae330](http://doi.org/10.1093/hr/uhae330)
* DIAMOND: Buchfink B, Xie C, Huson D H. Fast and sensitive protein alignment using DIAMOND [J]. Nat. Methods, 2015, 12 (1): 59–60 [10.1038/nmeth.3176](http://doi.org/10.1038/nmeth.3176)
* Newick utilities: Junier T, Zdobnov E M. The Newick utilities: high-throughput phylogenetic tree processing in the Unix shell [J]. Bioinformatics, 2010, 26 (13): 1669–1670 [10.1093/bioinformatics/btq243](http://doi.org/10.1093/bioinformatics/btq243)
