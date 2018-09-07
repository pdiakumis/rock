
# Rocking R at UMCCR

[![Travis build
status](https://travis-ci.org/umccr/rock.svg?branch=master)](https://travis-ci.org/umccr/rock)
[![Coverage
status](https://codecov.io/gh/umccr/rock/branch/master/graph/badge.svg)](https://codecov.io/github/umccr/rock?branch=master)

`rock` is an R package that (hopefully) helps with the day to day
bioinformatics life at UMCCR (UniMelb Centre for Cancer Research).

You can do the following:

  - Create circos plots using structural variant calls from
    [Manta](https://github.com/Illumina/manta), and copy number variant
    calls from [CNVkit](https://github.com/etal/cnvkit),
    [FACETS](https://github.com/mskcc/facets),
    [TitanCNA](https://github.com/gavinha/TitanCNA) or
    [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator).
    The OmicCircos R package is used. You can also use the Perl circos
    tool for Manta with CNVkit.

  - Create CNV profiles in horizontal facets for multiple samples or
    callers (piano plots). Can also zoom into specific chromosomes, and
    include an ideogram when specifying a single chromosome.

  - Generate bedgraph files for viewing the copy number segments in
    [IGV](http://software.broadinstitute.org/software/igv/) as a bar
    plot.

  - Generate IGV files for viewing SNP values in IGV as a scatter plot.

## Contents

<!-- vim-markdown-toc GFM -->

* [Installation](#installation)
* [Circos Plots](#circos-plots)
    * [Perl Circos](#perl-circos)
    * [OmicCircos](#omiccircos)
        * [Manta with CNVkit](#manta-with-cnvkit)
        * [Manta with PURPLE](#manta-with-purple)
* [Piano Plots](#piano-plots)
* [View CNV segments in IGV](#view-cnv-segments-in-igv)
* [View BED values in IGV](#view-bed-values-in-igv)

<!-- vim-markdown-toc -->

## Installation

You can install the development version of `rock` from
[GitHub](https://github.com/umccr/rock) with:

``` r
# install.packages("devtools") # if not pre-installed
devtools::install_github("umccr/rock") # master version
devtools::install_github("umccr/rock@v1.2.3") # release v1.2.3
devtools::install_github("umccr/rock@abcd") # commit abcd
```

There is no CRAN or conda version (yet).

Then just load with:

``` r
require(rock)
```

## Circos Plots

### Perl Circos

  - We can generate circos plots using the original
    [circos](http://circos.ca/) software package, written in Perl.
    **Note**: `circos` needs to be installed in your `PATH`.

  - Start by preparing the Manta and CNVkit calls. The required input
    files will be written to
`outdir`:

<!-- end list -->

``` r
manta <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
outdir <- "man/figures/perl_circos"
circos_prep(outdir = outdir, manta = manta, cnvkit = cnvkit)
#> Warning in dir.create(outdir, recursive = TRUE): 'man/figures/perl_circos'
#> already exists
#> Exporting Manta and CNVkit circos files to 'man/figures/perl_circos'.
#> Copying circos templates to 'man/figures/perl_circos'.
#> [1] TRUE
```

  - Then execute the following `circos` command either on the command
    line, or via the `plot_circos2`
function:

<!-- end list -->

``` bash
circos -nosvg -conf <outdir>/circos_simple.conf -outputdir <outdir> -outputfile foo_circos_cnvkit_manta.png
```

``` r
plot_circos2(outdir = outdir, name = "foo")
```

  - Result:

<!-- end list -->

``` r
knitr::include_graphics("man/figures/perl_circos/foo_circos_cnvkit_manta.png")
```

<img src="man/figures/perl_circos/foo_circos_cnvkit_manta.png" width="100%" />

### OmicCircos

  - We can generate circos plots using the functionality available in
    the
    [OmicCircos](https://bioconductor.org/packages/release/bioc/html/OmicCircos.html)
    Bioconductor R package.

  - Start by preparing the SV and CNV
calls.

<!-- end list -->

``` r
manta <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
facets <- system.file("extdata", "HCC2218_facets_cncf.tsv", package = "pebbles")
titan <- system.file("extdata", "HCC2218_titan.segs.tsv", package = "pebbles")
purple <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
truth <- system.file("extdata", "HCC2218_truthset_cnv_bcbio.tsv", package = "pebbles")

sv_manta <- prep_manta_vcf(manta)
cn_facets <- prep_facets_seg(facets)
cn_cnvkit <- prep_cnvkit_seg(cnvkit)
cn_purple <- prep_purple_seg(purple)
cn_truth <- prep_truth_seg(truth)
cn_titan <- prep_titan_seg(titan) # titan needs -1 for this case
cn_titan$cnv$tot_cn <- cn_titan$cnv$tot_cn - 1
```

  - Now we can generate a circos plot with Manta links and
    FACETS/CNVkit/TitanCNA/PURPLE segments (note that the `Warning`
    message below is due to a hack used in the OmicCircos code, where a
    matrix with numbers and characters is coerced to numeric. Just don’t
    worry about it..):

  - For the internal lines:
    
      - The \_inter\_chromosomal links take the chromosome colour of
        mate1 of each breakend pair.
      - The \_intra\_chromosomal lines are coloured according to the
        variant type:
          - Deletions: Red
          - Duplications: Green
          - Insertions: Purple
          - Inversions: Orange

#### Manta with CNVkit

``` r
plot_circos(sv = sv_manta, cnv = cn_cnvkit)
#> Warning in OmicCircos::circos(R = 260, cir = pebbles::circos_data$db, type
#> = "arc", : NAs introduced by coercion

#> Warning in OmicCircos::circos(R = 260, cir = pebbles::circos_data$db, type
#> = "arc", : NAs introduced by coercion
```

<img src="man/figures/README-circos-plot-manta-cnvkit-1.png" width="100%" />

#### Manta with PURPLE

``` r
plot_circos(sv = sv_manta, cnv = cn_purple)
#> Warning in OmicCircos::circos(R = 260, cir = pebbles::circos_data$db, type
#> = "arc", : NAs introduced by coercion

#> Warning in OmicCircos::circos(R = 260, cir = pebbles::circos_data$db, type
#> = "arc", : NAs introduced by coercion
```

<img src="man/figures/README-circos-plot-manta-purple-1.png" width="100%" />

## Piano Plots

  - We can generate ‘piano’ plots to compare CNV calls from multiple
    callers or
samples.

<!-- end list -->

``` r
cnv_list <- list(truth = cn_truth, cnvkit = cn_cnvkit, facets = cn_facets, purple = cn_purple, titan = cn_titan)
plot_piano(cnv_list = cnv_list)
```

<img src="man/figures/README-piano-plot-cnvkit-facets-purple-titan1-1.png" width="100%" />

  - You can also zoom into specific
chromosomes:

<!-- end list -->

``` r
plot_piano(cnv_list = cnv_list, chromosomes = c("1", "7", "8"), hide_x_lab = FALSE)
```

<img src="man/figures/README-piano-plot-chrom-1.png" width="100%" />

  - Change colours of the CNV segments:

<!-- end list -->

``` r
plot_piano(cnv_list = cnv_list, chromosomes = c("1", "7", "8"),
           seg.col = c("orange", "lightblue", "blue", "pink"), hide_x_lab = FALSE)
```

<img src="man/figures/README-piano-plot-chrom-colours-1.png" width="100%" />

  - And even plot an ideogram of the chromosome:

<!-- end list -->

``` r
require(patchwork)
#> Loading required package: patchwork
plot_ideogram(chrom = "13") +
  plot_piano(cnv_list = cnv_list,
             chromosomes = "13", hide_x_lab = FALSE) +
  plot_layout(ncol = 1, heights = c(1, 15))
```

<img src="man/figures/README-piano-plot-chrom-ideo-1.png" width="100%" />

## View CNV segments in IGV

``` r
cn_fname <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
cnv <- read_cnv(cn_fname)
cnv2igv(cnv, out_file = "~/Desktop/tmp/cnv_segs4igv.bedgraph", track_name = "cnv_segs2")
```

``` r
knitr::include_graphics("man/figures/README-cnv2igv_output.png")
```

<img src="man/figures/README-cnv2igv_output.png" width="100%" />

## View BED values in IGV

``` r
bed <- system.file("extdata", "HCC2218_baf.tsv", package = "pebbles")
bedval2igv(bed, out_file = "~/Desktop/tmp/baf1.igv", track_name = "baf", col = "purple")
```

``` r
# example for COLO829 whole-genome BAFs
knitr::include_graphics("man/figures/README-bedval2igv_output.png")
```

<img src="man/figures/README-bedval2igv_output.png" width="100%" />
