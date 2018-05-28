require(devtools)
require(GenomeInfoDb)

chr_info <- GenomeInfoDb::genomeStyles("Homo_sapiens")

# aut, X, Y, M
chr_info$length <-  c(249250621L, 243199373L, 198022430L, 191154276L,
                      180915260L, 171115067L, 159138663L, 146364022L,
                      141213431L, 135534747L, 135006516L, 133851895L,
                      115169878L, 107349540L, 102531392L,  90354753L,
                      81195210L,  78077248L,   59128983L,  63025520L,
                      48129895L,  51304566L,  155270560L,  59373566L, 16569L)

devtools::use_data(chr_info)
