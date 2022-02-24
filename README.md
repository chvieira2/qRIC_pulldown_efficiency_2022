# qRIC pulldown efficiency analysis code
<pre>
Cade associate with the publication in BioRXiv (https://www.biorxiv.org/content/10.1101/2021.07.12.452044v1.abstract)

To reproduce the main findings in the paper, simply download the "TXT" folders in the proteomics repository indicated in the paper and run this code in the order specified in the file name.

- 1_ProteinGroups.R will extract relevant protein level information and perform data filtering and manipulation for analysis. It also include some basic analysis.
- 2_PTMTable.R same as above, but for phosphorylation sites level information.
- 3_Peptides.R same as above, but for peptide level information.
- 4_workingtable.R will extract and aggregate the filtered protein, phosphorylation site and peptide level information into a unified table. This working table is further manipulated for calculation of the delta pull-down efficiencies.
- 5_PulldownEfficiency.R will extract the summarized information and perform some data analysis present in the final publication.
</pre>
