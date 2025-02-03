Multi-sequence Alignment and Principal Coordinate Analysis Mapping Program, Forms Clusters of Similar DNA Sequences
FullListBH4.1.edit.strip.R is original function, others are all derivatives used to break up and process in smaller chunks when making maps


Code first runs an alignment on extracted FASTA file sequences of the desired length, then creates a distance matrix based on similarities between sequences using principal coordinate analysis, which is then mapped out.

Contains code for the full version, and a version where it splits up the full data set and then has the maps overlaid.

Does not include any FASTA sequences.
