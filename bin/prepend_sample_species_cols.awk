#!/usr/bin/awk -f
BEGIN{
    OFS="\t";
	# Remove */ from input, e.g mlst/
	split(ARGV[1],file,"/");
	# Parse out sample name and species
	split(file[2],fn, "_");
	species = fn[1]" "fn[2];
	sample_name = fn[3];
	# Number of header lines
	header = 1;
	# What should be prepended on the header line
	header_prepends = "sample\tspecies"
}
{
    if (NR <= header) { 
		print header_prepends, $0;
		next;
	}
	print sample_name, species, $0
}