#!/usr/bin/env python3
import pathlib

out_dir = pathlib.Path("/home/vanbruggenmit/mit-ihh-pib/data/grch38")
out_dir.mkdir(parents=True, exist_ok=True)
out_file = out_dir / "contig_map_chr_to_plain.tsv"

mapping = (
    [(f"chr{i}", str(i)) for i in range(1, 23)]
    + [("chrX", "X")]  # add ("chrY","Y"), ("chrM","MT") if needed
)

with out_file.open("w") as f:
    for a, b in mapping:
        f.write(f"{a}\t{b}\n")
    f.write("# chrY\tY\n# chrM\tMT\n")

print(f"Wrote {out_file}")
