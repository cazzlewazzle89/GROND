#!/usr/bin/env python3
import sys
import re

def main():
    if len(sys.argv) < 4:
        print("Usage: replace_headers.py <input_fna> <identifiers_txt> <output_fna>")
        sys.exit(1)
    input_fna, id_file, output_fna = sys.argv[1:4]

    # Build lookup set of valid OperonIDs from column 2 of the identifiers file.
    # bedtools getfasta -name writes the BED name column (OperonID) directly into
    # the FASTA header, optionally followed by a "(strand)" suffix when -s is used.
    valid_ids = set()
    with open(id_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                valid_ids.add(parts[1])

    strand_suffix = re.compile(r'\([+-]\)$')

    with open(input_fna, 'r') as infile, open(output_fna, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                raw = line[1:].strip()
                operon_id = strand_suffix.sub('', raw)
                if operon_id not in valid_ids:
                    print(f"Warning: '{operon_id}' not found in identifiers file.", file=sys.stderr)
                outfile.write(f'>{operon_id}\n')
            else:
                outfile.write(line)

if __name__ == '__main__':
    main()
