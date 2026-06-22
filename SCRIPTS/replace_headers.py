#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) < 4:
        print("Usage: replace_headers.py <input_fna> <identifiers_txt> <output_fna>")
        sys.exit(1)
    input_fna, id_file, output_fna = sys.argv[1:4]
    
    # Read IDs (column 2 of id_file)
    ids = []
    with open(id_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                ids.append(parts[1])
                
    # Read FASTA, replace headers
    idx = 0
    with open(input_fna, 'r') as infile, open(output_fna, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                if idx < len(ids):
                    new_header = f">{ids[idx]}"
                    idx += 1
                else:
                    new_header = line.strip()
                outfile.write(new_header + '\n')
            else:
                outfile.write(line)
                
if __name__ == '__main__':
    main()
