import os
import argparse
import sys

def create_symlinks(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    with open(input_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith("#") or not line:
                continue  # Skip comments and empty lines

            try:
                s= line.split(",")[0]  # Split line into sample name and file path
                f = line.split(",")[1]  # Split line into sample name and file path
                genome_path = os.path.join(output_dir, f"{s}.fa")

                if not os.path.exists(genome_path):  # Avoid redundant links
                    os.symlink(os.path.abspath(f), genome_path)  # Create symlink
                    print(f"Created symlink: {genome_path} -> {f}")
                else:
                    print(f"Symlink already exists: {genome_path}", file=sys.stderr)
            except ValueError:
                print(f"Skipping malformed line: {line}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create symbolic links from an input CSV file.")
    parser.add_argument("input", help="Path to the input CSV file")
    parser.add_argument("-o", "--output", default="Genomes", help="Output directory for symlinks (default: Genomes)")

    args = parser.parse_args()
    create_symlinks(args.input, args.output)
