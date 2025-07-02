import re
import csv
import argparse

parser = argparse.ArgumentParser(description="Process input and generate output.")

parser.add_argument(
    "--input_file",
    type=argparse.FileType("r"),
    default="results/miranda_output.txt",
    help="raw miRanda output containing target predictions (default: results/miranda_output.txt)",
)
parser.add_argument(
    "--output_file",
    type=argparse.FileType("w"),
    default="results/miranda_output_parsed.csv",
    help="parsed and filtered results as CSV (default: results/miranda_output_parsed.csv)",
)

args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file

summary_data = []

with open(input_file, "r") as infile:
    for line in infile:
        line = line.strip()

        # miRanda target predictions start with '>>' followed by tab-separated fields
        if line.startswith(">>"):
            fields = line[2:].split("\t")
            if len(fields) >= 10:
                miRNA = fields[0]
                target = fields[1]
                max_score = float(fields[4])  # Alignment score (higher is better)
                max_energy = float(
                    fields[5]
                )  # Minimum free energy (more negative is better)
                len1 = int(fields[7])  # miRNA length
                len2 = int(fields[8])  # target length
                positions = fields[9].strip()  # Matching position string

                # Extract Ensembl Transcript ID without version
                target_parts = target.split("|")
                ensembl_transcript_id = (
                    target_parts[0].split(".")[0] if len(target_parts) > 0 else None
                )

                # Save parsed prediction to list
                summary_data.append(
                    {
                        "Seq1": miRNA,
                        "Seq2": ensembl_transcript_id,
                        "Max Score": max_score,
                        "Max Energy": max_energy,
                        "Len1": len1,
                        "Len2": len2,
                        "Positions": positions,
                    }
                )

# Write to CSV
with open(output_file, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=[
            "Seq1",
            "Seq2",
            "Max Score",
            "Max Energy",
            "Len1",
            "Len2",
            "Positions",
        ],
    )
    writer.writeheader()
    writer.writerows(summary_data)

print(f"✅ Parsed {len(summary_data)} predictions into {output_file}")
