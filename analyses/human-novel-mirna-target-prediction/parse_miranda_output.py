import re
import csv

input_file = "results/miranda_output.txt"
output_file = "results/miRNA_target_anno_v2.csv"

summary_data = []

with open(input_file, "r") as infile:
    for line in infile:
        line = line.strip()
        if line.startswith(">>"):
            fields = line[2:].split("\t")
            if len(fields) >= 10:
                miRNA = fields[0]
                target = fields[1]
                max_score = float(fields[4])
                max_energy = float(fields[5])
                len1 = int(fields[7])  # miRNA length
                len2 = int(fields[8])  # target length
                positions = fields[9].strip()

                # Extract Ensembl Transcript ID and Gene ID without version
                target_parts = target.split("|")
                ensembl_transcript_id = target_parts[0].split(".")[0] if len(target_parts) > 0 else None
                ensembl_gene_id = target_parts[1].split(".")[0] if len(target_parts) > 1 else None

                summary_data.append({
                    "Seq1": miRNA,
                    "Seq2": ensembl_transcript_id,
                    "Max Score": max_score,
                    "Max Energy": max_energy,
                    "Len1": len1,
                    "Len2": len2,
                    "Positions": positions,
                    "GeneID": ensembl_gene_id
                })

# Write to CSV
with open(output_file, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=[
        "Seq1", "Seq2", "Max Score", "Max Energy",
        "Len1", "Len2", "Positions", "GeneID"
    ])
    writer.writeheader()
    writer.writerows(summary_data)

print(f"✅ Parsed {len(summary_data)} predictions into {output_file}")

