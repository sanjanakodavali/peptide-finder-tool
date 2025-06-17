from flask import Flask, render_template, request, send_from_directory
from Bio import SeqIO
import os
import re
import csv
from collections import defaultdict

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = "uploads"
os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)

# Extract and highlight AQ...AQAQ patterns
def extract_and_highlight_all(protein_str):
    pattern = re.compile(r"(AQ)(.*?)(AQAQ)")
    matches = pattern.findall(protein_str)

    extracted = []
    highlighted = protein_str
    for m in matches:
        middle = m[1]
        extracted.append(middle)
        highlighted = re.sub(f"AQ{re.escape(middle)}AQAQ",
                             f"AQ<mark>{middle}</mark>AQAQ",
                             highlighted, count=1)
    return extracted, highlighted

# Translate DNA in 6 frames
def translate_6_frames(dna_seq):
    frames = []
    # Forward strand
    for i in range(3):
        frames.append((f"+{i+1}", dna_seq[i:].translate(to_stop=False)))
    # Reverse strand
    rev_seq = dna_seq.reverse_complement()
    for i in range(3):
        frames.append((f"-{i+1}", rev_seq[i:].translate(to_stop=False)))
    return frames

@app.route("/", methods=["GET", "POST"])
def index():
    matches = []

    if request.method == "POST":
        file = request.files["fastq_file"]
        file_path = os.path.join(app.config["UPLOAD_FOLDER"], file.filename)
        file.save(file_path)

        peptide_global_counts = defaultdict(int)

        for record in SeqIO.parse(file_path, "fastq"):
            dna_seq = record.seq
            all_frames = translate_6_frames(dna_seq)

            for frame_label, protein_seq in all_frames:
                protein_str = str(protein_seq)
                extracted, highlighted = extract_and_highlight_all(protein_str)

                if extracted:
                    match_counts = {}
                    for pep in extracted:
                        match_counts[pep] = match_counts.get(pep, 0) + 1
                        peptide_global_counts[pep] += 1

                    matches.append({
                        "id": record.id,
                        "frame": frame_label,
                        "highlighted_protein": highlighted,
                        "match_counts": match_counts
                    })

        # Write simplified output: Peptide, Count
        output_csv = os.path.join(app.config["UPLOAD_FOLDER"], "extracted_peptides.csv")
        with open(output_csv, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Peptide", "Count"])
            for pep, count in peptide_global_counts.items():
                writer.writerow([pep, count])

        record_count = len(matches)
        return render_template("index.html", matches=matches, record_count=record_count)

    return render_template("index.html", matches=None)

@app.route("/download/<path:filename>")
def download_file(filename):
    return send_from_directory("uploads", filename, as_attachment=True)

if __name__ == "__main__":
    print("🚀 Flask server is starting on http://127.0.0.1:5000")
    app.run(debug=True, port=5000)
