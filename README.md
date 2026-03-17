# virus_pipe

`virus_pipe` is a command-line workflow for virus discovery from sequencing data.

It supports input files in:
- `*.fastq`
- `*.fastq.gz`
- `*.fasta`

The pipeline orchestrates trimming, QC, assembly, mapping, unmapped read extraction, contig construction, `blastn`/`blastx`, and virus-focused reporting.

## Repository contents

- `virus_pipe.sh` – main analysis pipeline.
- `build_dictionary.sh` – helper script to build a BWA/Picard/SAMtools reference index.
- `virus_database.py` – adds ICTV/NCBI taxonomy information to candidate virus hits.
- `reverse.py` – helper used while generating adapter sequences.
- `ICTV_corrected.txt`, `ICTV_and_NCBI_Tax.txt` – taxonomy lookup tables used by `virus_database.py`.

## Requirements

The scripts assume these tools are available in your environment/path (or at the hard-coded locations in the scripts):

- Java
- Trimmomatic
- FastQC
- Trinity
- BWA
- SAMtools
- Picard
- CAP3
- BLAST+ (`blastn`, `blastx`)
- `fasta_formatter`
- Python (for `virus_database.py` and `reverse.py`)

> Note: `virus_pipe.sh` currently contains cluster-specific absolute paths (for example `/local/cluster/...`). Update those paths for your environment.

## 1) Build a reference dictionary

Before running the main pipeline, build a reference dictionary from a FASTA reference:

```bash
./build_dictionary.sh reference.fasta
```

This creates a directory named after `reference` with BWA, Picard dictionary (`.dict`), and `samtools faidx` outputs.

## 2) Configure `virus_pipe.sh`

In `virus_pipe.sh`, update the dictionary selection block to point to your reference index path.

Current logic expects a key and path placeholder:

```bash
if [ $2 = "genome_to_map_to" ]; then
    Dictionary="location"
fi
```

Replace:
- `genome_to_map_to` with the key you want to pass as the second argument.
- `location` with the actual reference index prefix path.

## 3) Run the pipeline

```bash
./virus_pipe.sh your_file.fastq genome_to_map_to
```

Where:
- `your_file.fastq` can also be `*.fastq.gz` or `*.fasta`.
- `genome_to_map_to` is the dictionary key configured in the script.

The script creates an output folder named after the input basename and writes intermediate files, BLAST reports, and a final `*_virus_report.txt` summary.

## Notes

- The default resource settings in `virus_pipe.sh` are:
  - `CPU=8`
  - `memory=64G`
- Intermediate files are compressed during processing to reduce disk usage.
- Legacy mention of `crawler.sh` has been removed from this README because that script is not present in this repository.
