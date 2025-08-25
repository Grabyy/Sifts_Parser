# Sifts_Parser

**Sifts_Parser** is a simple Python script to read and process **SIFTS XML files**.  
It was originally developed to compare sequences between **UniProt** and **PDB** entries,  
with the goal of building datasets for predicting **3D peptide structures**.

---

## Features

- Parse **SIFTS XML** files directly from the [EBI SIFTS database](https://www.ebi.ac.uk/pdbe/docs/sifts/).
- Save information in a readable tsv file that could be easily processed.
- Uses only **standard Python 3 libraries** (no external dependencies).
- Process files directly from `curl`/`gunzip` pipelines — no need to save intermediate files.

---

## Installation

Only **Python 3**, `curl`, and `gunzip` are required.  
You can set up a dedicated conda environment:

```bash
conda create -n sifts_parser python -y
conda activate sifts_parser
```

---

## Usage

The idea inspired to &Github& was to curl the xml without saving it, then extract it and feed it to my script :

```bash
curl -s https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/1xyz.xml.gz | gunzip | python sifts_parse.py
```
---

## License

This project is released under the MIT License.
Feel free to use, modify, and distribute.
