# sourmash_plugin_directsketch

[![PyPI](https://img.shields.io/pypi/v/sourmash_plugin_directsketch)](https://pypi.org/project/sourmash_plugin_directsketch/)
[![DOI](https://zenodo.org/badge/792101561.svg)](https://zenodo.org/doi/10.5281/zenodo.11165725)


tl;dr - download and sketch data directly

## About

Commands:

- `gbsketch` - download and sketch NCBI Assembly Datasets by accession
- `urlsketch` - download and sketch directly from a url

This plugin is an attempt to improve sourmash database generation by downloading files, checking md5sum if provided or accessible, and sketching to a sourmash zipfile. FASTA files can also be saved if desired. It's quite fast, but still very much at alpha level. Here be dragons.

## Installation

```
pip install sourmash_plugin_directsketch
```

## `gbsketch`
download and sketch NCBI Assembly Datasets by accession

### Create an input file

First, create a file, e.g. `acc.csv` with GenBank identifiers and sketch names.
```
accession,name,ftp_path
GCA_000961135.2,GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-45,
GCA_000175555.1,GCA_000175555.1 ACUK01000506.1 Saccharolobus solfataricus 98/2,
```
> Three columns must be present: `accession`, `name`, and `ftp_path`. The `ftp_path` column can be empty, but no additional columns may be present.

#### What is ftp_path?

If you do not provide an `ftp_path`, `gbsketch` will use the accession to find the `ftp_path` for you.

If you choose to provide it, `ftp_path` must be the `ftp_path` column from NCBI's assembly summary files.

For reference:

- example `ftp_path`: [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/600/915/GCA_036600915.1_ASM3660091v1](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/600/915/GCA_036600915.1_ASM3660091v1)
- bacteria assembly summary file: [https://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt](https://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt)

### Run:

To run the test accession file at `tests/test-data/acc.csv`, run:
```
sourmash scripts gbsketch tests/test-data/acc.csv -o test-gbsketch.zip -f out_fastas -k --failed test.failed.csv -p dna,k=21,k=31,scaled=1000,abund -p protein,k=10,scaled=100,abund -r 1
```

Full Usage:

```
usage:  gbsketch [-h] [-q] [-d] [-o OUTPUT] [-f FASTAS] [-k] [--download-only] [--failed FAILED] [-p PARAM_STRING] [-c CORES] [-r RETRY_TIMES] [-g | -m] input_csv

download and sketch GenBank assembly datasets

positional arguments:
  input_csv             a txt file or csv file containing accessions in the first column

options:
  -h, --help            show this help message and exit
  -q, --quiet           suppress non-error output
  -d, --debug           provide debugging output
  -o OUTPUT, --output OUTPUT
                        output zip file for the signatures
  -f FASTAS, --fastas FASTAS
                        Write fastas here
  -k, --keep-fasta      write FASTA files in addition to sketching. Default: do not write FASTA files
  --download-only       just download genomes; do not sketch
  --failed FAILED       csv of failed accessions and download links (should be mostly protein).
  -p PARAM_STRING, --param-string PARAM_STRING
                        parameter string for sketching (default: k=31,scaled=1000)
  -c CORES, --cores CORES
                        number of cores to use (default is all available)
  -r RETRY_TIMES, --retry-times RETRY_TIMES
                        number of times to retry failed downloads
  -g, --genomes-only    just download and sketch genome (DNA) files
```

## `urlsketch`
download and sketch directly from a url
### Create an input file

First, create a file, e.g. `acc-url.csv` with identifiers, sketch names, and other required info.
```
accession,name,moltype,md5sum,download_filename,url
GCA_000961135.2,GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-454,dna,47b9fb20c51f0552b87db5d44d5d4566,GCA_000961135.2_genomic.urlsketch.fna.gz,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/961/135/GCA_000961135.2_ASM96113v2/GCA_000961135.2_ASM96113v2_genomic.fna.gz
GCA_000961135.2,GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-454,protein,fb7920fb8f3cf5d6ab9b6b754a5976a4,GCA_000961135.2_protein.urlsketch.faa.gz,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/961/135/GCA_000961135.2_ASM96113v2/GCA_000961135.2_ASM96113v2_protein.faa.gz
GCA_000175535.1,GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14,dna,a1a8f1c6dc56999c73fe298871c963d1,GCA_000175535.1_genomic.urlsketch.fna.gz,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz
```
> Six columns must be present:
> - `accession` - an accession or unique identifier. Ideally no spaces.
> - `name` - full name for the sketch.
> - `moltype` - is the file 'dna' or 'protein'?
> - `md5sum` - expected md5sum (optional, will be checked after download if provided)
> - `download_filename` - filename for FASTA download. Required if `--keep-fastas`, but useful for signatures, too (saved in sig data).
> - `url` - direct link for the file

### Run:

To run the test accession file at `tests/test-data/acc-url.csv`, run:
```
sourmash scripts urlsketch tests/test-data/acc-url.csv -o test-urlsketch.zip -f out_fastas -k --failed test.failed.csv -p dna,k=21,k=31,scaled=1000,abund -p protein,k=10,scaled=100,abund -r 1
```

Full Usage:
```
usage:  urlsketch [-h] [-q] [-d] [-o OUTPUT] [-f FASTAS] [-k] [--download-only] [--failed FAILED] [-p PARAM_STRING] [-c CORES] [-r RETRY_TIMES] input_csv

download and sketch GenBank assembly datasets

positional arguments:
  input_csv             a txt file or csv file containing accessions in the first column

options:
  -h, --help            show this help message and exit
  -q, --quiet           suppress non-error output
  -d, --debug           provide debugging output
  -o OUTPUT, --output OUTPUT
                        output zip file for the signatures
  -f FASTAS, --fastas FASTAS
                        Write fastas here
  -k, --keep-fasta, --keep-fastq
                        write FASTA/Q files in addition to sketching. Default: do not write FASTA files
  --download-only       just download genomes; do not sketch
  --failed FAILED       csv of failed accessions and download links (should be mostly protein).
  -p PARAM_STRING, --param-string PARAM_STRING
                        parameter string for sketching (default: k=31,scaled=1000)
  -c CORES, --cores CORES
                        number of cores to use (default is all available)
  -r RETRY_TIMES, --retry-times RETRY_TIMES
                        number of times to retry failed downloads
```

## Code of Conduct

This project is under the [sourmash Code of Conduct](https://github.com/sourmash-bio/sourmash/blob/latest/CODE_OF_CONDUCT.rst).

## Support

We suggest filing issues in [the directsketch issue tracker](https://github.com/sourmash-bio/sourmash_plugin_directsketch/issues) or [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues).

## Authors

* N. Tessa Pierce-Ward

## Dev docs

`sourmash_plugin_directsketch` is developed at https://github.com/sourmash-bio/sourmash_plugin_directsketch.

### Testing

Run:
```
pytest tests
```

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
python -m build
```

followed by `twine upload dist/...`.
