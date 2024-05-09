# sourmash_plugin_directsketch

[![PyPI](https://img.shields.io/pypi/v/sourmash_plugin_directsketch)](https://pypi.org/project/sourmash_plugin_directsketch/)


tl;dr - download and sketch NCBI Assembly Datasets by accession

## About

This plugin is an attempt to improve database generation by downloading assemblies, checking md5sum, and sketching to a sourmash zipfile. FASTA files can also be saved if desired. It's quite fast, but still very much at alpha level. Here be dragons.

## Installation

```
pip install sourmash_plugin_directsketch
```

## Usage
First, create a file, e.g. `acc.csv` with GenBank identifiers and sketch names.
```
ident,name
GCA_000961135.2,GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-45
GCA_000175555.1,GCA_000175555.1 ACUK01000506.1 Saccharolobus solfataricus 98/2
```
> Extra columns are ok, as long as the first two columns contain the identifier and sketch name

Run:

To run the test accession file at `tests/test-data/acc.csv`, run:
```
sourmash scripts gbsketch tests/test-data/acc.csv -o test.zip -f out_fastas -k --failed test.failed.csv -p dna,k=21,k=31,scaled=1000,abund -p protein,k=10,scaled=100,abund -r 1
```

Full Usage:

```
usage:  gbsketch [-h] [-q] [-d] -o OUTPUT [-f FASTAS] [-k] [--download-only] [--failed FAILED] [-p PARAM_STRING] [-c CORES]
                 [-r RETRY_TIMES] [-g | -m]
                 input_csv

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
  -k, --keep-fastas     write FASTA files in addition to sketching. Default: do not write FASTA files
  --download-only       just download genomes; do not sketch
  --failed FAILED       csv of failed accessions and download links (should be mostly protein).
  -p PARAM_STRING, --param-string PARAM_STRING
                        parameter string for sketching (default: k=31,scaled=1000)
  -c CORES, --cores CORES
                        number of cores to use (default is all available)
  -r RETRY_TIMES, --retry-times RETRY_TIMES
                        number of times to retry failed downloads
  -g, --genomes-only    just download and sketch genome (DNA) files
  -m, --proteomes-only  just download and sketch proteome (protein) files
```

## Code of Conduct

This project is under the [sourmash Code of Conduct](https://github.com/sourmash-bio/sourmash/blob/latest/CODE_OF_CONDUCT.rst).

## Support

We suggest filing issues in [the directsketch issue tracker](https://github.com/bluegenes/sourmash_plugin_directsketch/issues) or [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues).

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
