# sourmash_plugin_directsketch

[![PyPI](https://img.shields.io/pypi/v/sourmash_plugin_directsketch)](https://pypi.org/project/sourmash_plugin_directsketch/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sourmash_plugin_directsketch.svg)](https://anaconda.org/conda-forge/sourmash_plugin_directsketch)
[![DOI](https://zenodo.org/badge/792101561.svg)](https://zenodo.org/doi/10.5281/zenodo.11165725)


tl;dr - download and sketch data directly

## About

Commands:

- `gbsketch` - download and sketch NCBI Assembly Datasets by accession
- `urlsketch` - download and sketch directly from a url

This plugin is an attempt to improve sourmash database generation by downloading files, checking md5sum if provided or accessible, and sketching to a sourmash zipfile. FASTA files can also be saved if desired. It's quite fast, but still very much at alpha level. Here be dragons.

## Installation

### Linux

Option 1 (recommended): Create a conda environment and install into it:
```
conda create -n directsketch sourmash_plugin_directsketch # create and install
conda activate directsketch # activate
```
Option 2: Install without creating an environment

```
conda install sourmash_plugin_directsketch
```

### Other Platforms

On other platforms, you can create a conda environment with requirements like so:
```
curl -JLO https://raw.githubusercontent.com/sourmash-bio/sourmash_plugin_directsketch/main/environment.yml
conda env create -f environment.yml
```

then activate the environment and install `sourmash_plugin_directsketch` via `pip`:
```
conda activate directsketch
pip install sourmash_plugin_directsketch
```

## Usage Considerations

If you're building large databases (over 20k files), we highly recommend you use batched zipfiles (v0.4+) to facilitate restart. If you encounter unexpected failures and are using a single zipfile output (default), `gbsketch`/`urlsketch` will have to re-download and re-sketch all files. If you instead set a batch size using `--batch-size`, e.g. 10000, then `gbsketch`/`urlsketch` can load any batched zips that finished writing, and avoid re-generating those signatures. For `gbsketch`, the batch size represents the number of accessions included in each zip, with all signatures associated with an accession grouped within a single `zip`. For `urlsketch`, the batch size represents the number of total signatures included in each zip. Note that batches will use the `--output` file to build batched filenames, so if you provided `output.zip`, your batches will be `output.1.zip`, `output.2.zip`, etc.


## Running the commands

## `gbsketch`
download and sketch NCBI Assembly Datasets by accession

### Create an input file

First, create a file, e.g. `acc.csv` with GenBank identifiers and sketch names.
```
accession,name,ftp_path
GCA_000961135.2,GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-45,
GCA_000175555.1,GCA_000175555.1 ACUK01000506.1 Saccharolobus solfataricus 98/2,
```
> Three columns must be present: `accession`, `name`, and `ftp_path`. The `ftp_path` column can be empty (as above), but no additional columns may be present.

#### What is ftp_path?

If you do not provide an `ftp_path`, `gbsketch` will use the accession to find the `ftp_path` for you.

If you choose to provide it, `ftp_path` must be the `ftp_path` column from NCBI's assembly summary files.

For reference:

- example `ftp_path`: [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/600/915/GCA_036600915.1_ASM3660091v1](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/600/915/GCA_036600915.1_ASM3660091v1)
- bacteria assembly summary file: [https://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt](https://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt)

### Run:

To test `gbsketch`, you can download a csv file and run:
```
curl -JLO https://raw.githubusercontent.com/sourmash-bio/sourmash_plugin_directsketch/main/tests/test-data/acc.csv
sourmash scripts gbsketch acc.csv -o test-gbsketch.zip -f out_fastas -k --failed test.failed.csv --checksum-fail test.checksum-failed.csv -p dna,k=21,k=31,scaled=1000,abund -p protein,k=10,scaled=100,abund -r 1
```
To check that the `zip` was created properly, you can run:
```
sourmash sig summarize test-gbsketch.zip
```
and you should see the following as output:

```
** loading from 'test-gbsketch.zip'
path filetype: ZipFileLinearIndex
location: /path/to/your/test-gbsketch.zip
is database? yes
has manifest? yes
num signatures: 5
** examining manifest...
total hashes: 10815
summary of sketches:
   2 sketches with dna, k=21, scaled=1000, abund      2884 total hashes
   2 sketches with dna, k=31, scaled=1000, abund      2823 total hashes
   1 sketches with protein, k=10, scaled=100, abund   5108 total hashes
```

Full Usage:

```
usage:  gbsketch [-h] [-q] [-d] [-o OUTPUT] [-f FASTAS] [--batch-size BATCH_SIZE] [-k] [--download-only] --failed FAILED --checksum-fail CHECKSUM_FAIL [-p PARAM_STRING] [-c CORES]
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
  --batch-size BATCH_SIZE
                        Write smaller zipfiles, each containing sigs associated with this number of accessions. This allows gbsketch to recover after unexpected failures, rather than needing to
                        restart sketching from scratch. Default: write all sigs to single zipfile.
  -k, --keep-fasta      write FASTA files in addition to sketching. Default: do not write FASTA files
  --download-only       just download genomes; do not sketch
  --failed FAILED       csv of failed accessions and download links (should be mostly protein).
  --checksum-fail CHECKSUM_FAIL
                        csv of accessions where the md5sum check failed or the md5sum file was improperly formatted or could not be downloaded
  -p PARAM_STRING, --param-string PARAM_STRING
                        parameter string for sketching (default: k=31,scaled=1000)
  -c CORES, --cores CORES
                        number of cores to use (default is all available)
  -r RETRY_TIMES, --retry-times RETRY_TIMES
                        number of times to retry failed downloads
  -g, --genomes-only    just download and sketch genome (DNA) files
  -m, --proteomes-only  just download and sketch proteome (protein) files
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
usage:  urlsketch [-h] [-q] [-d] [-o OUTPUT] [--batch-size BATCH_SIZE] [-f FASTAS] [-k] [--download-only] --failed FAILED [--checksum-fail CHECKSUM_FAIL] [-p PARAM_STRING] [-c CORES]
                  [-r RETRY_TIMES]
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
  --batch-size BATCH_SIZE
                        Write smaller zipfiles, each containing sigs associated with this number of accessions. This allows urlsketch to recover after unexpected failures, rather than needing to
                        restart sketching from scratch. Default: write all sigs to single zipfile.
  -f FASTAS, --fastas FASTAS
                        Write fastas here
  -k, --keep-fasta, --keep-fastq
                        write FASTA/Q files in addition to sketching. Default: do not write FASTA files
  --download-only       just download genomes; do not sketch
  --failed FAILED       csv of failed accessions and download links.
  --checksum-fail CHECKSUM_FAIL
                        csv of accessions where the md5sum check failed. If not provided, md5sum failures will be written to the download failures file (no additional md5sum information).
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

Bump version number in `Cargo.toml` and push.

Make a new release on github.

Then pull, and:

```
make sdist
```

followed by `make upload_sdist`.

> you may need to `pip install twine` if it is not available.
