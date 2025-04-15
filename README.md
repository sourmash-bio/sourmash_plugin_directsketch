# sourmash_plugin_directsketch

[![PyPI](https://img.shields.io/pypi/v/sourmash_plugin_directsketch)](https://pypi.org/project/sourmash_plugin_directsketch/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sourmash_plugin_directsketch.svg)](https://anaconda.org/conda-forge/sourmash_plugin_directsketch)
[![DOI](https://zenodo.org/badge/792101561.svg)](https://zenodo.org/doi/10.5281/zenodo.11165725)


tl;dr - download and sketch data directly

## About

Commands:

- `gbsketch` - download and sketch NCBI Assembly Datasets by accession
- `urlsketch` - download and sketch directly from a url

This plugin is an attempt to improve sourmash database generation by downloading files, checking md5sum if provided, and sketching to a sourmash zipfile. FASTA/Q files can also be saved if desired. It's under active development.

## Installation

Option 1 (recommended): Create a conda environment and install into it:
```
conda create -n directsketch sourmash_plugin_directsketch # create and install
conda activate directsketch # activate
```
Option 2: Install without creating an environment

```
conda install sourmash_plugin_directsketch
```

## Usage Considerations

### Allowing restart with batching

If you're building large databases, we highly recommend you use batched zipfiles (v0.4+) to facilitate restart. If you encounter unexpected failures and are using a single zipfile output (default), `gbsketch`/`urlsketch` will have to re-download and re-sketch all files. If you instead set a batch size using `--batch-size`, then `gbsketch`/`urlsketch` can load any batched zips that finished writing, and avoid re-generating those signatures. The batch size represents the number of files downloaded, so it is possible DNA and protein signatures of the same accession may be split across zipfiles. Note that batches will use the `--output` file to build batched filenames, so if you provided `output.sig.zip`, batches will be `output.1.sig.zip`, etc (or `output.zip` --> `output.N.zip`). For small genomes (e.g. microbes), you can keep batch sizes quite large, e.g. 1000s-10000s. For large eukaryotic genomes where download takes much longer, you may want to use smaller batch sizes.

To build a single database after batched sketching, you can use `sig cat` to build a single zipfile (`sourmash sig cat *.zip -o OUTPUT.sig.zip`) or `sig collect` to collect all the zips into a standalone manifest that can be used with sourmash and branchwater commands.

### Rerunning failures (with batching or not)

**With batching**, you can use the same exact command to restart: directsketch will detect already-generated signatures and downloaded FASTAs and skip them. It will start writing signature from the next batch index, so if you last completed zipfile was `output.2.sig.zip`, running the command again with the same `--output output.sig.zip` will start from `output.3.sig.zip`. Incomplete zipfiles (marked as `.incomplete`) will be deleted.

**Without batching**, the existing signature zipfile will be unreadable and marked as `.incomplete`. You can restart by using the same input csv file, but if you are using `--keep-fasta`, you can add the option `--no-overwrite-fasta` to avoid overwriting any completed FASTA files. Incomplete fasta files are marked with the `.incomplete` extension and will be re-downloaded. If sketching, the incomplete zipfile will be deleted upon re-run, and we will still re-download all fasta files needed to generate sketches.

### Memory Requirements

Directsketch v0.6+ streams the downloaded data, sketching and/or writing as it goes. For gzipped files, the library we use checks the internal `crc32` to make sure we obtained the full download. `urlsketch` can also verify a user-provided `md5sum`. While you don't need to hold entire files in memory, **you do need enough memory to hold chunks of downloaded data and signatures while sketching**. You can limit the number of concurrent downloads (`--n-simultaneous-downloads`) to avoid memory issues. While testing with 10 eukaroyotic genomes under 1Gb each (10 simultaneous downloads), we used a maximum of ~2.5Gb.


### Using an NCBI API Key (gbsketch only)

`gbsketch` uses the NCBI REST API to download a dehydrated file with direct download links for the genomes. An API key may be needed if you run `gbsketch` many times, as NCBI has altered its download limitations. If you are unable to download the dehydrated zipfile, try providing an API Key. To obtain one, follow the instructions [here](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317). Once you have a key, you can provide it via the command line or set the `NCBI_API_KEY` variable (`export NCBI_API_KEY=YOUR_KEY`), which `gbsketch` will check and use automatically.

## Running the commands

## `gbsketch`
download and sketch NCBI Assembly Datasets by accession

### Create an input file

First, create a file, e.g. `acc.csv` with GenBank identifiers and sketch names.
```
accession,name
GCA_000961135.2,GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-45
GCA_000175555.1,GCA_000175555.1 ACUK01000506.1 Saccharolobus solfataricus 98/2
```
> Two columns must be present: `accession`, and `name`.

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
usage:  gbsketch [-h] [-q] [-d] [-o OUTPUT] [-f FASTAS] [--batch-size BATCH_SIZE] [-k] [--download-only] [--failed FAILED] [--checksum-fail CHECKSUM_FAIL] [-p PARAM_STRING] [-c CORES] [-r RETRY_TIMES]
                 [-n {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}] [-a API_KEY] [-v] [--write-urlsketch-csv] [--no-overwrite-fasta] [-g | -m]
                 input_csv

download and sketch GenBank assembly datasets

positional arguments:
  input_csv             A txt file or csv file containing accessions in the first column.

options:
  -h, --help            show this help message and exit
  -q, --quiet           suppress non-error output
  -d, --debug           provide debugging output
  -o OUTPUT, --output OUTPUT
                        Output ZIP file for the signatures. Must end with '.zip'.
  -f FASTAS, --fastas FASTAS
                        Write fastas here
  --batch-size BATCH_SIZE
                        Write smaller zipfiles, each containing sigs associated with this number of accessions. This allows gbsketch to recover after unexpected failures, rather than needing to restart sketching from scratch. Default: write all
                        sigs to single zipfile.
  -k, --keep-fasta      Write FASTA files. Default: do not write FASTA files.
  --download-only       Download FASTAS but do not sketch. Requires '--keep-fasta'. By default this downloads both genomes and proteomes.
  --failed FAILED       CSV of failed accessions and download links (should be mostly protein).
  --checksum-fail CHECKSUM_FAIL
                        CSV of accessions where the md5sum check failed or the md5sum file was improperly formatted or could not be downloaded.
  -p PARAM_STRING, --param-string PARAM_STRING
                        Parameter string for sketching (default: k=31,scaled=1000).
  -c CORES, --cores CORES
                        Number of cores to use (default is all available).
  -r RETRY_TIMES, --retry-times RETRY_TIMES
                        Number of times to retry failed downloads (default=3).
  -n {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}, --n-simultaneous-downloads {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
                        Number of files to download simultaneously (1-30; default=10).
  -a API_KEY, --api-key API_KEY
                        API Key for NCBI REST API. Alternatively, set NCBI_API_KEY environmental variable. If provided, will be used when downloading the initial dehyrated file.
  -v, --verbose         print progress for every download.
  --write-urlsketch-csv
                        Write urlsketch-formatted csv with all direct download links. Will be '{input_csv}.urlsketch.csv'.
  --no-overwrite-fasta  Requires `--keep-fasta`. If set, do not overwrite existing FASTA files in the --fastas directory. Will still re-download those files if needed for sketching.
  -g, --genomes-only    Download and sketch genome (DNA) files only.
  -m, --proteomes-only  Download and sketch proteome (protein) files only.
```

## `urlsketch`
download and sketch directly from URL(s)

### Create an input file

First, create a file, e.g. `acc-url.csv` with identifiers, sketch names, and other required info.
```
accession,name,moltype,md5sum,download_filename,url,range
GCA_000961135.2,GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-454,dna,47b9fb20c51f0552b87db5d44d5d4566,GCA_000961135.2_genomic.urlsketch.fna.gz,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/961/135/GCA_000961135.2_ASM96113v2/GCA_000961135.2_ASM96113v2_genomic.fna.gz,
GCA_000961135.2,GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-454,protein,fb7920fb8f3cf5d6ab9b6b754a5976a4,GCA_000961135.2_protein.urlsketch.faa.gz,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/961/135/GCA_000961135.2_ASM96113v2/GCA_000961135.2_ASM96113v2_protein.faa.gz,
GCA_000175535.1,GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14,dna,a1a8f1c6dc56999c73fe298871c963d1,GCA_000175535.1_genomic.urlsketch.fna.gz,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz,
```
> Six columns must be present:
> - `accession` - an accession or unique identifier. Ideally no spaces.
> - `name` - full name for the sketch. Sourmash convention is to use '{accession} rest of name' as the signature name.
> - `moltype` - is the file 'dna' or 'protein'?
> - `md5sum` - expected md5sum(s). Optional, will be checked after download if provided.
> - `download_filename` - filename for FASTA download. Required if `--keep-fastas`, but useful for signatures, too (saved in sig data).
> - `url` - direct link(s) for the file(s)
> - `range` - if desired, include base pair range(s), e.g 500-10000. This range will be selected from the record(s) and sketched (and/or saved to the download_filename). If there are multiple records in a FASTA file, the range will be applied to each record.

#### Note: Merging Files into the same signature
As of v0.5.0, `urlsketch` allows specification of multiple URLs to be downloaded and sketched into a single signature. If providing multiple URLs for a single accession/name, you must either provide no `md5sum` or `range`, or the number of entries in these columns must match the number of URLs. In each case, separate the entries with ';' -- e.g. "abc;def" for two md5sums.

### Run:

To run after creating file above:
```
sourmash scripts urlsketch acc-url.csv -o test-urlsketch.zip -f out_fastas -k --failed test.failed.csv -p dna,k=21,k=31,scaled=1000,abund -p protein,k=10,scaled=100,abund -r 1
```

Full Usage:
```
usage:  urlsketch [-h] [-q] [-d] [-o OUTPUT] [--batch-size BATCH_SIZE] [-f FASTAS] [-k] [--download-only] [--failed FAILED] [--checksum-fail CHECKSUM_FAIL] [-p PARAM_STRING] [-c CORES] [-r RETRY_TIMES]
                  [-n {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}] [--force] [-v] [--no-overwrite-fasta] [-g | -m]
                  input_csv

download and sketch GenBank assembly datasets

positional arguments:
  input_csv             A txt file or csv file containing accessions in the first column.

options:
  -h, --help            show this help message and exit
  -q, --quiet           suppress non-error output
  -d, --debug           provide debugging output
  -o OUTPUT, --output OUTPUT
                        Output ZIP file for the signatures. Must end with '.zip'.
  --batch-size BATCH_SIZE
                        Write smaller zipfiles, each containing sigs associated with this number of urls. This allows urlsketch to recover after unexpected failures, rather than needing to restart sketching from scratch. Default: write all sigs
                        to single zipfile.
  -f FASTAS, --fastas FASTAS
                        Write fastas here.
  -k, --keep-fasta, --keep-fastq
                        Write FASTA/Q files. Default: do not write FASTA files.
  --download-only       Download FASTAS but do not sketch. Requires '--keep-fasta/--keep-fastq'.
  --failed FAILED       CSV of failed accessions and download links.
  --checksum-fail CHECKSUM_FAIL
                        CSV of accessions where the md5sum check failed. If not provided, md5sum failures will be written to the download failures file (no additional md5sum information).
  -p PARAM_STRING, --param-string PARAM_STRING
                        Parameter string for sketching (default: k=31,scaled=1000).
  -c CORES, --cores CORES
                        Number of cores to use (default is all available).
  -r RETRY_TIMES, --retry-times RETRY_TIMES
                        Number of times to retry failed downloads (default=3).
  -n {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}, --n-simultaneous-downloads {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
                        Number of files to download simultaneously (1-30; default=10).
  --force               Skip input rows with empty or improper URLs. Warning: these will NOT be added to the failures file.
  -v, --verbose         print progress for every download.
  --no-overwrite-fasta  Requires `--keep-fasta`. If set, do not overwrite existing FASTA files in the --fastas directory. Will still re-download those files if needed for sketching.
  -g, --genomes-only    Download and sketch genome (DNA) files only.
  -m, --proteomes-only  Download and sketch proteome (protein) files only.
```

## Code of Conduct

This project is under the [sourmash Code of Conduct](https://github.com/sourmash-bio/sourmash/blob/latest/CODE_OF_CONDUCT.rst).

## Support

We suggest filing issues in [the directsketch issue tracker](https://github.com/sourmash-bio/sourmash_plugin_directsketch/issues) or [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues).

## Authors

* N. Tessa Pierce-Ward
* C. Titus Brown

## Dev docs

`sourmash_plugin_directsketch` is developed at https://github.com/sourmash-bio/sourmash_plugin_directsketch.

### Testing

Run:
```
pytest tests
```

### Generating a release

Bump version number in `Cargo.toml` and `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
make sdist
```
> Make sure directory is clean to avoid pulling in additional files.

followed by `make upload_sdist`.

> you may need to `pip install twine` if it is not available.
