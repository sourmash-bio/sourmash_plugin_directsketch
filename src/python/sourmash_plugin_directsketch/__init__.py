#! /usr/bin/env python
import os
import sys
from sourmash.logging import notify
from sourmash.plugins import CommandLinePlugin
import importlib.metadata
import argparse

from . import sourmash_plugin_directsketch

__version__ = importlib.metadata.version("sourmash_plugin_directsketch")
def print_version():
    notify(f"=> sourmash_plugin_directsketch {__version__}")

def get_max_cores():
    try:
        if 'SLURM_CPUS_ON_NODE' in os.environ:
            return int(os.environ['SLURM_CPUS_ON_NODE'])
        elif 'SLURM_JOB_CPUS_PER_NODE' in os.environ:
            cpus_per_node_str = os.environ['SLURM_JOB_CPUS_PER_NODE']
            return int(cpus_per_node_str.split('x')[0])
        else:
            return os.cpu_count()
    except Exception:
        return os.cpu_count()


def set_thread_pool(user_cores):
    avail_threads = get_max_cores()  # Define how to get the maximum available cores
    num_threads = min(avail_threads, user_cores) if user_cores else avail_threads
    if user_cores and user_cores > avail_threads:
        notify(f"warning: only {avail_threads} threads available, using {avail_threads}")
    actual_tokio_cores = sourmash_plugin_directsketch.set_tokio_thread_pool(num_threads)
    return actual_tokio_cores

def non_negative_int(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError(f"Batch size cannot be negative (input value: {value})")
    return ivalue

class Download_and_Sketch_Assemblies(CommandLinePlugin):
    command = 'gbsketch'
    description = 'download and sketch GenBank assembly datasets'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('input_csv', help="a txt file or csv file containing accessions in the first column")
        p.add_argument('-o', '--output', default=None,
                       help="output zip file for the signatures. Must end with '.zip'")
        p.add_argument('-f', '--fastas',
                       help='Write fastas here', default = '.')
        p.add_argument('--batch-size', type=non_negative_int, default = 0,
                       help='Write smaller zipfiles, each containing sigs associated with this number of accessions. \
                            This allows gbsketch to recover after unexpected failures, rather than needing to \
                            restart sketching from scratch. Default: write all sigs to single zipfile.')
        p.add_argument('-k', '--keep-fasta', action='store_true',
                       help="write FASTA files in addition to sketching. Default: do not write FASTA files")
        p.add_argument('--download-only', help='just download genomes; do not sketch', action='store_true')
        p.add_argument('--failed', help='csv of failed accessions and download links (should be mostly protein).', required=True)
        p.add_argument('--checksum-fail', help="csv of accessions where the md5sum check failed or the md5sum file was improperly formatted or could not be downloaded", required=True)
        p.add_argument('-p', '--param-string', action='append', type=str, default=[],
                          help='parameter string for sketching (default: k=31,scaled=1000)')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')
        p.add_argument('-r', '--retry-times', default=1, type=int,
                       help='number of times to retry failed downloads')
        group = p.add_mutually_exclusive_group()
        group.add_argument('-g', '--genomes-only', action='store_true', help='just download and sketch genome (DNA) files')
        group.add_argument('-m', '--proteomes-only', action='store_true', help='just download and sketch proteome (protein) files')


    def main(self, args):
        print_version()
        if not args.param_string:
            args.param_string = ["k=31,scaled=1000"]
        notify(f"params: {args.param_string}")

        if args.download_only and not args.keep_fasta:
            notify("Error: '--download-only' requires '--keep-fasta'.")
            sys.exit(-1)
        if args.output is None and not args.download_only:
            notify("Error: output signature zipfile is required if not using '--download-only'.")
            sys.exit(-1)

        # convert to a single string for easier rust handling
        args.param_string = "_".join(args.param_string)
        # lowercase the param string
        args.param_string = args.param_string.lower()

        num_threads = set_thread_pool(args.cores)

        notify(f"downloading and sketching all accessions in '{args.input_csv} using {args.retry_times} retries and {num_threads} threads")

        super().main(args)
        status = sourmash_plugin_directsketch.do_gbsketch(args.input_csv,
                                                           args.param_string,
                                                           args.failed,
                                                           args.checksum_fail,
                                                           args.retry_times,
                                                           args.fastas,
                                                           args.keep_fasta,
                                                           args.genomes_only,
                                                           args.proteomes_only,
                                                           args.download_only,
                                                           args.batch_size,
                                                           args.output)
        
        if status == 0:
            notify("...gbsketch is done!")
            if args.output is not None:
                notify(f"Sigs in '{args.output}'.")
            if args.keep_fasta:
                notify(f"FASTAs in '{args.fastas}'.")

        return status


class Download_and_Sketch_Url(CommandLinePlugin):
    command = 'urlsketch'
    description = 'download and sketch GenBank assembly datasets'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('input_csv', help="a txt file or csv file containing accessions in the first column")
        p.add_argument('-o', '--output', default=None,
                       help='output zip file for the signatures')
        p.add_argument('--batch-size', type=non_negative_int, default = 0,
                       help='Write smaller zipfiles, each containing sigs associated with this number of accessions. \
                            This allows urlsketch to recover after unexpected failures, rather than needing to \
                            restart sketching from scratch. Default: write all sigs to single zipfile.')
        p.add_argument('-f', '--fastas',
                       help='Write fastas here', default = '.')
        p.add_argument('-k', '--keep-fasta', '--keep-fastq', action='store_true',
                       help="write FASTA/Q files in addition to sketching. Default: do not write FASTA files")
        p.add_argument('--download-only', help='just download genomes; do not sketch', action='store_true')
        p.add_argument('--failed',help='csv of failed accessions and download links.', required=True)
        # don't require checksum_fail here b/c users don't need to provide checksums
        p.add_argument('--checksum-fail', help="csv of accessions where the md5sum check failed. If not provided, md5sum failures will be written to the download failures file (no additional md5sum information).", default=None)
        p.add_argument('-p', '--param-string', action='append', type=str, default=[],
                          help='parameter string for sketching (default: k=31,scaled=1000)')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')
        p.add_argument('-r', '--retry-times', default=1, type=int,
                       help='number of times to retry failed downloads')


    def main(self, args):
        print_version()
        if not args.param_string:
            args.param_string = ["k=31,scaled=1000"]
        notify(f"params: {args.param_string}")

        if args.download_only and not args.keep_fasta:
            notify("Error: '--download-only' requires '--keep-fasta'.")
            sys.exit(-1)
        if args.output is None and not args.download_only:
            notify("Error: output signature zipfile is required if not using '--download-only'.")
            sys.exit(-1)

        # convert to a single string for easier rust handling
        args.param_string = "_".join(args.param_string)
        # lowercase the param string
        args.param_string = args.param_string.lower()

        num_threads = set_thread_pool(args.cores)

        notify(f"downloading and sketching all accessions in '{args.input_csv} using {num_threads} threads")

        super().main(args)
        status = sourmash_plugin_directsketch.do_urlsketch(args.input_csv,
                                                           args.param_string,
                                                           args.failed,
                                                           args.retry_times,
                                                           args.fastas,
                                                           args.keep_fasta,
                                                           args.download_only,
                                                           args.batch_size,
                                                           args.output,
                                                           args.checksum_fail)
        
        if status == 0:
            notify("...gbsketch is done!")
            if args.output is not None:
                notify(f"Sigs in '{args.output}'.")
            if args.keep_fasta:
                notify(f"FASTAs in '{args.fastas}'.")

        return status
