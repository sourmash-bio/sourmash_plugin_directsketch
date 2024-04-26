#! /usr/bin/env python
import os
from sourmash.logging import notify
from sourmash.plugins import CommandLinePlugin
import importlib.metadata

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
    avail_threads = get_max_cores()
    num_threads = min(avail_threads, user_cores) if user_cores else avail_threads
    if user_cores and user_cores > avail_threads:
        notify(f"warning: only {avail_threads} threads available, using {avail_threads}")
    actual_rayon_cores = sourmash_plugin_directsketch.set_global_thread_pool(num_threads)
    return actual_rayon_cores

class Download_and_Sketch_Assemblies(CommandLinePlugin):
    command = 'gbsketch'
    description = 'download and sketch GenBank assembly datasets'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('accessions_csv', help="a txt file or csv file containing accessions in the first column")
        p.add_argument('-o', '--output', required=True,
                       help='output zip file for the signatures')
        p.add_argument('-f', '--fasta-output',
                       help='Write fastas here')
        p.add_argument('-k', '--keep-all-fastas',
                       help="keep all fastas after sketching. Default: just keep genomes for failed protein downloads.")
        p.add_argument('--failed',help='csv of failed accessions and download links (should be mostly protein).')
        p.add_argument('-p', '--param-string', action='append', type=str, default=[],
                          help='parameter string for sketching (default: k=31,scaled=1000)')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')
        p.add_argument('-r', '--retry-times', default=1, type=int,
                       help='number of times to retry failed downloads)')
        # p.add_argument('-s', '--singleton', action="store_true",
        #                help='build one sketch per FASTA record, i.e. multiple sketches per FASTA file')

    def main(self, args):
        print_version()
        if not args.param_string:
            args.param_string = ["k=31,scaled=1000"]
        notify(f"params: {args.param_string}")

        # convert to a single string for easier rust handling
        args.param_string = "_".join(args.param_string)
        # lowercase the param string
        args.param_string = args.param_string.lower()

        num_threads = set_thread_pool(args.cores)

        notify(f"downloading and sketching all accessions in '{args.input_csv}' using {num_threads} threads")

        super().main(args)
        status = sourmash_plugin_directsketch.do_gbsketch(args.accessions_csv,
                                                           args.param_string,
                                                           args.failed,
                                                           args.output,
                                                           args.retry_times,
                                                           args.fasta_output,
                                                           args.keep_all_fastas)
        
        if status == 0:
            notify(f"...gbsketch is done! Sigs in '{args.output}'. Fastas in '{args.fasta_output}'")
        return status
