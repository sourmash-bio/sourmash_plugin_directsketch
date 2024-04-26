"""DirectSketch plugin description"""

usage="""
   sourmash scripts gbsketch
"""

epilog="""
See https://github.com/xyz for more examples.

Need help? Have questions? Ask at http://github.com/sourmash-bio/sourmash/issues!
"""

import argparse
import sourmash

from sourmash.logging import debug_literal, notify
from sourmash.plugins import CommandLinePlugin


class Download_and_Sketch_Assemblies(CommandLinePlugin):
    command = 'gbsketch'
    description = 'sourmash plugin to do Download and Sketch GenBank Assembly Datasets'

    def __init__(self, p):
        super().__init__(p)
        p.add_argument('accessions_csv', help="a txt file or csv file containing accessions in the first column")
        p.add_argument('-o', '--output', required=True,
                       help='output zip file for the signatures')
        p.add_argument('-f', '--fasta-output',
                       help='Write fastas here')
        p.add_argument('-k', '--keep-all-fastas',
                       help="keep all fastas after sketching. Default: just keep genomes for failed protein downloads.")
        p.add_argument('--failed', 'csv of failed accessions and download links (should be mostly protein).')
        p.add_argument('-p', '--param-string', action='append', type=str, default=[],
                          help='parameter string for sketching (default: k=31,scaled=1000)')
        p.add_argument('-c', '--cores', default=0, type=int,
                       help='number of cores to use (default is all available)')
        p.add_argument('-r', '--retry-times', default=1, type=int,
                       help='number of times to retry failed downloads)')
        # p.add_argument('-s', '--singleton', action="store_true",
        #                help='build one sketch per FASTA record, i.e. multiple sketches per FASTA file')

    def main(self, args):
        # print_version()
        if not args.param_string:
            args.param_string = ["k=31,scaled=1000"]
        notify(f"params: {args.param_string}")

        # convert to a single string for easier rust handling
        args.param_string = "_".join(args.param_string)
        # lowercase the param string
        args.param_string = args.param_string.lower()

        # num_threads = set_thread_pool(args.cores)

        # notify(f"sketching all files in '{args.fromfile_csv}' using {num_threads} threads")

        super().main(args)
        status = sourmash_plugin_directsketch.do_gbsketch(args.accessions_csv,
                                                           args.param_string,
                                                           args.output,
                                                           args.fasta_output,
                                                           args.keep_all_fastas)
        
        if status == 0:
            notify(f"...gbsketch is done! Sigs in '{args.output}'. Fastas in '{args.fasta_output}'")
        return status
