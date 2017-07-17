from __future__ import print_function
import matplotlib as mpl

mpl.use('pdf')
import argparse
import sys

from bracerlib.tasks import Assembler, Summariser, Tester, Builder


def launch():
    parser = argparse.ArgumentParser(
        description='BraCeR: reconstruction of BCR sequences from single-cell RNAseq data',
        usage=''' bracer <mode> [<args>]

              Modes are :

              - assemble: assemble BCR sequences from single-cell RNA-sequencing reads
              - summarise: summarise BCR sequences from set of cells, build clonotype networks
              - test : use a small dataset from three cells to test BraCeR installation
              - build : build resource files from gene segment sequences

              use bracer <mode> -h for specific help
              ''')
    parser.add_argument('mode', metavar="<MODE>", help='bracer mode (assemble, summarise, test or build)',
                        choices=['assemble', 'summarise', 'summarize', 'test', 'build'])
    args = parser.parse_args(sys.argv[1:2])

    task_mapper = {
        'assemble': Assembler,
        'summarise': Summariser,
        'summarize': Summariser,
        'test': Tester,
        'build': Builder
    }

    if args.mode not in task_mapper:
        print('Unrecognised mode')
        parser.print_help()
        exit(1)

    Task = task_mapper[args.mode]
    Task().run()
