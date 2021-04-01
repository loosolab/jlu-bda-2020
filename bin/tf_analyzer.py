"""
@author: Jasmin
"""
import pandas as pd
import argparse
import os
from argparse import RawTextHelpFormatter
import xmlrpc.client as xc
import subprocess
from natsort import natsorted
from tabulate import tabulate
import numpy as np
import regex
import sys


def main():
    """
    This function ties in all transcription factor analysis functions. It provides the interface for the user to pass
    parameters and run the analysis. The parameters are submitted via the command line.
    parameter  -h : help
    parameter -g / --genome= : genome
    parameter -b / --biosource= : one or more biosources divided by space
    parameter -t / --tf= : one or more transcription factors divided by space
    parameter -c / --chromosome: one or more chromosomes divided by space
    parameter -w / --width= : a parameter that determines the size of the areas to be analyzed on the chromosomes.
                              The start position is determined by subtracting the width from the summit of the peak,
                              the end position is determined by adding the width to the summit.
    parameter -o / --output_path: the path were the data and results should be stored
    parameter -cs / --component_size: single integer determining the component size for the analysis
    parameter --redo_analysis= : existing results will be recalculated
    parameter --visualize= : calls visualization for all existing results
    parameter --list_chromosomes: prints a list of all genomes with their associated chromosomes
    parameter --list_downloaded_data: prints a table containing all downloaded genomes, biosources, tfs and chromosomes
    parameter --offline: runs the program in offline mode
    parameter --redo_file_validation: downloaded files will be validated and sorted again
    parameter --check_local_files: This parameter takes the path of an external directory containing Deepblue data.
                                  The files matching query will be copied from here instead of downloaded.
    """

    # import score, author: Noah
    import scripts.score

    # import analyse_main, author: Jan
    import scripts.analyse_main

    # import generate_data, author: Jonathan
    import generate_data

    # links to deepbluer for possible genomes, biosources and tfs
    # given to the user if script is called with -h / --help
    genomelink = 'https://deepblue.mpi-inf.mpg.de/dashboard.php#ajax/deepblue_view_genomes.php'
    biosourcelink = 'https://deepblue.mpi-inf.mpg.de/dashboard.php#ajax/deepblue_view_biosources.php'
    tflink = 'https://deepblue.mpi-inf.mpg.de/dashboard.php#ajax/deepblue_view_epigenetic_marks.php'

    # parse arguments submitted by the user
    parser = argparse.ArgumentParser(description='Analysis of transcription factors',
                                     formatter_class=RawTextHelpFormatter)

    parser.add_argument('-g', '--genome', default='hg19', type=str,
                        help='Allowed values are genome names from\n' + genomelink,
                        metavar='GENOME')
    parser.add_argument('-b', '--biosource', default=['all'], type=str, nargs='+',
                        help='Allowed values are \'all\' or biosources from\n' + biosourcelink,
                        metavar='BIOSOURCE')
    parser.add_argument('-t', '--tf', default=['all'], type=str, nargs='+',
                        help='Allowed values are \'all\' or epigenetic marks from\n' + tflink,
                        metavar='TF')
    parser.add_argument('-c', '--chromosome', default=['all'], type=str, nargs='+',
                        help='Allowed values are \'all\' or or the chromosomes belonging to the selected genome. \n'
                             'To display the possible chromosomes, call the program with the parameter '
                             '--list_chromosomes.',
                        metavar='CHROMOSOME')
    parser.add_argument('-w', '--width', default=50, type=int,
                        help='A parameter that determines the size of the areas to be analyzed on the chromosomes. \n'
                             'The start position is determined by subtracting the width from the summit of the peak, \n'
                             'the end position is determined by adding the width to the summit.')
    parser.add_argument('-o', '--output_path', default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
                        type=str, nargs='?', help='The path were the downloaded data and the results will be stored.')
    parser.add_argument('-cs', '--component_size', type=int, nargs='?',
                        help='single integer determining the component size for the analysis')
    parser.add_argument('--redo_analysis', action='store_true',
                        help='existing results will be executed again and overwritten if argument is submitted')
    parser.add_argument('--offline', action='store_true',
                        help='runs the program in offline mode')
    parser.add_argument('--redo_file_validation', action='store_true',
                        help='downloaded files will be validated and sorted again')
    parser.add_argument('--check_local_files', type=str, nargs='?',
                        help='This parameter takes the path of an external directory containing Deepblue data. The '
                             'files matching query will be copied from here instead of downloaded.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--visualize', action='store_true',
                       help='calls visualization for all existing results')
    group.add_argument('--list_chromosomes', action='store_true',
                       help='prints a list of all genomes with their associated chromosomes')
    group.add_argument('--list_downloaded_data', action='store_true',
                       help='prints a table containing all downloaded genomes, biosources, tfs and chromosomes')

    args = parser.parse_args()

    try:
        lt = pd.read_csv(
            os.path.join(args.output_path, 'data', 'linking_table.csv'),
            sep=';', usecols=['genome', 'epigenetic_mark', 'biosource', 'chromosome'])
        lt.drop(lt[lt['epigenetic_mark'] == ('dnasei' or 'dna accessibility')].index, inplace=True)
    except FileNotFoundError:
        lt = None

    # test if img folder for visualization exists
    if not os.path.exists(
            os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'visualization', 'src', 'assets', 'img'))):
        os.mkdir(
            os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'visualization', 'src', 'assets', 'img')))

    # print a table containing the downloaded genomes, biosources, tfs and chromosomes from the linking_table
    if args.list_downloaded_data:
        try:
            tab = lt.groupby(["genome", "biosource", "epigenetic_mark"]).chromosome.unique().apply(
                lambda x: '\n'.join(x)).reset_index().to_dict(orient='list')
            print(tabulate(tab, headers=['genome', 'biosource', 'transcription factor', 'chromosome'],
                           tablefmt="fancy_grid"))
        except AttributeError:
            sys.exit('There is currently no downloaded data in your specified path.')

    # if -v / --visualize was stated, call visualisation for existing results
    elif args.visualize:
        subprocess.Popen(['python',
                          os.path.join(os.path.dirname(__file__), 'scripts', 'visualization_app_api_start.py')])
        subprocess.call(['ng', 'serve'],
                        cwd=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'visualization')))

    else:
        if args.offline:
            if lt is not None:
                genome_choices = list(set(lt['genome']))
                biosource_choices = list(set(lt['biosource']))
                tf_choices = list(set(lt['epigenetic_mark']))
                chromosomes = lt.groupby(["genome"]).chromosome.unique().to_dict()
            else:
                sys.exit('There is currently no downloaded data in your specified path.')
        else:
            # call deepblue to get all available genomes, biosources and tfs
            url = "http://deepblue.mpi-inf.mpg.de/xmlrpc"
            server = xc.Server(url, encoding='UTF-8', allow_none=True)
            user_key = "anonymous_key"

            try:
                genome_choices = [genome[1].lower() for genome in server.list_genomes(user_key)[1]]
                genome_choices.sort()

                biosource_choices = [biosource[1].lower().replace("'", "\\'") for biosource in
                                     server.list_biosources(None, user_key)[1]]
                biosource_choices.append('all')
                biosource_choices.sort()

                tf_choices = [tf[1].lower() for tf in
                              server.list_epigenetic_marks({"category": "Transcription Factor Binding Sites"},
                                                           user_key)[1]]
                tf_choices.append('all')
                tf_choices.sort()

                chromosomes = {}
                for genome in genome_choices:
                    chrom = [x[0].lower() for x in server.chromosomes(genome, user_key)[1]]
                    chromosomes[genome] = natsorted(chrom)
            except Exception as e:
                sys.exit(e + ': No connection could be established to deepblue. If you want to run an analysis with '
                             'already downloaded data, please call the program with the parameter \'--offline\'.')

        # print a list of all genomes and their associated chromosomes
        if args.list_chromosomes:
            if not args.offline:
                for genome in genome_choices:
                    print(genome + ':\n')
                    print('\n'.join([' '.join(x.ljust(len(max(chromosomes[args.genome], key=len))) for x in g) for g in
                                     np.split(chromosomes[genome], range(5, len(chromosomes[genome]), 5))]))
                    print('-' * 125)
            else:
                sys.exit('You need a connection to Deepblue to list all possible chromosomes.')

        # compute analysis
        else:
            # test if submitted chromosomes, biosources and tfs are valid
            if args.genome not in genome_choices:
                parser.error('argument -g/--genome: invalid choice: \'' + args.genome + '\', choose from:\n' +
                             '\t'.join(x.ljust(len(max(genome_choices, key=len))) for x in genome_choices))
            for chrom in args.chromosome:
                if chrom != 'all' and chrom not in chromosomes[args.genome]:
                    parser.error(
                        'argument -c/--chromosome: invalid choice: \'' + chrom + '\' , choose from:\n' + '\n'.join(
                            ['\t'.join(x.ljust(len(max(chromosomes[args.genome], key=len))) for x in c) for c in
                             np.split(chromosomes[args.genome], range(5, len(chromosomes[args.genome]), 5))]))
            for bs in args.biosource:
                if bs != 'all' and bs not in biosource_choices:
                    possible_biosources = list(filter(regex.compile('.*' + bs + '.*').match, biosource_choices))[0:30]
                    if len(possible_biosources) == 0:
                        possible_biosources = biosource_choices[0:30]
                    parser.error('argument -b/--biosource: invalid choice: \'' + bs + '\', possible arguments are:\n' +
                                 '\n'.join(
                                    ['\t'.join(x.ljust(len(max(possible_biosources, key=len))) for x in b) for b in
                                        np.split(possible_biosources, range(2, len(possible_biosources), 2))]))
            for one_tf in args.tf:
                if one_tf != 'all' and one_tf not in tf_choices:
                    possible_tfs = list(filter(regex.compile('.*' + one_tf + '.*').match, tf_choices))[0:30]
                    if len(possible_tfs) == 0:
                        possible_tfs = tf_choices[0:30]
                    parser.error(
                        'argument -t/--transcription_factor: invalid choice: \'' + one_tf + '\', possible '
                                                                                            'arguments are:\n' +
                        '\n'.join(['\t'.join(x.ljust(len(max(possible_tfs, key=len))) for x in t) for t in
                                   np.split(possible_tfs, range(4, len(possible_tfs), 4))]))

            # test if biosource, tf or chromosome equals 'all'
            # if yes, set the value to all possible values from deepbluer
            if 'all' in args.biosource:
                args.biosource = biosource_choices
            else:
                args.biosource = [x.replace("'", "\\'") for x in args.biosource]

            if 'all' in args.tf:
                args.tf = tf_choices

            if 'all' in args.chromosome:
                args.chromosome = chromosomes[args.genome]

            # download data from download_dict
            requested_data = generate_data.DataConfig([args.genome], args.chromosome, args.biosource, args.tf,
                                                      args.output_path, 'linking_table.csv', 'bigwig',
                                                      args.check_local_files, args.redo_file_validation, args.offline)
            requested_data.pull_data()

            # run the script score.py and store the calculated scores in the dictionary 'scores'
            scores, exist = scripts.score.findarea(args.width, args.genome.lower(), [x.lower() for x in args.biosource],
                                                   [x.lower() for x in args.tf], args.chromosome, args.output_path,
                                                   args.redo_analysis, args.component_size)

            # test if 'scores' is an empty dictionary
            # if not, generate plots with the script analyse_main.py
            # if yes and exist is True, pass (the plots are already exist in result)
            # if yes and exist is False, notify that there is no data for the submitted combination of genome, biosource
            # and transcription factor and exit the program
            if scores:
                scripts.analyse_main.TF_analyser(n_comps=args.component_size, genome=args.genome, width=args.width,
                                                 path=args.output_path, chromosome=args.chromosome).mainloop(
                    data=scores)

            if scores or exist:

                subprocess.Popen(['python',
                                  os.path.join(os.path.dirname(__file__), 'scripts', 'visualization_app_api_start.py')])
                subprocess.call(['ng', 'serve'],
                                cwd=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'visualization')))

            else:
                raise Exception(
                    'there is no data for your entered combination of genome, biosource and transcription factor')


if __name__ == '__main__':
    main()
