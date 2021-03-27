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
    parameter -r / --redo_analysis= : existing results will be recalculated
    parameter --visualize= : calls visualization for all existing results
    parameter --list_chromosomes: prints a list of all genomes with their associated chromosomes
    parameter --list_downloaded_data: prints a table containing all downloaded genomes, biosources, tfs and chromosomes
    """

    # import score, author: Noah
    import scripts.score

    # import analyse_main, author: Jan
    import scripts.analyse_main

    # import generate_data, author: Jonathan
    import generate_data

    # call deepblue to get all available genomes, biosources and tfs
    url = "http://deepblue.mpi-inf.mpg.de/xmlrpc"
    server = xc.Server(url, encoding='UTF-8', allow_none=True)
    user_key = "anonymous_key"

    try:
        genome_choices = [genome[1].lower() for genome in server.list_genomes(user_key)[1]]
        genome_choices.sort()

        biosource_choices = [biosource[1].lower().replace("'","\\'") for biosource in server.list_biosources(None, user_key)[1]]
        biosource_choices.append('all')
        biosource_choices.sort()

        tf_choices = [tf[1].lower() for tf in
                      server.list_epigenetic_marks({"category": "Transcription Factor Binding Sites"}, user_key)[1]]
        tf_choices.append('all')
        tf_choices.sort()

        chromosomes = {}
        chr_choices = []
        for genome in genome_choices:
            chrom = [x[0] for x in server.chromosomes(genome, user_key)[1]]
            chr_choices += chrom
            chromosomes[genome] = natsorted(chrom)
    except xc.ProtocolError:
        raise Exception('DeepBlue is currently offline.')

    # links to deepbluer for possible genomes, biosources and tfs
    # given to the user if script is called with -h / --help
    genomelink = 'https://deepblue.mpi-inf.mpg.de/dashboard.php#ajax/deepblue_view_genomes.php'
    biosourcelink = 'https://deepblue.mpi-inf.mpg.de/dashboard.php#ajax/deepblue_view_biosources.php'
    tflink = 'https://deepblue.mpi-inf.mpg.de/dashboard.php#ajax/deepblue_view_epigenetic_marks.php'

    # parse arguments submitted by the user
    parser = argparse.ArgumentParser(description='Analysis of transcription factors',
                                     formatter_class=RawTextHelpFormatter)

    parser.add_argument('-g', '--genome', default='hg19', type=str, choices=genome_choices,
                        help='Allowed values are genome names from\n' + genomelink,
                        metavar='GENOME')
    parser.add_argument('-b', '--biosource', default=['all'], type=str, nargs='+', choices=biosource_choices,
                        help='Allowed values are \'all\' or biosources from\n' + biosourcelink,
                        metavar='BIOSOURCE')
    parser.add_argument('-t', '--tf', default=['all'], type=str, nargs='+', choices=tf_choices,
                        help='Allowed values are \'all\' or epigenetic marks from\n' + tflink,
                        metavar='TF')
    parser.add_argument('-c', '--chromosome', default=['all'], type=str, nargs='+', choices=chr_choices,
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
    parser.add_argument('-r', '--redo_analysis', action='store_true',
                        help='existing results will be executed again and overwritten if argument is submitted')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--visualize', action='store_true',
                       help='calls visualization for all existing results')
    group.add_argument('--list_chromosomes', action='store_true',
                       help='prints a list of all genomes with their associated chromosomes')
    group.add_argument('--list_downloaded_data', action='store_true',
                       help='prints a table containing all downloaded genomes, biosources, tfs and chromosomes')

    args = parser.parse_args()

    # test if img folder for visualization exists
    if not os.path.exists(
            os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'visualization', 'src', 'assets', 'img'))):
        os.mkdir(
            os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'visualization', 'src', 'assets', 'img')))

    # remove images in visualization folder if analysis is executed again
    if args.redo_analysis:
        for f in os.listdir(os.path.abspath(
                os.path.join(os.path.dirname(__file__), '..', 'visualization', 'src', 'assets', 'img'))):
            os.remove(f)

    # print a list of all genomes and their associated chromosomes
    if args.list_chromosomes:
        for genome in genome_choices:
            print(genome + ":\n")
            l = len(chromosomes[genome]) % 5
            i = 0
            while i < len(chromosomes[genome]) - l:
                print(" ".join(x.ljust(25) for x in chromosomes[genome][i:i + 5]))
                i += 5
            if l > 0:
                print(" ".join(x.ljust(25) for x in chromosomes[genome][i:i + l]))
            print("-" * 125)

    # print a table containing the downloaded genomes, biosources, tfs and chromosomes from the linking_table
    elif args.list_downloaded_data:

        try:
            lt = pd.read_csv(
                os.path.join(args.output_path, 'data', 'linking_table.csv'),
                sep=';', usecols=['genome', 'epigenetic_mark', 'biosource', 'chromosome'])
            lt.drop(lt[lt['epigenetic_mark'] == ('dnasei' or 'dna accessibility')].index, inplace=True)
            tab = lt.groupby(["genome", "biosource", "epigenetic_mark"]).chromosome.unique().apply(
                lambda x: '\n'.join(x)).reset_index().to_dict(orient='list')
            print(tabulate(tab, headers=['genome', 'biosource', 'transcription factor', 'chromosome'],
                           tablefmt="fancy_grid"))

        except FileNotFoundError:
            raise Exception('You have not yet downloaded any data in this directory')

    # if -v / --visualize was stated, call visualisation for existing results
    elif args.visualize:
        subprocess.Popen(['python',
                          os.path.join(os.path.dirname(__file__), 'scripts', 'visualization_app_api_start.py')])
        subprocess.call(['ng', 'serve'],
                        cwd=os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'visualization')))
    # compute analysis
    else:
        # test if biosource, tf or chromosome equals 'all'
        # if yes, set the value to all possible values from deepbluer
        if 'all' in args.biosource:
            args.biosource = biosource_choices
            args.biosource.remove('all')
        else:
            args.biosource = [x.replace("'","\\'") for x in args.biosource]

        if 'all' in args.tf:
            args.tf = [x for x in tf_choices if x != "all"]

        if 'all' in args.chromosome:
            args.chromosome = chromosomes[args.genome]

        # download data from download_dict
        requested_data = generate_data.DataConfig([args.genome], args.chromosome, args.biosource, args.tf,
                                                  args.output_path, 'linking_table.csv', 'bigwig')
        requested_data.pull_data()

        # run the script score.py and store the calculated scores in the dictionary 'scores'
        scores, exist = scripts.score.findarea(args.width, args.genome.lower(), [x.lower() for x in args.biosource],
                                               [x.lower() for x in args.tf], args.chromosome, args.output_path, args.redo_analysis)

        # test if 'scores' is an empty dictionary
        # if not, generate plots with the script analyse_main.py
        # if yes and exist is True, pass (the plots are already exist in result)
        # if yes and exist is False, notify that there is no data for the submitted combination of genome, biosource and
        # transcription factor and exit the program
        if scores:
            scripts.analyse_main.TF_analyser(n_comps=args.component_size, genome=args.genome, width=args.width,
                                             path=args.output_path, chromosome=args.chromosome).mainloop(data=scores)

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
