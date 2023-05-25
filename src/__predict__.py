import os
import sys
import time
from .predictor import predictor
from multiprocessing import cpu_count

from argparse import ArgumentParser, RawTextHelpFormatter
from .__init__ import __version__

# CONSTANTS

CPU_COUNT = cpu_count()

# Description strings for command line parameters
DESCRIPTION = """Classifies whether bacterial proteins are secreted by the Type III secretion system, 
based on information contained in the protein sequence. Sequence- and amino acid property-based features"""

CPU_CORES_HELP = """(Optional) The number of CPU-cores that you want to use for prediction. By default all available CPU cores are used.
If your selected number of cores is above the available number of cores you will be provided with the
number of accessible CPU-cores on your operating system, to provide a valid choice the this parameter."""


TRUE_LABELS_HELP = """"(Optional) Path to the file containing the comma-separated true labels encoded as integers (0 for False (not-secreted) and 1 for True (secreted)).
The labels must be in the same order as the protein sequences contained in the input fasta-file. If this command line parameter is set, i.e. not None
then an additional file ("{output_file_name_as_set_in_the_command_line_argument}_metrics.json") containing all kinds of evaluation metrics for the prediction will be saved.
So the json-file containing the evaluation metrics will be saved inside the same path as the outputfile name containing only the predictions but with "_metrics.json"
concatenated to the filename.

E.g. if you have 10 protein sequences the first then your labels file should contain the following:

1,1,0,1,0,0,0,1,0,0

assuming you have 4 secreted and 6 non-secreted proteins in the exact same order.

There must be no spaces in between the commas and trailing commas

"""

# Original model was trained on this sequence region
# so unless you have not trained a new model on a new sequence region
# do not change this variable
SEQ_RANGE = (1, 26)


def parse_args():
    parser = ArgumentParser(description=DESCRIPTION,
                            formatter_class=RawTextHelpFormatter)

    # Fasta file containing sequences to predict
    parser.add_argument('-f', '--file', required=True, type=str,
                        help="(Required) Path to the input fasta-file, e.g. 'your_folder/your_file.fasta'")

    # Output file path
    parser.add_argument('-o', '--ofile', required=False, type=str, default='results.txt',
                        help='(Required) Provide the file path for the output file. --ofile path/{file_name}.json to '
                        + 'save it in json-format or path/{file_name}.txt to save it in txt-format')

    # Number of cores to use for prediction
    parser.add_argument('-c', '--cores', choices=list(range(1, CPU_COUNT+1)), required=False, type=int, default=CPU_COUNT,
                        help=CPU_CORES_HELP)

    # True labels
    parser.add_argument('-l', '--truelabels', required=False, type=str, default=None,
                        help=TRUE_LABELS_HELP)

    # # Sequence range
    # parser.add_argument('-r', '--range', required=False, type=str, default=SEQ_RANGE,
    #                     help="(Optional) The range of amino acid sequences to use for prediction. "
    #                          "Default is the first 25 amino acids of the protein sequence. "
    #                          "Provide a tuple of integers, e.g. (1, 26) to use the first 25 amino acids (if first aa, i.e. methionine is excluded), "
    #                          "or (10, 26) to use the 10th to the 26th amino acid. "
    #                          "If you want to use the full sequence, provide (1, 0).")

    # Program version
    parser.add_argument('-v', '--version', action='version', version='bastion3clone ' + __version__,
                        help="(Optional) Show program's version number and exit")

    return parser.parse_args()


def convert_seconds(seconds):
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    return f"{days}d {hours}h {minutes}m {seconds}s"


def start(pargs):
    start = time.time()
    predictor(fasta_file=pargs.file, num_cores=pargs.cores, ofile_path=pargs.ofile,
              seq_range=SEQ_RANGE, true_labels_file_name=pargs.truelabels)
    print(f"\nPrediction took {convert_seconds(time.time() - start)}\n")


def main():
    # try:
    # Start program
    args = parse_args()
    start(args)
    file_path = os.path.join(os.getcwd(), args.ofile)
    print('Successful execution of the program!')
    print('\n--> Please find the results here: ' + file_path)
    sys.exit(0)
    # except Exception as e:
    #    print('Program ran into an error: ', str(e))
    #    sys.exit(0)


if __name__ == '__main__':
    main()
