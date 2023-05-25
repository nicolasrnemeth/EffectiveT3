import os
import sys
import json
import time
import pickle
from typing import Any
from .trainer import Trainer
from argparse import ArgumentParser, RawTextHelpFormatter

"""
    Once the training has finished add model_G1.bin, model_G2.bin and model_G3.bin
    inside the folder 'SAVED_MODELS_FOLDER' (see variable below) inside src/models.
    src/models contains the models that will be used by the prediction script.
"""

# Constants
# Folder where to save models and parameters
SAVED_MODELS_FOLDER = "src/training/models_and_parameters"

DESCRIPTION = """Trains the model and saves the optimized hyperparameters and the model 
to the folder 'src/training/models_and_parameters/'."""

# Original model was trained on this sequence region
# if you train and new model and change the sequence region
# make sure to adjust this variable at the top of the file __predict__.py
# such that it has the same value
SEQ_RANGE = (1, 26)


def parse_args():
    parser = ArgumentParser(description=DESCRIPTION,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('-p', '--pos', required=True, type=str,
                        help="(Required) Path to the fasta-file containing positive protein sequences.")

    parser.add_argument('-n', '--neg', required=True, type=str,
                        help="(Required) Path to the fasta-file containing negative protein sequences.")

    return parser.parse_args()


def convert_seconds(seconds):
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    return f"{days}d {hours}h {minutes}m {seconds}s"


def save_model(model_path: str, param_path: str, model: Any, parameters: dict) -> None:
    with open(model_path, 'wb') as ofile:
        pickle.dump(model, ofile)
    with open(param_path, 'w') as ofile:
        json.dump(parameters, ofile)


def start(pargs: dict) -> None:
    start = time.time()

    print("\nLoading models and computing encodings ...\n")
    trainer = Trainer(pargs.pos, pargs.neg, seq_range=SEQ_RANGE)
    print("Training Model ...\n")
    model, parameters = trainer.train()
    save_model(SAVED_MODELS_FOLDER+"/model.bin", SAVED_MODELS_FOLDER+"/optimized_parameters.json",
               model, parameters)
    print("Done! Model and parameters saved!")

    print(f"\nTraining took {convert_seconds(time.time() - start)}\n")


def main():
    # try:
    # Parse command line arguments
    args = parse_args()
    # Start the training program
    start(args)
    folder_path = os.path.join(
        os.getcwd(), "src/training/models_and_parameters")
    destination_path = os.path.join(os.getcwd(), "models")
    print('\nSuccessful execution of training!')
    print('\n--> Please find the saved models and optimized hyperparameters here: ' + folder_path)
    print('\n\n move the models model_G1.bin, model_G2.bin, model_G3.bin from this folder into ',
          destination_path, "if you want to use the newly trained models for prediction.")
    sys.exit(0)
    # except Exception as e:
    #    print('Program ran into an error: ', str(e))
    #    sys.exit(0)


if __name__ == '__main__':
    main()
