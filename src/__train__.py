import os
import sys
import json
import time
import pickle
from typing import Any

import numpy as np
from .trainer import Trainer
from argparse import ArgumentParser, RawTextHelpFormatter

"""
    Once the training has finished add model.bin,
    inside the folder 'SAVED_MODELS_FOLDER' (see variable below) inside src/models.
    src/model contains the model that will be used by the prediction script.
"""

# Constants
# Folder where to save models and parameters
SAVED_MODELS_FOLDER = os.path.join("src", "training", "model_and_parameters")

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

    parser.add_argument('-i', '--featureimportance', action="store_true",
                        help="Set this flag to save feature importances and their "
                        + "labels to a file with the same filepath as the metrics file, but named 'feature_importances.json'")

    return parser.parse_args()


def convert_seconds(seconds):
    """
        Convert seconds to a human readable string
    """
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    return f"{days}d {hours}h {minutes}m {seconds}s"


def save_feature_importance(feature_importance: np.ndarray, feat_imp_path: str):
    """
        Save feature importances to a json file
    """
    with open(os.path.join("src", "training", "feature_names.json"), 'r') as ifile:
        feature_labels = np.array(json.load(ifile))
    # Sort from highest to lowest
    sorted_indices = np.argsort(feature_importance)[::-1]
    label_importance_dict = {str(lab): float(imp) for lab, imp in zip(
        feature_labels[sorted_indices], feature_importance[sorted_indices])}
    # Save
    with open(feat_imp_path, 'w') as ofile:
        json.dump(label_importance_dict, ofile, indent=4)


def save_model(model_path: str, param_path: str, feat_imp_path: str, model: Any,
               parameters: dict, feature_importance: np.ndarray, save_feat_imp: bool = False) -> None:
    """
        Save the model, the parameters and the feature importances (if specified)
    """
    if save_feat_imp:
        # Save feature importances
        save_feature_importance(feature_importance, feat_imp_path)
    # Save model
    with open(model_path, 'wb') as ofile:
        pickle.dump(model, ofile)
    # Convert to bool because np.bool is not json serializable
    if parameters.get("extra_trees") is not None:
        parameters["extra_trees"] = bool(parameters["extra_trees"])
    # Save parameters
    with open(param_path, 'w') as ofile:
        json.dump(parameters, ofile, indent=4)


def start(pargs: dict) -> None:
    """
        Start the training program
    """
    start = time.time()

    print("\nLoading models and computing encodings ...\n")
    trainer = Trainer(pargs.pos, pargs.neg, seq_range=SEQ_RANGE)
    print("Training Model ...\n")
    model, parameters = trainer.train()
    save_model(os.path.join(SAVED_MODELS_FOLDER, "model.bin"),
               os.path.join(SAVED_MODELS_FOLDER, "optimized_parameters.json"),
               os.path.join(SAVED_MODELS_FOLDER, "feature_importances.json"),
               model, parameters, model.feature_importances_,
               # Whether to save feature importances
               save_feat_imp=pargs.featureimportance)
    print("Done! Model and parameters saved!")

    print(f"\nTraining took {convert_seconds(time.time() - start)}\n")


def main():
    try:
        # Parse command line arguments
        args = parse_args()
        # Start the training program
        start(args)
        folder_path = os.path.join(
            os.getcwd(), "src/training/models_and_parameters")
        destination_path = os.path.join(os.getcwd(), "models")
        print('\nSuccessful execution of training!')
        print('\n--> Please find the saved models and optimized hyperparameters here: ' + folder_path)
        print('\n\n move the model model.bin, from this folder into ',
              destination_path, "if you want to use the newly trained models for prediction.")
        sys.exit(0)
    except Exception as e:
        print('Program ran into an error: ', str(e))
        sys.exit(0)


if __name__ == '__main__':
    main()
