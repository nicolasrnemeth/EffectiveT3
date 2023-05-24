import os
import re
import sys
import uuid
import json
import pickle
import numpy as np
# Parallelize feature computation and prediction
from multiprocessing import Process, Manager
from typing import Tuple, List, Callable
from sklearn import metrics

from .sequtils import read_fasta
from .encoders.encode import encode


# Number of protein sequences required for multiprocessing
# to work if activated, regardless of the input parameter to num_core
PARALLELIZATION_THRESHOLD: int = 100
DECISION_THRESHOLD: float = 0.5


def predictor(fasta_file: str, num_cores: int, ofile_path: str = "results.txt",
              seq_range: Tuple[int, int] = None, true_labels_file_name: str = None) -> None:
    """
        Computes the prediction for protein sequences and writes the results to a .txt file

        Args:
            fasta_file (str): input fasta file containing the protein sequences
            ofile_path (str): path for the output file containing the results
            seq_range (Tuple[int, int]): the sequence range to use for prediction,
                                         if None then full-length is chosen
            num_cores (int): number of cores to use for prediction, NB: the number of cores
                             only takes effect when the number of protein sequences is significantly
                             higher than the number of input protein sequences. This is due to the fact
                             that the overhead of initializing and running multiple processes only
                             justifies when the data is large enough. Is only used when number of
                             protein sequences exceeds PARALLELIZATION_THRESHOLD, otherwise this parameter is ignored

        Returns:
            None
    """
    # Read in data
    fastas = read_fasta.read_fasta(fasta_file)
    # Split data into ´num_cores´-folds and run on ´num_cores´ cores
    if len(fastas) > PARALLELIZATION_THRESHOLD and len(fastas) > num_cores and num_cores > 1:
        data_splits = list()
        size = len(fastas) // num_cores
        for i in range(num_cores-1):
            data_splits.append(fastas[i*size:(i+1)*size])
        data_splits.append(fastas[(i+1)*size:])
        # Encode protein sequences into numerical features
        predict_args = [seq_range]
        # Run processes in parallel
        probabilities = parallelize(
            predict, predict_args, data_splits, num_cores)
    # run on a single core
    else:
        probabilities = predict(fastas, seq_range)
    # Write results to file
    true_labels = None
    if true_labels_file_name is not None:
        with open(true_labels_file_name, 'r') as ifile:
            true_labels = [int(l) for l in ifile.read().split(",")]
    # Write results to file
    write_results(ofile_path, np.array(probabilities), true_labels=true_labels)


def parallelize(worker: Callable, args: List,
                data_split: List[np.ndarray], num_cores: int) -> List:
    """
        Parallelizes the execution of a task on several data splits, such
        that they can run on several cpu cores

        Args:
            worker (Callable): the task to be executed on each different cpu process
            args (List): list of arguments passed to the worker in the Process
            data_split (List[np.ndarray]): the data splits to be executed on different cpu cores
            num_cores (int): number of parallel tasks to execute

        Returns:
            List containing the concatenated results of each
            tasks in the original order they were input in
    """
    with Manager() as manager:
        # Dictionary containing the result of each process
        results_processes = manager.dict()

        process_ids = [uuid.uuid4() for _ in range(num_cores)]
        processes = list()
        for pid, split in zip(process_ids, data_split):
            p = Process(target=worker, args=tuple(
                [split] + args + [results_processes, pid]))
            processes.append(p)
            p.start()

        # Ensure processes complete before continuing execution
        for p in processes:
            p.join()

        # Join the results preserving the original order
        output = list()
        for pid in process_ids:
            output += list(results_processes.get(pid))
    # Return accumulated result from
    return output


def predict(fastas: np.ndarray, seq_range: Tuple[int, int],
            results_dict: dict = None, process_id: str = None) -> List[float]:
    """
        Encodes protein sequences, computes the prediction probability
        and determines the ensemble prediction for all chosen models

        Args: 
            fastas (np.ndarray): array containing the protein sequences
            seq_range (Tuple[int, int]): the sequence range to use for prediction
            results_dict (dict): the dictionary the results from the
                                 multiple processes are written to
            process_id (str): uuid of the process, to concatenate the results in the right order

        Returns:
            None
    """
    names, features = encode(fastas, seq_range)

    model = load_model()

    # Probability for positive label, i.e. secreted protein
    probas = model.predict_proba(features)[:, 1][:, None]

    if results_dict is not None and process_id is not None:
        results_dict[process_id] = probas

    return probas


def load_model() -> object:
    # Load trained model
    try:
        with open(os.path.join("src", "model", "model.bin"), 'rb') as ifile:
            model = pickle.load(ifile)
    except Exception as e:
        print("The EffectiveT3 model could not be loaded. Make sure it is present in the directory: ",
              + os.path.join(os.getcwd(), "src", "model"))
        print("If not it needs to be generated by running the training script or downloaded it again ",
              "from the GitHub-repository.")
        print("\nERROR MESSAGE:", str(e))
        sys.exit(0)
    return model


def write_results(ofile_path: str, probabilities: np.ndarray, model_: str, true_labels: List[int] = None) -> None:
    """
        Write prediction results to output file (either .json or .txt depending 
        on the file extension in the provided file name)
        Decision threshold for prediction is DECISION_THRESHOLD: 
            larger-equal than DECISION_THRESHOLD -> positive (secreted)
            smaller than DECISION_THRESHOLD -> negative (not secreted)

        Args:
            ofile_path (str): path of the output file containing the prediction results
            probabilities (np.ndarray): probabilities of protein being secreted
            model_ (str): the selected combination of models, by default it is 
                          the ensemble model consisting of model_G1, model_G2 and model_G3
            true_labels (List[int]): list containing the true labels of the input protein sequences

        Returns:
            None
    """
    # np.ndarray containing the boolean labels for each protein sequence
    labels = (probabilities >= DECISION_THRESHOLD).astype(bool)

    file_path, file_ext = os.path.splitext(ofile_path)

    # If true labels are present then compute all kinds of evaluation metrics
    try:
        if true_labels is not None:
            evaluation_metrics(file_path, y_pred=labels.astype(int).tolist(),
                               y_true=true_labels, y_probas=probabilities)
    except Exception as e:
        print("\n\nThere seems to be an error with your file containing the comma-separated labels")
        print("The evaluation metrics therefore could not be computed!")
        print("Please check the syntax of your file and whether there are as many labels as there are protein sequences.")
        print("ERROR MESSAGE: ", str(e), "\n\n")

    # Write prediction results to file
    if ".json" not in file_ext and ".txt" not in file_ext:
        ofile_path = file_path + ".txt"
    with open(ofile_path, 'w') as ofile:
        if ".json" in re.split(r'[/\\]', ofile_path)[-1]:
            results_dict = {str(key): dict(label=str(lab), probability=round(
                prob, 3)) for key, lab, prob in zip(range(len(labels)), labels, probabilities)}
            json.dump(results_dict, ofile, indent=4)
        else:
            dashes = "-"*30
            results = f"PREDICTION RESULTS\n{dashes}\n"
            results += str(sum(labels)) + " pos. / " + \
                str(len(labels)-sum(labels))
            results += " neg.  --> " + \
                str(round(sum(labels)/len(labels)*100, 4)) + \
                f" % positives\n{dashes}\n\n"
            results += "sequence number, prediction by " + \
                "Effective T3" + f", probability\n{dashes}\n"
            for seqNo, lab, proba in zip(range(len(labels)), labels, probabilities):
                results += '> ' + str(seqNo) + ', ' + str(lab) + \
                    ', ' + str(round(proba, 3)) + '\n'
            ofile.write(results[:-1])


def evaluation_metrics(file_path: str, y_pred: List[int], y_true: List[int], y_probas: List[float]) -> None:
    """
        Compute all kinds of metrics and save them to .json-file

        Args:
            file_path (str): path where to save metrics file
            y_pred (List[int]): predicted labels
            y_true (List[int]): true labels
            y_probas (List[float]): probabilities instead of actual labels

        Returns:
            None
    """
    # Compute metrics
    metrics_dict = {
        'accuracy': metrics.accuracy_score(y_true, y_pred),
        'precision': metrics.precision_score(y_true, y_pred),
        'recall': metrics.recall_score(y_true, y_pred),
        'f1_score': metrics.f1_score(y_true, y_pred),
        'roc_auc': metrics.roc_auc_score(y_true, y_pred),
        'log_loss': metrics.log_loss(y_true, y_pred),
        'confusion_matrix': metrics.confusion_matrix(y_true, y_pred).tolist()
    }

    # Compute ROC curve
    fpr, tpr, _ = metrics.roc_curve(y_true, y_probas)
    roc_curve_dict = {'fpr': fpr.tolist(), 'tpr': tpr.tolist()}
    metrics_dict['roc_curve'] = roc_curve_dict

    # Compute Precision-Recall curve
    precision, recall, _ = metrics.precision_recall_curve(y_true, y_probas)
    pr_curve_dict = {'precision': precision.tolist(),
                     'recall': recall.tolist()}
    metrics_dict['precision_recall_curve'] = pr_curve_dict

    # Save metrics to a JSON file
    with open(file_path + '_metrics.json', 'w') as ofile:
        json.dump(metrics_dict, ofile, indent=4)
