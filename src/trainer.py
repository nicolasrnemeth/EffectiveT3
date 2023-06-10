# Standard packages
import os
import json
from typing import Tuple

# External packages / libraries
import yaml  # used to load configurations in 'training_config.yaml' for training model
import numpy as np
import lightgbm as lgbm  # classifier
from sklearn_genetic import GASearchCV  # heuristic optimization
from sklearn_genetic.space import Categorical, Integer, Continuous
# deterministic optimization
from sklearn.model_selection import GridSearchCV, StratifiedKFold

# Own package imports
from .sequtils import read_fasta
from .encoders.encode import encode


class Trainer(object):
    """
        Used to train a model (e.g. Light-Gradient-Boosting Machine) to predict 
        the secretion of proteins by the Type III secretion system.

        Uses the training parameters configured in 'training_config.yaml'
    """
    # The sequence region to use for prediction
    seq_range: Tuple[int, int] = None
    # Neg. : pos. class ratio
    scale_pos_weight: float = None
    # A list of protein sequences and their identifiers
    protein_sequences: np.ndarray = None
    # The true labels
    labels: np.ndarray = None
    # encoded features of the protein sequences
    features: np.ndarray = None
    # Training parameters / configuration
    TRAINING_CONFIG: dict = dict()
    try:
        with open("training_config.yaml", 'r') as configfile:
            TRAINING_CONFIG = yaml.safe_load(configfile)
    except Exception as e:
        print(f"\nTraining configuration file {os.path.join(os.getcwd(),'training_config.yaml')}"
              + " could not be found or loaded properly.")
        print("Therefore default training parameters are used instead.\n")
        print("ERROR MESSAGE: ", str(e))

    def __init__(self, pos_fasta_file: str, neg_fasta_file: str,
                 seq_range: Tuple[int, int] = None) -> None:
        """
            Creates new instance.
        """
        positive_sequences = read_fasta.read_fasta(pos_fasta_file)
        negative_sequences = read_fasta.read_fasta(neg_fasta_file)

        self.protein_sequences = np.vstack((
            positive_sequences,
            negative_sequences
        ))
        self.labels = np.hstack((
            np.ones(len(positive_sequences)),
            np.zeros(len(negative_sequences))
        ))
        # Sequence region to use for prediction
        self.seq_range = seq_range
        # The hyperparameter space to optimize over
        with open('src/training/hyperparameter_space.json', 'r') as ifile:
            self.hyperparameter_space = json.load(ifile)
        # Compute protein encodings
        self.features = encode(self.protein_sequences, seq_range)[1]
        # Weight the positive class based on the actual neg. : pos. class ratio
        y = self.labels
        # neg count divided by pos count
        self.scale_pos_weight = (len(y)-sum(y)) / sum(y)
        # Uncomment below to display feature dimensions during training
        print("Number of samples | feature dimensions:",
              self.features.shape, "\n")

    def train(self) -> Tuple[lgbm.LGBMClassifier, dict]:
        """
            Train the classifier.

            Returns:
                Tuple[lgbm.LGBMClassifier, dict]: the optimized classifier and the optimized hyperparameters
        """
        print("\nOne-by-one parameter optimization ...")
        optimized_parameters = self.__one_by_one_parameter_optimization(
            **self.__collect_params_from_config_file(step_two=False)
        )
        print("\nHeuristic optimization of one-by-one optimized parameters ...")
        final_classifier, optimized_hyperparameters = self.__heuristic_optimization(
            initial_params=optimized_parameters,
            **self.__collect_params_from_config_file(step_two=True)
        )
        return final_classifier, optimized_hyperparameters

    def __collect_params_from_config_file(self, step_two: bool) -> dict:
        """
            Collect the parameters from the training_config.yaml file.
            Parameters which are not set or missspelled will be assigned a default value.

            Args:
                step_one: whether to return parameters for first or second optimization step

            Returns:
                dict: containing training parameter values
        """
        parameters = dict()
        missingParams = list()

        if self.TRAINING_CONFIG["params"].get('n_estimators') is not None:
            parameters['n_estimators'] = self.TRAINING_CONFIG["params"]["n_estimators"]
        else:
            parameters['n_estimators'] = 30
            missingParams.append('n_estimators')

        if self.TRAINING_CONFIG["params"].get('k_fold') is not None:
            parameters['k_fold'] = self.TRAINING_CONFIG["params"]["k_fold"]
        else:
            parameters['k_fold'] = 4
            missingParams.append('k_fold')

        if self.TRAINING_CONFIG["params"].get('evaluation_metric') is not None:
            parameters['evaluation_metric'] = self.TRAINING_CONFIG["params"]["evaluation_metric"]
        else:
            parameters['evaluation_metric'] = 'balanced_accuracy'
            missingParams.append('evaluation_metric')

        if step_two:
            # Add parameters for heuristic optimization on top of common parameters
            if self.TRAINING_CONFIG["params"].get("population_size") is not None:
                parameters['population_size'] = self.TRAINING_CONFIG["params"]["population_size"]
            else:
                parameters['population_size'] = 30
                missingParams.append('population_size')

            if self.TRAINING_CONFIG["params"].get("generations") is not None:
                parameters['generations'] = self.TRAINING_CONFIG["params"]["generations"]
            else:
                parameters['generations'] = 7
                missingParams.append('generations')

        # If any parameters are missing inform the user
        if len(missingParams) != 0:
            print("The following parameters from the config filed could not be found or are missspelled: ",
                  ", ".join(missingParams))
            print("There the default values are taken instead.")

        return parameters

    def __heuristic_optimization(self, initial_params: dict, population_size: int = 30,
                                 generations: int = 7, n_estimators: int = 30, k_fold: int = 4,
                                 evaluation_metric: str = 'roc_auc') -> Tuple[lgbm.LGBMClassifier, dict]:
        """
            Performs heuristic parameter optimization using genetic algorithm.

            Args:
                initial_params (dict): parameters to optimize over
                population_size (int): how many species to generate in each generation
                generations (int): how many generations to optimize over
                n_estimators (int): number of estimators to use in boosting machine
                k_fold (int): fold of the cross-validation

            Returns:
                lgbm.LGBMClassifier: the optimized/trained classifier
        """
        param_grid = self.__parameter_grid(initial_params)
        model = lgbm.LGBMClassifier(n_jobs=-1,
                                    n_estimators=n_estimators,
                                    verbose=-1,
                                    scale_pos_weight=self.scale_pos_weight).set_params(**initial_params)

        sk_fold = StratifiedKFold(n_splits=k_fold, shuffle=True)
        clf_GA = GASearchCV(model, cv=sk_fold, param_grid=param_grid, scoring=evaluation_metric,
                            n_jobs=-1, population_size=population_size,
                            generations=generations,
                            verbose=2, error_score='raise')
        clf_GA.fit(self.features, self.labels)

        # Optimize again but over all protein sequences, because the best estimator
        # was only optimized over 'k_fold'-1 folds but not all folds due to cross-validation step
        final_classifier = clf_GA.best_estimator_
        final_classifier.fit(self.features, self.labels)

        return final_classifier, clf_GA.best_params_

    def __one_by_one_parameter_optimization(self, n_estimators: int = 30, k_fold: int = 4,
                                            evaluation_metric: str = 'roc_auc') -> dict:
        """
            Performs one by one parameter optimization.

            Args:
                n_estimators (int): number of estimators to use for gradient boosting machine
                k_fold (int): k-fold cross validation

            Returns:
                dict: the optimized parameters
        """
        best_params = dict()
        for key in self.hyperparameter_space:
            if "boosting_type" in best_params and best_params["boosting_type"] == "goss":
                if key in ["subsample", "subsample_freq"]:
                    continue
            print("Optimization with respect to " + key.upper() + ':', end=" ")
            model = lgbm.LGBMClassifier(n_jobs=-1,
                                        n_estimators=n_estimators,
                                        verbose=-1,
                                        scale_pos_weight=self.scale_pos_weight).set_params(**best_params)
            sk_fold = StratifiedKFold(n_splits=k_fold, shuffle=True)
            clf = GridSearchCV(model, {key: self.hyperparameter_space[key]}, n_jobs=-1,
                               cv=sk_fold, scoring=evaluation_metric, error_score='raise')
            clf.fit(self.features, self.labels)
            best_params.update(clf.best_params_)
            print("Done!")
        return best_params

    def __parameter_grid(self, params: dict) -> dict:
        """
            Turns the optimized parameters from ´__one_by_one_parameter_optimization´
            into the suitable format to use them for the heuristic optimization
            with the genetic algorithm.

            Args:
                params (dict): the parameters to bring into the right format

            Returns
                dict: the dictionary with the parameters in the right format
        """
        grid = dict()
        for key in params:
            if key == "boosting_type" and self.hyperparameter_space[key] == "goss":
                if key in ["subsample", "subsample_freq"]:
                    continue
            if key == 'boosting_type':
                grid.update({'boosting_type': Categorical(
                    self.hyperparameter_space[key])})
                continue
            if key == 'extra_trees':
                grid.update({'extra_trees': Categorical(
                    self.hyperparameter_space[key])})
                continue
            val = params[key]
            if key == "min_child_weight" and val == 0:
                val = 0.001
            lb = val * \
                0.2 if self.param_type(key) == 'Continuous' else int(val*0.2)
            ub = val * \
                1.8 if self.param_type(key) == 'Continuous' else int(val*1.8)+1
            if key == 'colsample_bytree' and ub > 1:
                ub = 1
            if key == "drop_rate" and ub > 1:
                ub = 1
            if key == "max_drop" and lb <= 0:
                lb = 1
            if self.param_type(key) == 'Continuous':
                grid.update({key: Continuous(lb, ub)})
            else:
                grid.update({key: Integer(lb, ub)})
        return grid

    @staticmethod
    def param_type(key: str) -> str:
        """
            Hardcoded method that returns the type of a hyperparameter.

            Args:
                key (str): name of the hyperparameter to optimize

            Returns:
                str: type of the variable
        """
        if key in ['num_leaves', 'max_depth', 'min_child_samples', 'max_bin', 'max_drop', 'subsample_freq']:
            return 'Integer'
        if key == 'boosting_type' or key == 'extra_trees':
            return 'Categorical'
        return 'Continuous'
