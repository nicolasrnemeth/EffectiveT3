# Effective T3 Version 3.0

Binary classification for predicting whether proteins are secreted by the bacterial Type III secretion system based on
a light gradient boosting machine model that was trained on sequence- and amino acid property based feature groups.

## INTRODUCTION

Type 3 secretion systems (T3SS) are critical components of many Gram-negative bacterial pathogens, playing a vital role in the pathogenesis of various infections. Accurate prediction of T3SS secreted proteins, also known as effectors, is crucial for understanding the molecular mechanisms underlying bacterial pathogenicity and developing novel therapeutic strategies. We present a new version of Effective T3, a machine learning approach for predicting T3SS effectors, surpassing the current state-of-the-art (SOTA) software Bastion3 in terms of gener-alization and computation time while achieving similar prediction performance.
Unlike Bastion3, which heavily relies on position-specific scoring matrices (PSSMs) for pre-dictions, Effective T3 employs a more flexible framework that solely incorporates N-terminal protein-sequence-based information which we should better predict evolutionary distinct T3SS effectors - more precisely those that deviate from experimentally verified ones. Furthermore, EffectiveT3 offers faster performance and utilizes evaluation metrics more suitably for imbalanced datasets, contrasting with Bastion3's heavy reliance on the ROC-AUC metric, which can be overly optimistic when we are dealing with imbalanced datasets.
To ensure robust and reliable predictions, Effective T3 was trained on a dataset, which was compiled using more rigorous criteria than the Bastion3 training dataset. Key differences include a stricter reduction of redundancy in protein sequences and the incorporation of a negative dataset solely consisting of protein sequences from bacterial species with experi-mentally verified T3SS effectors.
By developing a more rigorous training da¬¬taset and employing a novel approach that es-chews PSSMs, Effective T3 sets a new standard for T3SS effector prediction.

## Install dependencies using the following command

### Note: they will be automatically installed via pip using the install command below

> pip install -r requirements.txt

## Install using pip

install from PyPI:

> pip install EffectiveT3

install locally:

> pip install .

## Commands available by the program (check help description by running below commands)

> effectivet3 --help

> effectiveTrain --help

The 1. command is used for prediction.

The 2. command is used for training a model.

## Configure training parameters used for training script

Inside the file 'training_config.yaml' you can change the existing training parameters.
If configuration file cannot be loaded due to an error or there are missspellings for certain keys,
then the default values will be taken. In this case you will be informed via console-output.
They contain the original parameter settings as were used to obtain the model behind Effective T3 Version 3.0

## LICENSE

see LICENSE.txt

### Original papers of the previous Effective T3 versions:

- Version 2.0:

Eichinger V, Nussbaumer T, Platzer A, Jehl MA, Arnold R, Rattei T.
EffectiveDB--updates and novel features for a better annotation of bacterial secreted proteins and Type III, IV, VI secretion systems.
Nucleic Acids Res. 2016 Jan 4;44(D1):D669-74. doi: 10.1093/nar/gkv1269.
Epub 2015 Nov 20. PMID: 26590402; PMCID: PMC4702896.

- Version 1.0:

Arnold R, Brandmaier S, Kleine F, Tischler P, Heinz E, Behrens S, Niinikoski A, Mewes HW, Horn M, Rattei T.
Sequence-based prediction of type III secreted proteins. PLoS Pathog. 2009 Apr;5(4):e1000376.
doi: 10.1371/journal.ppat.1000376. Epub 2009 Apr 24. Erratum in: PLoS Pathog. 2009 Apr;5(4).
doi: 10.1371/annotation/78659a32-7869-4b14-91a6-b301a588d937. PMID: 19390696; PMCID: PMC2669295.
