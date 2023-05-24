# Effective T3 Version 3.0

Binary classification for predicting whether proteins are secreted by the bacterial Type III secretion system based on
a light gradient boosting machine model that was trained on sequence- and amino acid property based feature groups.

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
