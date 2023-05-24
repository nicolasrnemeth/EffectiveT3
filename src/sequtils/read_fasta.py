import re, os, sys
import numpy as np
from typing import List

def read_fasta(file: str) -> np.ndarray:
    """
        Parses input fasta file.
        
        Args:
            file (str): input fasta-file name
            
        Returns:
            numpy-array of dimension n x m where n is the number
            of protein sequences and m being 2 for 2-sized list
            containing the protein identifier (str) and sequence (str)
    """
    if os.path.exists(file) == False:
        print('Error: "' + file + '" does not exist.')
        sys.exit(1)

    with open(file, 'r') as ifile:
        records = ifile.read()

    if re.search('>', records) == None:
        print('The input file is not in a valid fasta format.')
        sys.exit(1)

    records = records.split('>')[1:]
    myFasta = list()
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append([name, sequence])
    return np.array(myFasta)