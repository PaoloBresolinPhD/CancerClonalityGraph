import numpy as np
import random as rd

def get_number_trees(dataset):
    """
    Counts the number of trees in the input dataset.

    Parameters:
    - dataset: numpy array with patients. Each patient contains a set of plausible trees.

    Returns:
    - num_trees: number of trees in the input dataset.
    """

    # number of trees
    num_trees = 0

    # sum the number of trees over patients
    for patient in dataset:
        num_trees += len(patient)

    return num_trees

def get_mutations(tree):
    """
    Returns the set of mutations that label nodes in the input dataset of trees.

    Parameters:
    - tree: list of edges. Each edge is a list with the two mutations it connects.

    Returns:
    - mutation_set: set with all mutation names that label nodes of the input trees.
    """

    # list with mutation names
    mutation_set = set()

    # iterate through all mutations in the input dataset
    for edge in tree:
        for mut in edge:
            if mut not in mutation_set:
                mutation_set.add(mut)
    
    return mutation_set

def collapse_gene_level(dataset):
    """
    Collapses mutations at gene level for all trees in the input dataset.
    Mutations are collapsed at gene level by truncating at the first occurrence of an underscore ('_') 
    or a period ('.').

    Parameters:
        - dataset (list): ndarray of patients, where each patient is represented as a list of trees
                          and each tree is a list of edges, which is a list of mutation pairs.

    Returns:
        - list: new version of the input dataset with mutations collapsed at gene level.
    """

    # new version of the dataset that will be returned
    new_dataset = []

    # iterate through patients
    for patient in dataset:
        
        # new version of the current patient
        new_patient = []

        # iterate through trees of the current patient
        for tree in patient:

            # new version of the current tree of the current patient
            new_tree = []

            # iterate through edges of the current tree
            for edge in tree:

                # new version of the edge
                new_edge = []

                # iterate through the pair of mutations contained in the edge
                for mutation in edge:

                    # new version of the current mutation
                    new_mutation = ""

                    # collapse the mutation at gene level
                    for i in range(len(mutation)):
                        if mutation[i] == "_" or mutation[i] == ".":
                            break
                        else:
                            new_mutation += mutation[i]
                    
                    # add the collapsed mutation to the new edge
                    new_edge.append(new_mutation)
                
                # add new_edge to new_tree
                new_tree.append(new_edge)
            
            # add new_tree to new_patient
            new_patient.append(new_tree)
        
        # add new_patient to new_dataset
        new_dataset.append(new_patient)
    
    return np.array(new_dataset, dtype='object')
