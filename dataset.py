import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class GraphDataset:
    """
    Class to handle a dataset of patients with graphs.
    """

    def __init__(self, dataset, collapse_gene_level=False):
        """
        Creates an instance of Dataset class from an input numpy array.

        Parameters:
        - dataset: numpy array with patients. Each patient contains an array of graphs. Each graph stores edges into a list. An edge is a list with the two mutations it links.
        - collapse_gene_level: boolean that indicates if the mutations should be collapsed at gene level before counting the frequency of each edge.
        """
        
        # list of patients
        self.__patients = dataset

        # collapse mutations at gene level, if required
        if collapse_gene_level:
            self._collapse_gene_level()
        
        # graph representing the whole dataset
        self.__dataset_graph = self._compute_dataset_graph()

    def _collapse_gene_level(self):
        """
        Collapses mutations at gene level for all graphs in the input dataset.
        Mutations are collapsed at gene level by truncating at the first occurrence of an underscore ('_') 
        or a period ('.').
        Notice that this method modifies the dataset irreversibly.
        """

        # new version of the dataset that will be returned
        new_dataset = []

        # iterate through patients
        for patient in self.__patients:
            
            # new version of the current patient
            new_patient = []

            # iterate through graphs of the current patient
            for graph in patient:

                # new version of the current graph of the current patient
                new_graph = []

                # iterate through edges of the current graph
                for edge in graph:

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
                    
                    # add new_edge to new_graph
                    new_graph.append(new_edge)
                
                # add new_graph to new_patient
                new_patient.append(new_graph)
            
            # add new_patient to new_dataset
            new_dataset.append(new_patient)

        # update the dataset
        self.__patients = np.array(new_dataset, dtype='object')

    def _compute_dataset_weighted_edges(self):
        """
        Returns a dictionary with all edges in the input dataset with their respective frequencies.
        The frequency of an edge is incremented by 1/(n_patients * n_graphs) for each occurrence in a patient with n_graphs graphs, where the dataset has n_patients patients.

        Returns:
        - edge_frequencies: dictionary with the frequency of each edge in the dataset.
        """

        # dictionary with the frequency of each edge
        edge_frequencies = {}

        # number of patients in the dataset
        n_patients = self.get_number_patients()

        # iterate through all edges in the input dataset and count their frequency
        for patient in self.__patients:

            # get the normalization term for the frequency of edges found in the current patient
            curr_normalization = n_patients * len(patient)

            # iterate through all edges in the current patient and update their frequency
            for graph in patient:
                for edge in graph:
                    if tuple(edge) in edge_frequencies:
                        edge_frequencies[tuple(edge)] += 1/curr_normalization
                    else:
                        edge_frequencies[tuple(edge)] = 1/curr_normalization
        
        return edge_frequencies

    def _compute_dataset_graph(self):
        """
        Creates a graph with a node for each mutation appearing in graphs in the input dataset and an edge for each edge appearing in graphs in the input dataset.
        Each edge has a weight representing the frequency of the edge in the dataset.
        Remind that a patient potentially has multiple graphs, due to uncertainty, so each edge frequency is normalized also by the number of graphs in each patient.
        
        Returns:
        - dataset_graph: networkx DiGraph with all mutations present in graphs in the input dataset as nodes and all edges appearing in the dataset.
                         Each edge has the 'frequency' attribute storing the frequency of the edge in the dataset.
        """

        # get a dictionary storing all edges in the dataset with their respective frequencies
        edge_frequencies = self.compute_dataset_weighted_edges()

        # create an empty DiGraph
        dataset_graph = nx.DiGraph()

        # add all edges to the created graph
        for edge in edge_frequencies:
            dataset_graph.add_edge(edge[0], edge[1], frequency=edge_frequencies[edge])

        return dataset_graph

    def get_patient(self, index):
        """
        Returns the array of graphs related to the patient in position index in the dataset.

        Parameters:
        - index: index of the patient in the dataset.

        Returns:
        - self.__patients[index]: array of graphs related to the patient in position index in the dataset.
        """

        # return the patient in position index in the dataset
        return self.__patients[index]

    def get_graph(self, patient_index, graph_index):
        """
        Returns the graph in position graph_index of the patient in position patient_index in the dataset.

        Parameters:
        - patient_index: index of the patient in the dataset.
        - graph_index: index of the graph in the patient.

        Returns:
        - self.__patients[patient_index][graph_index]: graph in position graph_index of the patient in position patient_index in the dataset.
        """

        # return the graph in position graph_index of the patient in position patient_index in the dataset
        return self.__patients[patient_index][graph_index]

    def get_number_patients(self):
        """
        Counts the number of patients in the input dataset.

        Returns:
        - len(self.patients): number of patients in the input dataset.
        """

        # number of patients
        return len(self.__patients)
    
    def get_number_graphs(self):
        """
        Counts the number of graphs in the input dataset.

        Returns:
        - n_graphs: number of graphs in the input dataset.
        """

        # number of graphs
        n_graphs = 0

        # sum the number of graphs over patients
        for patient in self.__patients:
            n_graphs += len(patient)

        return n_graphs
    
    def min_number_graphs(self):
        """
        Computes the minimum number of graphs in a patient in the dataset.

        Returns:
        - min_graphs: minimum number of graphs in a patient in the dataset.
        """

        # minimum number of graphs in a patient
        min_graphs = len(self.__patients[0])

        # iterate through all patients and update the minimum number of graphs
        for patient in self.__patients:
            curr_n_graphs = len(patient)
            if curr_n_graphs < min_graphs:
                min_graphs = curr_n_graphs
        
        return min_graphs

    def max_number_graphs(self):
        """
        Computes the maximum number of graphs in a patient in the dataset.

        Returns:
        - max_graphs: maximum number of graphs in a patient in the dataset.
        """

        # maximum number of graphs in a patient
        max_graphs = len(self.__patients[0])

        # iterate through all patients and update the maximum number of graphs
        for patient in self.__patients:
            curr_n_graphs = len(patient)
            if curr_n_graphs > max_graphs:
                max_graphs = curr_n_graphs
        
        return max_graphs

    def avg_number_graphs(self):
        """
        Computes the average number of graphs in a patient in the dataset.

        Returns:
        - avg_graphs: average number of graphs in a patient in the dataset.
        """

        # average number of graphs in a patient in the dataset
        return self.get_number_graphs()/self.get_number_patients()

    def median_number_graphs(self):
        """
        Computes the median number of graphs in a patient in the dataset.

        Returns:
        - median_graphs: median number of graphs in a patient in the dataset.
        """

        # list with the number of graphs in each patient
        n_graphs_list = []

        # iterate through all patients and store the number of graphs in each patient
        for patient in self.__patients:
            n_graphs_list.append(len(patient))
        
        # return the median number of graphs in a patient
        return np.median(n_graphs_list)

    def get_mutations(self):
        """
        Returns the set of mutations that label nodes in the input dataset of graphs.

        Returns:
        - mutation_set: set with all mutation names that label nodes of the input graphs.
        """

        # list with mutation names
        mutation_set = set()

        # iterate through all mutations in the input dataset and add the new ones
        for patient in self.__patients:
            for graph in patient:
                for edge in graph:
                    for mut in edge:
                        if mut not in mutation_set:
                            mutation_set.add(mut)
        
        return mutation_set

    def get_digraph(self, patient_index, graph_index):
        """
        Returns the networkx digraph related to the graph in position graph_index of the patient in position patient_index in the dataset.

        Parameters:
        - patient_index: index of the patient in the dataset.
        - graph_index: index of the graph in the patient.

        Returns:
        - nx.DiGraph(self.__patients[patient_index][graph_index]): networkx digraph related to the graph in position graph_index of the patient in position patient_index in the dataset.
        """

        # create and return a DiGraph with the edges of the graph in the desired position of the dataset
        return nx.DiGraph(self.__patients[patient_index][graph_index])
    
    def get_dataset_graph(self):
        """
        Returns the graph representing the whole dataset.

        Returns:
        - self.__dataset_graph: graph representing the whole dataset.
        """

        # return the graph representing the whole dataset
        return self.__dataset_graph

    def plot_digraph(self, patient_index, graph_index, with_labels=False):
        """
        Plots the graph in position graph_index of the patient in position patient_index in the dataset.

        Parameters:
        - patient_index: index of the patient in the dataset.
        - graph_index: index of the graph in the patient.
        - with_labels: boolean that indicates if the graph should be plotted with labels or not.
        """

        # create a DiGraph with the graph in the desired position of the dataset
        graph = self.get_digraph(patient_index, graph_index)

        # plot the graph
        nx.draw(graph, with_labels=with_labels)
        plt.show()

    def plot_dataset_graph(self):
        """
        Plots the graph representing the whole dataset.
        """

        # plot the graph representing the whole dataset
        nx.draw(self.__dataset_graph, with_labels=True)
        plt.show()
      