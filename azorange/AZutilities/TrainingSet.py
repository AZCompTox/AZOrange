"""
code to do with parsing training set files, and some functions to help you use them.
"""

import string, sys


class TrainingSet:
    """Stores the information from a training set

    The two fields are:
       id_list -- list of compound identifiers, one for each row
       descr_names -- names for each descriptor column
       data_table -- training values for each compound, stored as a
          2D array where the ith row is for the ith id_list compound
          and the elements in the row are the descriptor values
    """
    def __init__(self, smiles_list, id_list, measured_list, descr_names, data_table):
        assert len(smiles_list) == len(id_list)
        assert len(id_list) == len(data_table)
        if measured_list is not None:
            assert len(measured_list) == len(data_table)
        assert len(descr_names) > 0
        self.smiles_list = smiles_list
        self.id_list = id_list
        self.descr_names = descr_names
        self.data_table = data_table
        self.measured_list = measured_list
        
    def __len__(self):
        return len(self.data_table)
    
    def find_compound(self, id, smiles):
        i = -1
        try:
            i = self.id_list.index(id)
        except ValueError:
            pass
        
        if i == -1:
            try:
                i = self.smiles_list.index(smiles)
            except ValueError:
                pass
        return i

                

def read_training_set(infile):
    """training file, num_cols -> returns a TrainingSet

    The second input parameter, 'num_cols' tells how many descriptor
    values to use in a row.  The first column is the id and is not
    counted as a descriptor value.  If num_cols is None then all columns
    will be used.

    """
    line = infile.readline().strip()  # descriptions in the first line
    if "\t" not in line:
        name = getattr(infile, "name", "<unknown>")
        raise TypeError("training set file %r is not tab-delimited" % name)
    
    # First line is:
    # SMILES\t id \t word1\tword2\t...
    # and I would like it to have a \t instead of a " " after the 'id'
    # id\tword1\tword2\t...
    # (it could also be mnumber")
    #if line[:3] == "id ":
        #line = "id\t" + line[3:]
    #elif line[:8] == "mnumber ":
        #line = "mnumber\t" + line[8:]

    words = string.split(line, "\t")
    assert words[0] in ("SMILES"), \
          "First word of first line of the training set must be 'SMILES', not %s" % (repr(words[0]), )
    assert words[1] in ("id", "mnumber", "ID"), \
          "Second word of first line of the training set must be 'id' or 'mnumber' or 'ID', not %s" % (repr(words[1]), )
    measured = False
    if words[-1] in ("measured", "target", "Measured", "Target"):
        measured = True
    
    if measured:
        descr_names = words[2:-1]    
    else:
        descr_names = words[2:]
    
    assert len(descr_names) > 0, "training set file contains no descriptor names in first row %s" % line

    # All the rest of the fields look like:
    # SMILES word float float ...
    # where the 'word' is the identifier
    smiles_list = []
    id_list = []
    measured_list = None
    if measured:
        measured_list = []
    table = []
    linenumber = 0
    for line in infile.xreadlines():
        line = line.strip()
        if (line == ""):
            linenumber += 1
            continue # skip blank lines
        words = line.split('\t')
        smiles_list.append(words[0]) # first column is SMILES
        id_list.append(words[1])    # second column is the identifiers
        try:
            if measured:
                data = words[2:-1]
                measured_list.append(float(words[-1]))
            else:
                data = words[2:]  # get the requested data fields
            data = map(float, data)     # convert those fields to floats
        except ValueError, details:
            # indicates error in training set
            raise ValueError("error %s in training set at line %s: %s" % (details, linenumber, line))
             
        table.append(data)
        linenumber += 1

    import numpy
    table = numpy.array(table, numpy.float)
    return TrainingSet(smiles_list, id_list, measured_list, descr_names, table)

#### Some utility functions for searching the TrainingSet
def compute_scaled_norm(v1, v2, sd):
    assert len(v1) == len(v2) == len(sd), "vectors must be the same length"
    
    dist = 0.0
    for i in range(len(v1)):
        dist = dist + ((v2[i]-v1[i])/sd[i])**2
    return (dist ** (0.5)) * ((15.0/len(sd))**(0.5))

def find_closest(training_set, descr_values, scaling):
    # Find the closest value in the training set
    best_dist = None
    best_id = None
    for i in range(len(training_set)):
        row = training_set.data_table[i]

        dist = compute_scaled_norm(row, descr_values, scaling)
        if best_dist is None or dist < best_dist:
            best_id = training_set.id_list[i]
            best_dist = dist
    return best_id, best_dist



               
