import unittest
import orange

class AZorngTestUtil(unittest.TestCase):

    def randSamp(self, inData, trainingFrac):
        """Use random sampling to partition inData into a training and a test set.
           The seed seem to be reset and the same data partitioning is obtained each time
           the same data is partitioned. """

        selection = orange.MakeRandomIndices2(inData, trainingFrac)
        train_data = inData.select(selection, 0)
        test_data = inData.select(selection, 1)

        return train_data, test_data

##scPA Removed because we are using the evalUtilities.getClassificationAccuracy
#    def getAccuracy(self, learner, test_data):
#        """Return the accuracy (# of correctly classified examples/# of examples) of the learner on test_data. """
#        correct = 0.0
#        for example in test_data:
#            if learner(example) == example.getclass():
#               correct += 1
#        acc = correct/len(test_data)
#        return acc
##ecPA
