This dataset focus on the task of ranking documents by their reading difficulty. It is composed of 490 documents. Using the CrowdFlower crowdsourcing platform,annotators in the United States and Canada
were shown representative passages from randomly selected
pairs of these documents, and asked to decide which of the
two texts was more challenging to read and understand.

If you want to test your method on this dataset, please load the paired comparison data (i.e. readdataset.mat).  The first column in readdataset.mat is the document ID which is more challenging to read and the second column is the document ID that is less challenging.

groundtruthdata.mat is the gold-standard reading difficulty
level from 1 to 12 of each documents.