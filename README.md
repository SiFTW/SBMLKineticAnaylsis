#Michaelis Menten in SBML

This script will analyze a folder of Systems Biology Markup Language (SBML) computational models and identify how many of them contain Michaelis Menten reactions. I have run this on the entire of BioModels to generate statistics for the usage of Michaelis Menten equations.

Adapted from one of the libSBML examples.

Requires [LibSBML](http://sbml.org/Software/libSBML/)

printMM.py will analyze every SBML model in a folder. printMMFromList.py takes a list of file names and will only analyze those files (useful for pre-processing to remove all models without enzymes).