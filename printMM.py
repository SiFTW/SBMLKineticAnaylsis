#!/usr/bin/env python
## 
## @file    printMath.py
## @brief   Prints whether or not a model contains michaelis mention kinetics for an entire folder of SBML models
## @author  Ben Bornstein
## @author  Sarah Keating
## @author  Simon Mitchell
## 
## <!--------------------------------------------------------------------------
## This sample program is distributed under a different license than the rest
## of libSBML.  This program uses the open-source MIT license, as follows:
##
## Copyright (c) 2013-2014 by the California Institute of Technology
## (California, USA), the European Bioinformatics Institute (EMBL-EBI, UK)
## and the University of Heidelberg (Germany), with support from the National
## Institutes of Health (USA) under grant R01GM070923.  All rights reserved.
##
## Permission is hereby granted, free of charge, to any person obtaining a
## copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
##
## Neither the name of the California Institute of Technology (Caltech), nor
## of the European Bioinformatics Institute (EMBL-EBI), nor of the University
## of Heidelberg, nor the names of any contributors, may be used to endorse
## or promote products derived from this software without specific prior
## written permission.
## ------------------------------------------------------------------------ -->
## 


import sys
import os.path
from libsbml import *

#returns true is a function contains michaelis menten operations
def isMichaelisMenten(fd):
    if fd.isSetMath():


        math = fd.getMath();

        # Print function body. 
        if not (math.getNumChildren() == 0):
            math = math.getChild(math.getNumChildren() - 1);
            
            #function searches for a MM style equation defined as being a division node with a plus and a multiplication node as its leaves with the same terms on one side of the plus as one side of the multiplication.
            def dfs(formula):
                numberOfChildren=formula.getNumChildren()
                for i in range(numberOfChildren):
                    #If this is a division node check whether the ones below it are a plus and a multiplication
                    if(formula.getType()==AST_DIVIDE):
                    #if one of the children is a plus
                    
                        if(formula.getChild(i).getType()==AST_PLUS):
                        #check if the other is a multiplication
                            if(formula.getChild((i+1)%2).getType()==AST_TIMES):
                                #if it is a multiplication chcek whether with a plus
                                leftPlusVar=formulaToString(formula.getChild(i).getLeftChild())
                                rightPlusVar=formulaToString(formula.getChild(i).getRightChild())
                                leftMultiVar=formulaToString(formula.getChild((i+1)%2).getLeftChild())
                                rightMultiVar=formulaToString(formula.getChild((i+1)%2).getRightChild())
                                varsInBoth=set([leftPlusVar,rightPlusVar])&set([leftMultiVar,rightMultiVar])
                                if not (len(varsInBoth)==0):

                                    return True
                                    dfs(formula.getChild(i))
                                    
            if dfs(math):
                print("MM Function Found in " + fd.getId()+".\n"+formulaToString(math));
                return True


def checkForMMInModel(m):
    try:
        for n in range(0,m.getNumFunctionDefinitions()):
            if isMichaelisMenten(m.getFunctionDefinition(n)):
                return True
        for n in range(m.getNumReactions()):
            if isMichaelisMenten(m.getReaction(n).getKineticLaw()):
                return True
        return False
    except:
        print "error with file."
        


def main (args):
  """Usage: printMath filename
  """
  if (len(args) != 2):
      print("\n" + "Usage: printMath filename" + "\n" + "\n");
      return 1;
  
  foldername = args[1];
  outfile = open('MMAnalysis'+foldername+'.txt','w')
  outfile.write("Filename,Contains MM,numReactions,NumSpecies,numParameters\n")
  count=0
  countOfMM=0
  for filename in os.listdir(foldername):

    fullFileName=foldername+os.path.sep+filename
    print fullFileName
    
    document = readSBML(fullFileName);
    count+=1
    if (document.getNumErrors() > 0):
        #print the errors but keep going in case they are warnings
        #TODO: better handling of warnings vs errors
        print("Encountered the following SBML errors:" + "\n");
        print("In filename: "+filename)
        document.printErrors();
        
    model = document.getModel();
          
    if (model == None):
      print("No model present." + "\n");
    else:

        if checkForMMInModel(model):
            outfile.write(filename+",1,"+str(model.getNumSpecies())+","+str(model.getNumReactions())+","+str(model.getNumParameters())+"\n")
            countOfMM+=1
        else:
            outfile.write(filename+",0,"+str(model.getNumSpecies())+","+str(model.getNumReactions())+","+str(model.getNumParameters())+"\n")
  outfile.write('total:'+str(countOfMM)+'/'+str(count))
  outfile.close()
if __name__ == '__main__':
  main(sys.argv)  
