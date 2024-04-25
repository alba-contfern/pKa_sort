
from rdkit import Chem
import pubchempy as pcp
from pubchempy import get_compounds
from rdkit.Chem import Draw
import pandas as pd

from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

#function that gets the smile for each molecule
def get_test(compound):
    results = pcp.get_compounds(compound, 'name')
    for compound in results:
        smiles= compound.isomeric_smiles
        mol=Chem.MolFromSmiles(smiles)
        return mol
#function that allows all molecules to be represented at the same time
list=['glucose','aspirin','formaldehyde']
def generate_image(list:list):
    mss=[]
    for index,value in enumerate(list):
        mss.append(get_test(value))
    return Draw.MolsToGridImage(mss)


generate_image(list)