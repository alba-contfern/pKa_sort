from rdkit import Chem
from PubChemPy import get_cids
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

get_cids(c('Aspirin'))
    
    
from rdkit import Chem
import pubchempy as pcp
from pubchempy import get_compounds

from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

def get_smiles(compound):
    c=pcp.get_compounds (compound,'name')
    d=pcp.get_compounds (compound,'smiles')
    return c,d
    

print (get_smiles('aspirin'))