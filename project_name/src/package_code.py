

from rdkit import Chem
import pubchempy as pcp
from pubchempy import get_compounds

from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

def get_test(compound):
    results = pcp.get_compounds(compound, 'name')
    print (results)
    for compound in results:
        smiles= compound.isomeric_smiles
        mol=Chem.MolFromSmiles(smiles)
        return mol

        
       

print (get_test('glucose'))