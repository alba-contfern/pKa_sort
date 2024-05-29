import PySimpleGUI as sg
from rdkit import Chem
import pubchempy as pcp
from pubchempy import get_compounds
from rdkit.Chem import Draw
import pandas as pd
import re
import sys
import traceback
import xml.etree.ElementTree as ET
from typing import Optional
import pandas as pd
import requests

from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True


#function that gets the smile for each molecule
def get_test(compound):
    results = pcp.get_compounds(compound, 'name')
    for compound in results:
        smiles= compound.isomeric_smiles
        mol=Chem.MolFromSmiles(smiles)
        return mol



#We need an entry in Cas or the name but need to specify it in the second entry of the function
# either pka_lookup_pubchem("acetic acid", "Name") or pka_lookup_pubchem("'64-19-7' ","cid")

debug = False

#(c) 2020 khoivan88
def pka_lookup_pubchem(identifier, namespace=None, domain='compound') -> Optional[str]:
    global debug

    if len(sys.argv) == 2 and sys.argv[1] in ['--debug=True', '--debug=true', '--debug', '-d']:
        debug = True

    # if debug:
    #     print(f'In DEBUG mode: {debug}')

    # Identify lookup source (Pubchem in this case)
    lookup_source = 'Pubchem'

    try:
        headers = {
            'user-agent': 'Mozilla/5.0 (X11; CentOS; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/73.0.3683.75 Safari/537.36'}
        # Using pubchem api for python
        # Getting CID number, the result of this, by default is exact match. The result is returned as a list.
        cids = []
        identifier_type = ''

        if not namespace:
            identifier_type = classify(identifier)

            # If the input is inchi, inchikey or smiles (this could be a false smiles):
            if identifier_type in ['smiles', 'inchi', 'inchikey']:
                lookup = pcp.get_cids(identifier, namespace=identifier_type)
                if lookup:
                    cids.append(lookup[0])
            else:
                lookup = pcp.get_cids(identifier, namespace='name')
                if lookup:
                    cids.append(lookup[0])
                    # print(f'namespace from pubchem lookup is: {namespace}')
        elif namespace == 'cas':
            cids = pcp.get_cids(identifier, namespace='name')
        else:
            cids = pcp.get_cids(identifier, namespace=namespace)

        if not cids:
            lookup = pcp.get_cids(identifier, namespace='name')
            if lookup:
                cids.append(lookup[0])

            # cids = pcp.get_cids(identifier, namespace=namespace)
            identifier_type = namespace

        if len(cids) > 0:
            # if Pubchem found the result, get the first result of the list
            cid = cids[0]

            exact_match = True

            # synonyms = []
            synonyms = pcp.get_synonyms(cid)[0]['Synonym'] or []
            
            # Extract CAS number from the list of synonyms
            returned_cas = ''
            for synonym in synonyms:
                cas_nr = re.search(r'^\d{2,7}-\d{2}-\d$', synonym)
                if cas_nr:
                    cas_nr = cas_nr.group()
                    returned_cas = cas_nr
                    break

            lookup_result = pcp.get_properties(['inchi', 'inchikey',
                                        'canonical_smiles', 'isomeric_smiles',
                                        'iupac_name'],
                                cid)

            if identifier_type == 'cas':
                # To double check if the CAS number is correct:
                # using pubchem api, get a list of synonym. The result is a list of dict.
                # choose the first result and check all values for 'Synonym' key:
                exact_match = identifier in synonyms

            elif identifier_type in ['inchi', 'inchikey']:

                if identifier_type == 'inchi':
                    # print(lookup_result[0].get('InChI', False))
                    # print(f'input:\n{identifier}')
                    exact_match = (identifier == lookup_result[0].get('InChI', False))
                
                elif identifier_type == 'inchikey':
                    exact_match = (identifier == lookup_result[0].get('InChIKey', False))

            if not exact_match:
                if debug:
                    print(f'Exact match between input and Pubchem return value? {identifier in synonyms}')
                raise ValueError('This is not an exact match on Pubchem!')


            pka_lookup_result_xml = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/XML?heading=Dissociation+Constants'.format(cid)

            # Get the html request info using CID number from pubchem
            r = requests.get(pka_lookup_result_xml, headers=headers, timeout=15)
            # Check to see if give OK status (200) and not redirect
            if r.status_code == 200 and len(r.history) == 0:
                # print(r.text)
                # Use python XML to parse the return result
                tree = ET.fromstring(r.text)
            
                # Get the XML tree of  only
                info_node = tree.find('.//*{http://pubchem.ncbi.nlm.nih.gov/pug_view}Information')

                # Get the pKa reference:
                original_source = info_node.find('{http://pubchem.ncbi.nlm.nih.gov/pug_view}Reference').text
                # Get the pKa result:
                pka_result = info_node.find('.//*{http://pubchem.ncbi.nlm.nih.gov/pug_view}String').text
                pka_result = re.sub(r'^pKa = ', '', pka_result)    # remove 'pka = ' part out of the string answer
                # print(pka_result)
                # print(original_source)
                # print(lookup_result)

                core_result = {
                    'source': lookup_source,
                    'Pubchem_CID': str(cid),
                    'pKa': pka_result,
                    'reference': original_source,
                    'Substance_CASRN': returned_cas,
                }
                extra_info = lookup_result[0]
                extra_info.pop('CID', None)    # Remove 'CID': ... from lookup_result[0]

                # Merge 2 dict: https://treyhunner.com/2016/02/how-to-merge-dictionaries-in-python/
                result = {**core_result, **extra_info}
                # Rename some keys in the dict
                s = pd.Series(result)
                s = s.rename({
                    'CanonicalSMILES': 'Canonical_SMILES',
                    'IsomericSMILES': 'Isomeric_SMILES',
                    'IUPACName': 'IUPAC_Name'
                })
                result = s.to_dict()            
                return result

            else:
                raise RuntimeError('pKa not found in Pubchem.')
    
        else:
            raise RuntimeError('Compound not found in Pubchem.')

    except Exception as error:
        if debug:
            traceback_str = ''.join(traceback.format_exception(etype=type(error), value=error, tb=error.__traceback__))
            #print(traceback_str)

        return None



#function that sort the list of molecule and return a list of set in order for the pka
def pka_increasing(list):
    dict={}
    for i in range (len(list)):
        pka=pka_lookup_pubchem(list[i],'name')
        dict[pka['pKa'][0:4]]=list[i]
    molecule_list_pka=dict.items()
    sorted_list = sorted(molecule_list_pka, key=lambda x: float(x[0]))
    return sorted_list
#print(pka_increasing(["aspirin","acetic acid"]))


#generate the name if the molecule from a list: [('3.47', 'aspirin'), ('4.75', 'acetic acid')] (output of pka increasing)
def generate_name(list_of_molecule):
    g=0
    while g ==0:
        name_list = [item[1] for item in list_of_molecule]
        #print(name_list)
        g=1
    return name_list

#get a list of the image in png: [<IPython.core.display.SVG object>, <IPython.core.display.SVG object>]
def molecule_list_image(molecule_list):
    image_list=[]
    for i in range(len(molecule_list)):
        molecule = get_test(molecule_list[i])
        #print(molecule)
        core=Chem.MolFromSmarts('[OH]')
        #print(core)
    #return Draw.MolsToGridImage(mss,highlightAtomLists=[mol.GetSubstructMatch(core) for mol in mss])
        image= Draw.MolsToGridImage([molecule], highlightAtomLists=[molecule.GetSubstructMatch(core)], useSVG=False, returnPNG=True)
        image_list.append(image)
    return image_list
#molecule_list = ["aspirin","acetic acid", "ethanol"]
#print(molecule_list_image(molecule_list))

#function to check the number of molecule to be compare if the entry is right
def number(number):
    try:
        number=int(number)
        test="a"
    except ValueError:
        test="b"
    return test

def molecule_test_pka(molecule):

    molecule = pka_lookup_pubchem(molecule,"name")
    if molecule is None:
        c=1
    else:
        c=0
    return c

#get the number of molecule
number_entry = sg.popup_get_text(("How many molecules do you want to compare?"),title="Number of molecules", background_color="grey")
if number_entry is None:
    sys.exit()
while number(number_entry)=="b": #check if it can be int
    sg.popup_error("The entry is not valid, please enter an integer number") #if it is not an valide int
    number_entry = sg.popup_get_text(("How many molecules do you want to compare?"),title="Number of molecules", background_color="grey")
    if number_entry is None:
        sys.exit()

if int(number_entry)>3:
    answer=sg.popup_yes_no(("If you enter more than 3 molecules, the graphical representation won't work and the program may be slow. Are you sure of your choice?"), title="Attention"
     , background_color="grey")
    #print(answer)
    if answer == "No": #if the number was wrongly choosen
        number_entry= sg.popup_get_text(("How many molecules do you want to compare?"), title="Number of molecules", background_color="grey")
        while number(number_entry)=="b": #recheck if it is an int
            sg.popup_error("The entry is not valid, please enter an integer number")
            number_entry = sg.popup_get_text(("How many molecules do you want to compare?"), title="Number of molecules", background_color="grey")
            if number_entry is None:
                sys.exit()

#get the name of the molecule
molecule_list = []

for i in range(int(number_entry)):
    text= sg.popup_get_text(('What is the name of the molecule?'),title="Name of the molecule", background_color="grey")
    if text is None:
        sys.exit()
    #print(text)
    compound=get_compounds(text, "name")
    #print(compound)
    while not compound:
        sg.popup_error("The entered molecule is not in pubchem database", title="error messsage")
        text= sg.popup_get_text(('Enter a valid name for the molecule'), title="Valid name", background_color="grey")
        compound=get_compounds(text, "name")
        if text is None:
            sys.exit()
    #print(text)
    
    #########################################        
    c = molecule_test_pka(text)
    while c==1:
        sg.popup_error("Sorry, the pKa of this molecule is not in pubchem database", title="error message")
        text=sg.popup_get_text(('Enter the name of another molecule'), title="Valid name", background_color="grey")
        if text is None:
            sys.exit()
        c = molecule_test_pka(text)
    molecule_list.append(text)
    

molecule_in_order = pka_increasing(molecule_list) #that's a list of set
#print(molecule_in_order)

list_in_order = generate_name(molecule_in_order) #same list but only the name
#print(list_in_order)

list_of_image = molecule_list_image(list_in_order)
#print(list_of_image)

def order(list_in_order):
    order = list_in_order[0]
    for i in range(len(list_in_order)-1):
        order = order + " " + ">" + " " + list_in_order[i+1]
    return order



if len(list_in_order)<=3:
    def make_window():
        image_list = list_of_image
        for_layout = []
        text_layout_3=[]
        
        def frame_with_text(text):
            return sg.Frame("", [[text]], pad=(5, 3), expand_x=True, expand_y=True, background_color='white', border_width=0)
        
        sg.theme('DarkGrey4')

        for i in range(len(image_list)):
            text= [sg.Text(f"The pKa of {list_in_order[i]} is {molecule_in_order[i][0]} ", text_color='black', background_color="white", auto_size_text=14)]
            text_layout_3.append(text)
        
        layout_frame1 =  [[ frame_with_text(text_layout_3)], 
        [sg.Frame("Ranking of the acidity", [[frame_with_text(sg.Text(order(list_in_order), text_color='black', background_color="white"))] ], pad=(5, 3), expand_x=True, expand_y=True, title_location=sg.TITLE_LOCATION_TOP,)]]
        
        for i in range(len(image_list)):
            text= [sg.Text(list_in_order[i]+":", text_color='white', font=("Times New Roman", 17), background_color="grey32", auto_size_text=14)]
            for_layout.append(text)
            image = [sg.Image(image_list[i].data,expand_y=True, expand_x=True, size=(600,200), background_color='grey32')]
            for_layout.append(image)


        layout_frame2 = [[ for_layout]]

        layout = [
            [sg.Frame("pKa of the molecule", layout_frame1, size=(700, 250), title_location=sg.TITLE_LOCATION_TOP),
            sg.Frame("Graphic representation of the molecules", layout_frame2, size=(700, 900), title_location=sg.TITLE_LOCATION_TOP)],]


        window = sg.Window("The molecule ranked in order of increasing pKa",layout,size=(1400,900), background_color='grey32', finalize=True)

        return window

else:
    def make_window():
        image_list = list_of_image
        text_layout_3=[]

        def frame_with_text(text):
            return sg.Frame("", [[text]], pad=(5, 3), expand_x=True, expand_y=True, background_color='white', border_width=0)

        sg.theme('DarkGrey4')

        for i in range(len(image_list)):
            text= [sg.Text(f"The pKa of {list_in_order[i]} is {molecule_in_order[i][0]} ", text_color='black', background_color="white", auto_size_text=14, justification='center')]
            text_layout_3.append(text)

        layout_frame2 =  [[ frame_with_text(text_layout_3)]]
        


        layout_frame1 = [[sg.Frame("", [[frame_with_text(sg.Text(order(list_in_order), text_color='black', background_color="white", justification="center"))] ], pad=(2, 5), expand_x=True, expand_y=True, title_location=sg.TITLE_LOCATION_TOP,)]]


        layout = [
            [sg.Frame("Ranking of the molecules with decreasing acidity", layout_frame1, size=(700, 250), title_location=sg.TITLE_LOCATION_TOP),
            sg.Frame("pKa of the molecules", layout_frame2, size=(700, 900), title_location=sg.TITLE_LOCATION_TOP)],]


        window = sg.Window("Acidic comparaison of molecules",layout,size=(1400,900), background_color='grey32', finalize=True)

        return window

def main():

    window = make_window()

    while True:
        event, values = window.read()     # wake every hour
        #print(event, values)
        if event == sg.WIN_CLOSED or event == 'Exit':
            break
        if event == 'Add Item':
            window.metadata += 1
            window.extend_layout(window['-TRACKING SECTION-'], [item_row(window.metadata)])
        elif event == 'Edit Me':
            sg.execute_editor(__file__)
        elif event == 'Version':
            sg.popup_scrolled(__file__, sg.get_versions(), location=window.current_location(), keep_on_top=True, non_blocking=True)
        elif event[0] == '-DEL-':
                window[('-ROW-', event[1])].update(visible=False)

if __name__ == '__main__':
    main()