import pandas as pd
from bioservices import KEGG
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import io
import matplotlib.pyplot as plt
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
from PIL import ImageOps

def get_pubchem_id_from_kegg(kegg_cid):
    '''
    Retrieves the PubChem Substance ID corresponding to a given KEGG compound ID.
    This function uses the bioservices package to access the KEGG database,
    fetches compound information, and extracts the PubChem ID from the available database links.

    :param kegg_cid: KEGG compound ID
    :return: PubChem ID as a string
    '''
    kegg = KEGG()
    compound_info = kegg.get(kegg_cid)
    parsed_info = kegg.parse(compound_info)
    
    # Extracting the PubChem ID
    for db_entry in parsed_info.get('DBLINKS', {}):
        if 'PubChem' in db_entry:
            return parsed_info['DBLINKS']['PubChem']

    return None


def get_smiles_from_pubchem(pubchem_id):
    '''
    Obtains the SMILES (Simplified Molecular Input Line Entry System) representation 
    of a chemical compound using its PubChem Substance ID to retreive the CID, then using that 
    to retrieve the SMILES string, which is a textual 
    representation of the compound's structure.

    :param pubchem_id: PubChem Substance ID of the compound
    :return: SMILES string
    '''
    
    substance = pcp.Substance.from_sid(pubchem_id)
    compound = pcp.Compound.from_cid(substance.cids[0])
    return compound.isomeric_smiles

def get_compounds_from_module(module_id):
    '''
    Fetches a list of compounds associated with a specific KEGG module ID.
    This function communicates with the KEGG database to obtain information 
    about a module and extracts the IDs and names of compounds involved in that module.

    :param module_id: KEGG module ID
    :return: List of tuples containing compound names and their KEGG IDs
    '''    
    kegg = KEGG()
    module_info = kegg.get(module_id)
    parsed_module = kegg.parse(module_info)
    
    # Extracting the compound IDs and names
    compound_data = []
    if 'COMPOUND' in parsed_module:
        for cid, name in parsed_module['COMPOUND'].items():
            compound_data.append((name, cid))
    
    return compound_data

def draw_molecule_on_plot(ax, smiles, x, y, scale=0.1):
    '''
    Renders a molecule's image on a given Matplotlib plot (`ax`) at specified coordinates (`x`, `y`).
    The molecule is represented by its SMILES string. This function uses RDKit to convert the SMILES 
    string to a molecular image and then plots this image on the provided Matplotlib Axes object at 
    the designated location and scale.

    :param ax: Matplotlib Axes object
    :param smiles: SMILES string of the molecule
    :param x: X-coordinate on the plot
    :param y: Y-coordinate on the plot
    :param scale: Scale for the size of the molecule image
    '''
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    # Create an image of the molecule
    mol_img = Draw.MolToImage(mol)

    # Convert the RDKit image to a format compatible with Matplotlib
    buf = io.BytesIO()
    mol_img.save(buf, format='PNG')
    buf.seek(0)
    im = Image.open(buf)

    # Define the extent of the image on the plot
    im_extent = (x - scale, x + scale, y - scale, y + scale)

    # Overlay the molecule image on the plot
    ax.imshow(im, aspect='auto', extent=im_extent, zorder=1)
