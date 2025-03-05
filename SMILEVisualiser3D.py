from rdkit import Chem #chem refers to the rdkit chemistry module, has functions for opearting with molecules
from rdkit.Chem import Draw, AllChem
import argparse
import py3Dmol
import sys
from IPython.display import display, HTML

def is_run_in_nb(): #checking to see if we are running in notebook
    try:
        from IPython import get_ipython
        return get_ipython() is not None
    except ImportError:
        return False

def smiles_chem_png(smiles_mol, image_file='chem_name', show=False, view3d=False): #smiles is the smile notation, this function takes the string into an image
    molecule = Chem.MolFromSmiles(smiles_mol) #MolFromSmiles is a function inside the Chem module. Converts the string to object in rdkit
    if molecule is None:
        print("Check SMILE for invalid notation")
        return

#Visualing in 2D
    molecule_image = Draw.MolToImage(molecule, size = (300,300)) #producing image of the molecule, size represents pixels

    molecule_image.save(image_file) #saves the image using the chem name given to it
    print(f"Molecule Image saved as {image_file}")

    if show:
        molecule_image.show() #produces the image

#Visualising in 3D
    if view3d:
        molecule_3d = Chem.AddHs(molecule) #adds hydrogen atoms for the 3d structure
        result = AllChem.EmbedMolecule(molecule_3d, AllChem.ETKDG()) #Embed computes 3D atomic positions, ETKDG ensures realistic conformations
        if result !=0:
            print("Failed to generate 3D coords, try different molecule")
            return

        molecule_block = Chem.MolToMolBlock(molecule_3d) #converting the molecule into a MolBlock format
        molecule_pdb = Chem.MolToPDBBlock(molecule_3d) #PDB Format

        if is_run_in_nb(): #if it is running in notebook
            viewer = py3Dmol.view(width=500, height=500) #creates a 3d canvas, with window dimensions 500 by 500 in pixels, initially an empty 3D space
            viewer.addModel(molecule_block, "mol") #molecule_block contains the 3d structure of the molecule in MolBlock format, "mol" is file format
            #this loads the molecule into the 3D Canvas
            viewer.setStyle({"stick":{}}) #style of molecule display, others include sphere, cartoon, line
            viewer.zoomTo() #allows molecule to automatically fit in the window
            display(HTML(viewer._make_html())) # displays the visualisations
        else:
            pdb_file_name = image_file.replace(".png", ".pdb")
            if molecule_pdb:
                with open(pdb_file_name, "w") as pdb_file:
                    pdb_file.write(molecule_pdb)
                print(f"3D structure saved as {pdb_file_name}. Use another viewer")
            else:
                print("Error: Could not generate PDB file")

if __name__ == "__main__":
    molecular_parser = argparse.ArgumentParser(description = "SMILES Notation Visualiser") #creates a command line interface
    molecular_parser.add_argument("smiles_mol", type=str, help="SMILES string for the molecule")
    molecular_parser.add_argument("--image", type=str, default="SMILE.png", help="file name")
    molecular_parser.add_argument("--show", action="store_true", help="Display image")
    molecular_parser.add_argument("--view3d", action="store_true", help="3D Visualisation")

    every_arg = molecular_parser.parse_args() #collects all the arguements and stores them in every_arg
    smiles_chem_png(every_arg.smiles_mol, every_arg.image, every_arg.show, every_arg.view3d)

    #So we have created an interface, in this interface whatever the user puts first is taken as the smiles string, whatever they put second
    #is taken as the name of the image that is produced from the smiles string, and whatever they put third will determine whether to show the
    #image or not. Python then takes these three arguments into our smiles_chem_png user defined function,
    #so in this user defined function smiles = args.smiles, output_file = args.output, show = args.show

    #For example, the smiles_mol = args.smiles_mol --> "CCO" (an example SMILE notation that could be written in the interface)
    #image_file = args.image --> "ethanol.png" (or molecule.png if user doesn't specify name of file)
    #sho = args.show --> True if --show is included in the interface

    ###IMPORTANT - WHEN RUNNING THIS IN JUPYTER, KERNERL MIGHT BE STUCK, TRY "jupyter notebook --no-browser" in terminal and open link in chrome
    ###IMPORTANT - If you want the 2D image to show in Jupyter, change if show:
                                                                            #molecule_image.show() #produces the image
                                                                #to     from IPython.display import display
                                                                        #if show:
                                                                            #display(molecule_image) This will show 2D image in Jupyter aswell 
    

