from rdkit import Chem #chem refers to the rdkit chemistry module, has functions for opearting with molecules
from rdkit.Chem import Draw
import argparse

def smiles_chem_png(smiles_mol, image_file='chem_name', show=False): #smiles is the smile notation, this function takes the string into an image
    molecule = Chem.MolFromSmiles(smiles_mol) #MolFromSmiles is a function inside the Chem module. Converts the string to object in rdkit
    if molecule is None:
        print("Check SMILE for invalid notation")
        return

    molecule_image = Draw.MolToImage(molecule, size = (300,300)) #producing image of the molecule, size represents pixels

    molecule_image.save(image_file) #saves the image using the chem name given to it
    print(f"Molecule Image saved as {image_file}")

    if show:
        molecule_image.show() #produces the image

if __name__ == "__main__":
    molecular_parser = argparse.ArgumentParser(description = "SMILES Notation Visualiser") #creates a command line interface
    molecular_parser.add_argument("smiles_mol", type=str, help="SMILES string for the molecule")
    molecular_parser.add_argument("--image", type=str, default="SMILE.png", help="file name")
    molecular_parser.add_argument("--show", action="store_true", help="Display image")

    every_arg = molecular_parser.parse_args() #collects all the arguements and stores them in every_arg
    smiles_chem_png(every_arg.smiles_mol, every_arg.image, every_arg.show)

    #So we have created an interface, in this interface whatever the user puts first is taken as the smiles string, whatever they put second
    #is taken as the name of the image that is produced from the smiles string, and whatever they put third will determine whether to show the
    #image or not. Python then takes these three arguments into our smiles_chem_png user defined function,
    #so in this user defined function smiles = args.smiles, output_file = args.output, show = args.show

    #For example, the smiles_mol = args.smiles_mol --> "CCO" (an example SMILE notation that could be written in the interface)
    #image_file = args.image --> "ethanol.png" (or molecule.png if user doesn't specify name of file)
    #sho = args.show --> True if --show is included in the interface
    
