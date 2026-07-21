import sys,os
sys.path.append(os.path.dirname(os.path.realpath('__file__')))
sys.path.append('../../')
from rxnrecer.config import config as cfg
from rxnrecer.utils import file_utils as fileTool
import json
from rxnrecer.lib.rxn.Molecule import Molecule
    

class Reaction:
    def __init__(self, rxn_smiles, rxn_equation='', rxn_equation_ref_chebi='', rxn_id='', rxn_ec=''):
        self.reaction_id = rxn_id
        self.reaction_smiles = rxn_smiles
        self.reaction_equation = rxn_equation
        self.reaction_equation_ref_chebi = rxn_equation_ref_chebi
        self.reaction_ec = rxn_ec
        
        self.reactants = []
        self.products = []
        self.parse_reaction()
        
    
    
    def get_cpd_coef(self, cpd_string):
        """Return the stoichiometric coefficient for a compound string."""
        

        cpd_string = cpd_string.strip()


        cpd_string_list = cpd_string.split(' ')


        if not cpd_string_list:
            print(f'Error: cpd_string: {cpd_string}')
            return 1


        coef_str = cpd_string_list[0].strip()


        if coef_str.isdigit():
            return int(coef_str)
        

        if coef_str.lower() in ['a', 'an']:
            return 1


        return 1


    
    def parse_reaction(self):
        """Parse reactants and products with defensive error handling."""
        try:

            if not self.reaction_smiles or self.reaction_smiles.strip() == '':
                raise ValueError("SMILES string is empty")
            
            if '>>' not in self.reaction_smiles:
                raise ValueError("Invalid SMILES string: missing '>>' separator")
            

            reactants_smiles, products_smiles = self.reaction_smiles.split('>>')
            reactants_smiles_list = reactants_smiles.split('.')
            products_smiles_list = products_smiles.split('.')
            

            if not self.reaction_equation or ' = ' not in self.reaction_equation:

                reactants_names_list = [f'reactant_{i+1}' for i in range(len(reactants_smiles_list))]
                products_names_list = [f'product_{i+1}' for i in range(len(products_smiles_list))]
            else:
                reactants_names, products_names = self.reaction_equation.split(' = ')
                reactants_names_list = reactants_names.split(' + ')
                products_names_list = products_names.split(' + ')
            

            reactants_coef_list = [self.get_cpd_coef(item) for item in reactants_names_list]
            products_coef_list = [self.get_cpd_coef(item) for item in products_names_list]
            

            if not self.reaction_equation_ref_chebi or ' = ' not in self.reaction_equation_ref_chebi:

                reactants_ref_chebi_list = ['' for _ in range(len(reactants_smiles_list))]
                products_ref_chebi_list = ['' for _ in range(len(products_smiles_list))]
            else:
                reactants_ref_chebi, products_ref_chebi = self.reaction_equation_ref_chebi.split(' = ')
                reactants_ref_chebi_list = reactants_ref_chebi.split(' + ')
                products_ref_chebi_list = products_ref_chebi.split(' + ')
            

            self.reactants = []
            for smiles, name, ref_chebi, ref_coef in zip(reactants_smiles_list, reactants_names_list, reactants_ref_chebi_list, reactants_coef_list):
                try:
                    mol = Molecule(cpd_smiles=smiles.strip(), cpd_name=name.strip(), cpd_ref_chebi=ref_chebi.strip(), cpd_num=ref_coef)
                    self.reactants.append(mol)
                except Exception as e:
                    print(f"Warning: could not create reactant molecule {name} (SMILES: {smiles}): {e}")
            
            self.products = []
            for smiles, name, ref_chebi, ref_coef in zip(products_smiles_list, products_names_list, products_ref_chebi_list, products_coef_list):
                try:
                    mol = Molecule(cpd_smiles=smiles.strip(), cpd_name=name.strip(), cpd_ref_chebi=ref_chebi.strip(), cpd_num=ref_coef)
                    self.products.append(mol)
                except Exception as e:
                    print(f"Warning: could not create product molecule {name} (SMILES: {smiles}): {e}")
                    
        except Exception as e:

            print(f"Failed to parse reaction {self.reaction_id}: {e}")
            self.reactants = []
            self.products = []
    

    def get_balanced_equation(self, res_type='ref_chebi'):    
        """Return the balanced reaction equation."""
        
        res = ''
        if res_type =='ref_chebi':
            reactants_ref_chebi = [f'{reactant.mol_num} {reactant.cpd_ref_chebi}' if reactant.mol_num > 1 else reactant.cpd_ref_chebi for reactant in self.reactants]
            products_ref_chebi = [f'{product.mol_num} {product.cpd_ref_chebi}' if product.mol_num > 1 else product.cpd_ref_chebi for product in self.products]
            res = f'{" + ".join(reactants_ref_chebi)} = {" + ".join(products_ref_chebi)}'
        else:
            print(f'Unsupported res_type: {res_type}')
            
            res = f'{self.reaction_equation_ref_chebi}'    
            #TODOss
            
            
        return res
        
        
    
    def to_html(self):
        """Render reactants and products as an HTML image block."""
        html_output = "<div style='display: flex; align-items: center; font-size:40px;'>"


        for reactant in self.reactants:
            lb_coef = reactant.mol_num if reactant.mol_num > 1 else ''
            
            html_output += f"<h2 style='font-sze:100px;'>{lb_coef}</h2><img src='{cfg.DIR_PROJECT_ROOT}/{reactant.mol_svg}' style='display:inline-block; margin-right: 10px;'/>"
            html_output += " + "

        html_output = html_output[:-3]
        html_output += " = "


        for product in self.products:
            
            lb_coef = product.mol_num if product.mol_num > 1 else ''
            html_output += f"<h2 style='font-sze:100px;'>{lb_coef}</h2><img src='{cfg.DIR_PROJECT_ROOT}/{product.mol_svg}' style='display:inline-block; margin-right: 10px;'/>"
            html_output += " + "

        html_output = html_output[:-3]
        html_output += "</div>"

        return html_output
    
    def to_dict(self):
        """Convert this reaction to a dictionary."""
        return {
            'reaction_id': self.reaction_id,
            'reaction_smiles': self.reaction_smiles,
            'reaction_equation': self.reaction_equation,
            'reaction_equation_ref_chebi': self.reaction_equation_ref_chebi,
            'reaction_ec': self.reaction_ec,
            'reactants': [reactant.to_dict() for reactant in self.reactants],
            'products': [product.to_dict() for product in self.products]
    }
    
    def to_json(self):
        """Serialize this reaction to JSON."""
        str_json = json.dumps(self.to_dict(), indent=4)
        return str_json
    
    
    def save_json_file(self, file_path, overwrite: bool = True):
        """Serialize this reaction to a JSON file."""
        fileTool.write_json_file(self, file_path, overwrite=overwrite)
        
    
    
if __name__ == '__main__':
    rxn_id = 'rxn1'
    rxn_smiles = '[H]O[H].[NH3+][C@H](COP([O-])([O-])=O)C([O-])=O>>[NH3+][C@H](CO)C([O-])=O.OP([O-])([O-])=O'
    rxn_equation = 'H2O + O-phospho-D-serine = D-serine + phosphate'
    rxn_equation_ref_chebi = 'CHEBI:15377 + CHEBI:58680 = CHEBI:35247 + CHEBI:43474'
    
    reaction = Reaction(rxn_smiles, rxn_equation, rxn_equation_ref_chebi, rxn_id=rxn_id, rxn_ec='1.1.1.1')
    reaction.save_json_file(f'{cfg.SAMPLE_DIR}/rxn_sample1.json')
    print(reaction.to_json())
