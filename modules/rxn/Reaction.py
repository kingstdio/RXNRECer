import sys,os
sys.path.append(os.path.dirname(os.path.realpath('__file__')))
sys.path.append('../../')
import tools.filetool as fileTool
import base64
import json
from modules.rxn.Molecule import Molecule
    
        


class Reaction:
    def __init__(self, rxn_smiles, rxn_equation='', rxn_equation_ref_chebi='', rxn_id=''):
        self.reaction_id = rxn_id
        self.reaction_smiles = rxn_smiles
        self.reaction_equation = rxn_equation
        self.reaction_equation_ref_chebi = rxn_equation_ref_chebi
        
        self.reactants = []
        self.products = []
        
        
        self.parse_reaction()
        
    
    def parse_reaction(self):
        """解析反应物和生成物"""
        reactants_smiles, products_smiles = self.reaction_smiles.split('>>')
        reactants_smiles_list = reactants_smiles.split('.')
        products_smiles_list = products_smiles.split('.')
        
        reactants_names, products_names = self.reaction_equation.split(' = ')
        reactants_names_list = reactants_names.split(' + ')
        products_names_list = products_names.split(' + ')
        
        reactants_ref_chebi, products_ref_chebi = self.reaction_equation_ref_chebi.split(' = ')
        reactants_ref_chebi_list = reactants_ref_chebi.split(' + ')
        products_ref_chebi_list = products_ref_chebi.split(' + ')

        
        self.reactants = [Molecule(cpd_smiles=smiles, cpd_name=name, cpd_ref_chebi=ref_chebi) 
                          for smiles, name, ref_chebi in zip(reactants_smiles_list, reactants_names_list, reactants_ref_chebi_list)]
        self.products = [Molecule(cpd_smiles=smiles, cpd_name=name, cpd_ref_chebi=ref_chebi) 
                         for smiles, name, ref_chebi in zip(products_smiles_list, products_names_list, products_ref_chebi_list)]
        
    

    

    def to_html(self):
        """生成 HTML 来展示反应物和生成物的图像及其链接"""
        html_output = "<div style='display: flex; align-items: center; font-size:40px;'>"

        # 添加反应物图片
        for reactant in self.reactants:
            svg_data_reactant = base64.b64encode(reactant.mol_svg.encode('utf-8')).decode('utf-8')
            html_output += f"<img src='data:image/svg+xml;base64,{svg_data_reactant}' style='display:inline-block; margin-right: 10px;'/>"
            html_output += " + "

        html_output = html_output[:-3]  # 移除最后的加号
        html_output += " = "

        # 添加生成物图片
        for product in self.products:
            svg_data_product = base64.b64encode(product.mol_svg.encode('utf-8')).decode('utf-8')
            html_output += f"<img src='data:image/svg+xml;base64,{svg_data_product}' style='display:inline-block; margin-right: 10px;'/>"
            html_output += " + "

        html_output = html_output[:-3]  # 移除最后的加号
        html_output += "</div>"

        return html_output
    
    def to_dict(self):
        """将 Reaction 对象转化为字典形式"""
        return {
            'reaction_id': self.reaction_id,
            'reaction_smiles': self.reaction_smiles,
            'reaction_equation': self.reaction_equation,
            'reaction_equation_ref_chebi': self.reaction_equation_ref_chebi,
            'reactants': [reactant.to_dict() for reactant in self.reactants],
            'products': [product.to_dict() for product in self.products]
    }
    
    def to_json(self):
        """将 Reaction 对象序列化为 JSON"""
        
        str_json = json.dumps(self.to_dict(), indent=4)
        
        return str_json
    
    
    def save_json_file(self, file_path):
        """将 Reaction 对象序列化为 JSON 文件"""
        
        if fileTool.checkFileExists_with_dir_make(file_path) == True:
            print(f'File exists: {file_path}')
        with open(file_path, 'w') as f:
            f.write(self.to_json())
    
    
    
if __name__ == '__main__':
    rxn_id = 'rxn1'
    rxn_smiles = '[H]O[H].[NH3+][C@H](COP([O-])([O-])=O)C([O-])=O>>[NH3+][C@H](CO)C([O-])=O.OP([O-])([O-])=O'
    rxn_equation = 'H2O + O-phospho-D-serine = D-serine + phosphate'
    rxn_equation_ref_chebi = 'CHEBI:15377 + CHEBI:58680 = CHEBI:35247 + CHEBI:43474'
    
    reaction = Reaction(rxn_smiles, rxn_equation, rxn_equation_ref_chebi, rxn_id=rxn_id)
    print(reaction.to_json())