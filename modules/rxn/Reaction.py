import sys,os
sys.path.append(os.path.dirname(os.path.realpath('__file__')))
sys.path.append('../../')
from config import conf as cfg
import tools.filetool as fileTool
import base64
import json
from modules.rxn.Molecule import Molecule
    
        


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
        """获取某个化合物的系数"""
        
        # 去除输入字符串的多余空格
        cpd_string = cpd_string.strip()

        # 拆分化合物字符串为列表
        cpd_string_list = cpd_string.split(' ')

        # 如果列表为空，返回 1 并输出错误信息
        if not cpd_string_list:
            print(f'Error: cpd_string: {cpd_string}')
            return 1

        # 获取可能的系数部分
        coef_str = cpd_string_list[0].strip()

        # 判断第一个字符串是否为数字
        if coef_str.isdigit():
            return int(coef_str)
        
        # 处理特殊情况，比如 'a' 或 'an'
        if coef_str.lower() in ['a', 'an']:
            return 1

        # 默认返回 1，保持扩展性
        return 1


    
    def parse_reaction(self):
        """解析反应物和生成物"""
        reactants_smiles, products_smiles = self.reaction_smiles.split('>>')
        reactants_smiles_list = reactants_smiles.split('.')
        products_smiles_list = products_smiles.split('.')
        
        reactants_names, products_names = self.reaction_equation.split(' = ')
        
        reactants_names_list = reactants_names.split(' + ')
        products_names_list = products_names.split(' + ')
        
        
        # print(f'reactants_names_list: {reactants_names_list}, products_names_list: {products_names_list}')
        
        reactants_coef_list = [self.get_cpd_coef(item) for item in reactants_names_list]
        products_coef_list = [self.get_cpd_coef(item)  for item in products_names_list]
        
        reactants_ref_chebi, products_ref_chebi = self.reaction_equation_ref_chebi.split(' = ')
        reactants_ref_chebi_list = reactants_ref_chebi.split(' + ')
        products_ref_chebi_list = products_ref_chebi.split(' + ')

        
        self.reactants = [Molecule(cpd_smiles=smiles, cpd_name=name, cpd_ref_chebi=ref_chebi, cpd_num = ref_coef) 
                          for smiles, name, ref_chebi, ref_coef in zip(reactants_smiles_list, reactants_names_list, reactants_ref_chebi_list, reactants_coef_list)]
        self.products = [Molecule(cpd_smiles=smiles, cpd_name=name, cpd_ref_chebi=ref_chebi, cpd_num = ref_coef) 
                         for smiles, name, ref_chebi, ref_coef in zip(products_smiles_list, products_names_list, products_ref_chebi_list, products_coef_list)]
    

    def get_balanced_equation(self, res_type='ref_chebi'):    
        """获取配平后的反应方程式"""
        
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
        """生成 HTML 来展示反应物和生成物的图像及其链接"""
        html_output = "<div style='display: flex; align-items: center; font-size:40px;'>"

        # 添加反应物图片
        for reactant in self.reactants:
            lb_coef = reactant.mol_num if reactant.mol_num > 1 else '' # 显示分子系数
            
            html_output += f"<h2 style='font-sze:100px;'>{lb_coef}</h2><img src='{cfg.DIR_PROJECT_ROOT}/{reactant.mol_svg}' style='display:inline-block; margin-right: 10px;'/>"
            html_output += " + "

        html_output = html_output[:-3]  # 移除最后的加号
        html_output += " = "

        # 添加生成物图片
        for product in self.products:
            
            lb_coef = product.mol_num if product.mol_num > 1 else ''
            html_output += f"<h2 style='font-sze:100px;'>{lb_coef}</h2><img src='{cfg.DIR_PROJECT_ROOT}/{product.mol_svg}' style='display:inline-block; margin-right: 10px;'/>"
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
            'reaction_ec': self.reaction_ec,
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