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
        """解析反应物和生成物，添加错误处理"""
        try:
            # 检查必要的输入数据
            if not self.reaction_smiles or self.reaction_smiles.strip() == '':
                raise ValueError("SMILES字符串为空")
            
            if '>>' not in self.reaction_smiles:
                raise ValueError("SMILES字符串格式不正确，缺少'>>'分隔符")
            
            # 解析SMILES
            reactants_smiles, products_smiles = self.reaction_smiles.split('>>')
            reactants_smiles_list = reactants_smiles.split('.')
            products_smiles_list = products_smiles.split('.')
            
            # 解析反应方程式
            if not self.reaction_equation or ' = ' not in self.reaction_equation:
                # 如果没有反应方程式，创建默认名称
                reactants_names_list = [f'reactant_{i+1}' for i in range(len(reactants_smiles_list))]
                products_names_list = [f'product_{i+1}' for i in range(len(products_smiles_list))]
            else:
                reactants_names, products_names = self.reaction_equation.split(' = ')
                reactants_names_list = reactants_names.split(' + ')
                products_names_list = products_names.split(' + ')
            
            # 获取系数
            reactants_coef_list = [self.get_cpd_coef(item) for item in reactants_names_list]
            products_coef_list = [self.get_cpd_coef(item) for item in products_names_list]
            
            # 解析ChEBI参考
            if not self.reaction_equation_ref_chebi or ' = ' not in self.reaction_equation_ref_chebi:
                # 如果没有ChEBI参考，创建默认值
                reactants_ref_chebi_list = ['' for _ in range(len(reactants_smiles_list))]
                products_ref_chebi_list = ['' for _ in range(len(products_smiles_list))]
            else:
                reactants_ref_chebi, products_ref_chebi = self.reaction_equation_ref_chebi.split(' = ')
                reactants_ref_chebi_list = reactants_ref_chebi.split(' + ')
                products_ref_chebi_list = products_ref_chebi.split(' + ')
            
            # 创建分子对象，添加错误处理
            self.reactants = []
            for smiles, name, ref_chebi, ref_coef in zip(reactants_smiles_list, reactants_names_list, reactants_ref_chebi_list, reactants_coef_list):
                try:
                    mol = Molecule(cpd_smiles=smiles.strip(), cpd_name=name.strip(), cpd_ref_chebi=ref_chebi.strip(), cpd_num=ref_coef)
                    self.reactants.append(mol)
                except Exception as e:
                    print(f"警告: 无法创建反应物分子 {name} (SMILES: {smiles}): {e}")
            
            self.products = []
            for smiles, name, ref_chebi, ref_coef in zip(products_smiles_list, products_names_list, products_ref_chebi_list, products_coef_list):
                try:
                    mol = Molecule(cpd_smiles=smiles.strip(), cpd_name=name.strip(), cpd_ref_chebi=ref_chebi.strip(), cpd_num=ref_coef)
                    self.products.append(mol)
                except Exception as e:
                    print(f"警告: 无法创建产物分子 {name} (SMILES: {smiles}): {e}")
                    
        except Exception as e:
            # 如果解析失败，创建空的反应物和产物列表
            print(f"解析反应失败 {self.reaction_id}: {e}")
            self.reactants = []
            self.products = []
    

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
    
    
    def save_json_file(self, file_path, overwrite: bool = True):
        """将 Reaction 对象序列化为 JSON 文件"""
        fileTool.write_json_file(self, file_path, overwrite=overwrite)
        
    
    
if __name__ == '__main__':
    rxn_id = 'rxn1'
    rxn_smiles = '[H]O[H].[NH3+][C@H](COP([O-])([O-])=O)C([O-])=O>>[NH3+][C@H](CO)C([O-])=O.OP([O-])([O-])=O'
    rxn_equation = 'H2O + O-phospho-D-serine = D-serine + phosphate'
    rxn_equation_ref_chebi = 'CHEBI:15377 + CHEBI:58680 = CHEBI:35247 + CHEBI:43474'
    
    reaction = Reaction(rxn_smiles, rxn_equation, rxn_equation_ref_chebi, rxn_id=rxn_id, rxn_ec='1.1.1.1')
    reaction.save_json_file(f'{cfg.SAMPLE_DIR}/rxn_sample1.json')
    print(reaction.to_json())