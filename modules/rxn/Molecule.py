import base64
from rdkit import Chem
import json
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')    


class Molecule:
    
    def __init__(self, cpd_smiles, cpd_name='', cpd_ref_chebi='', cpd_link='#', cpd_id=''):
        
        self.cpd_id = cpd_id
        self.cpd_smiles = cpd_smiles  # SMILES 字符串
        self.cpd_name = cpd_name      # 化合物名称
        self.cpd_link = cpd_link        # 外部链接
        self.cpd_ref_chebi = cpd_ref_chebi
        self.mol = Chem.MolFromSmiles(self.cpd_smiles)
        self.mol_svg = self.draw_mol_simple()

        
        
    
    def get_cpd_id(self):
        return self.cpd_id
    
    def get_cpd_smiles(self):
        return self.cpd_smiles
    
    def get_cpd_name(self):
        return self.cpd_name
    
    def get_cpd_link(self):
        
        if self.cpd_ref_chebi!='':
            self.cpd_link = f'https://www.ebi.ac.uk/chebi/searchId.do?chebiId={self.cpd_ref_chebi}'
        else:
            self.cpd_link = '#'
        
        return self.cpd_link
    
    def get_mol(self):
        return self.mol
    
    def get_mol_svg(self):
        return self.draw_mol_simple()
    
    
    def write_mol_svg(self, file_path):
        with open(file_path, 'w') as f:
            f.write(self.mol_svg)
    
    
    #region 获取 RDKit 的 SVG 绘图器
    def get_drawer(self, mol_pic_width, mol_pic_height, tag_box_height):
        """获取 RDKit 的 SVG 绘图器"""
        drawer = rdMolDraw2D.MolDraw2DSVG(width=mol_pic_width, 
                                    height=mol_pic_height + tag_box_height + 20 ,    # 加上底部方框的高度
                                    panelWidth=-1,
                                    panelHeight=mol_pic_height + tag_box_height # 加上底部方框的高度
                                    )
        
        
        do = drawer.drawOptions()
        do.addStereoAnnotation = False
        do.explicitMethyl = True
        do.addStereoAnnotation = True  # 显示立体化学注释

        # 控制字体大小
        do.annotationFontScale = 0.8  # 调整注释字体大小
        do.baseFontSize = 0.6         # 调整原子标签字体大小
        do.legendFontSize = 24       # 调整化合物名称的字体大小

        # 设置原子颜色
        do.setAtomPalette({
            6: (204/255, 153/255, 51/255, 1),  # 碳原子颜色
            8: (0/255, 153/255, 102/255, 1), # 氧原子颜色
            7: (51/255, 102/255, 153/255, 1),     # 氮原子颜色
            15: (255/255, 153/255, 51/255, 1) # 磷原子颜色 (例子颜色)
        })
        
        return drawer
    #endregion
    
    def add_mol_name_labele(self, mol_svg, mol_pic_width, mol_pic_height, mol_pic_size, tag_box_height):
        """在 SVG 图片上添加化合物名称及超链接"""
        box_width = int(mol_pic_width*0.8)
        box_height = 40
        box_up_padding = 20
        
        x_box = int((mol_pic_width - box_width)/2)
        y_box = int(mol_pic_height + tag_box_height/2 + box_up_padding)
        
        x_text = int(x_box + box_width/2)
        y_text = int(y_box  +box_height/2 + 5)
        
        cpd_name_len = len(self.cpd_name)
        label_font_size = int(1.6*box_width/cpd_name_len)
        
        if label_font_size > 12:
            label_font_size = 12

        self.cpd_link = self.get_cpd_link()
        # 调整底部方框的位置和大小
        additional_svg = f'''
        <a href="{self.cpd_link}" target="_blank" style="cursor:pointer;">
        <g>
            <!-- 底部方框 -->
            <rect x="{x_box}" y="{y_box}" width="{box_width}" height="{box_height}" rx="10" ry="10" style="fill:#0066CC;" />
            <!-- 化合物名称及超链接 -->
            <text x="{x_text }" y="{y_text}" font-size="{label_font_size}"  font-weight="bold" text-anchor="middle" fill="#FFFFFF">{self.cpd_name}</text>
        </g>
        </a>
        '''

        # 将方框及文字加到生成的 SVG 中，放在关闭标签 `</svg>` 之前
        mol_svg = mol_svg.replace('</svg>', additional_svg + '</svg>')
        
        return mol_svg
    
    
    def draw_mol_simple(self):
        mol = self.mol

        mol_pic_width = 200
        mol_pic_height = 150
        tag_box_height = 130
        mol_pic_size = mol.GetNumAtoms() + mol.GetNumBonds()
        
        
        if mol_pic_size <50:
            scale_param = 2.5
        else:
            scale_param = 2
            
        mol_pic_width = int(mol_pic_width + mol_pic_size *scale_param )
        
        # 初始化 RDKit 的 SVG 绘图器
        drawer = self.get_drawer(mol_pic_width=mol_pic_width, 
                                 mol_pic_height=mol_pic_height, 
                                 tag_box_height=tag_box_height
                                 
                                 )
        
        # 显示化合物名称到反应结构图上
        Chem.Draw.PrepareAndDrawMolecule(drawer, mol, kekulize=False)

        # 完成绘制并获取 SVG 文本
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        if self.cpd_name != '':
            svg = self.add_mol_name_labele(mol_svg=svg,
                                        mol_pic_width=mol_pic_width,
                                        mol_pic_height=mol_pic_height,
                                        mol_pic_size=mol_pic_size, 
                                        tag_box_height= tag_box_height)
            
            #TODO 将svg写入文件，返回链接    
            
            

        return svg

    
    def to_html(self):
        """生成 HTML 来展示反应物和生成物的图像及其链接"""
        html_output = "<div style='display: flex; align-items: center;'>"
        svg_data_reactant = base64.b64encode(self.mol_svg.encode('utf-8')).decode('utf-8')
        html_output += f"<img src='data:image/svg+xml;base64,{svg_data_reactant}' style='display:inline-block'/>"
        html_output += "</div>"

        return html_output
    
    def to_dict(self):
        return {
            'cpd_id'     : self.cpd_id,
            'cpd_smiles': self.cpd_smiles,
            'cpd_name': self.cpd_name,
            'cpd_ref_chebi': self.cpd_ref_chebi,
            'cpd_link':self.cpd_link,
            'mol_svg': self.mol_svg
        }
    

    def to_json(self):
        """将 Molecule 对象序列化为 JSON"""
        return json.dumps(self.to_dict(), indent=4)
    
    
if __name__ == '__main__':
    smiles = 'C1=CC=C(C=C1)C(C)C(=O)O'
    name = 'ethanol'
    ref = 'CHEBI:63637'
    mol = Molecule(smiles, name, ref)
    print(mol.to_json())