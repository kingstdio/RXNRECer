import sys,os
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, f'{project_root}/../')
from rxnrecer.config import config as cfg
from rdkit import Chem
import json
import hashlib
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')    


class Molecule:
    
    def __init__(self, cpd_smiles, cpd_name='', cpd_ref_chebi='', cpd_link='#', cpd_id='', cpd_num=1):
        
        self.cpd_id = cpd_id
        self.cpd_smiles = cpd_smiles
        self.cpd_name = cpd_name
        self.cpd_link = cpd_link
        self.cpd_ref_chebi = cpd_ref_chebi
        self.mol = Chem.MolFromSmiles(self.cpd_smiles) if cpd_smiles.strip() else None
        self.mol_num = cpd_num
        self.box_up_padding = 20
        

        if self.mol is None:
            raise ValueError(f"Could not parse SMILES string: {cpd_smiles}")
        
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
    
    
    def _get_stable_filename_hash(self):
        """Return a deterministic hash for filename based on CHEBI ID or canonical SMILES.

        Using the rendered SVG text causes filename changes when labels/styles vary.
        To keep a stable mapping for the same chemical structure, prefer CHEBI ID if
        present; otherwise use RDKit canonical SMILES.
        """
        try:
            if self.cpd_ref_chebi:
                key = self.cpd_ref_chebi.strip()
            else:
                # Canonical, isomeric SMILES for stability
                key = Chem.MolToSmiles(self.mol, isomericSmiles=True, canonical=True)
            return hashlib.md5(key.encode('utf-8')).hexdigest()
        except Exception:
            # Fallback to hashing the raw SMILES string
            return hashlib.md5(self.cpd_smiles.encode('utf-8')).hexdigest()

    def write_mol_svg(self, file_path):
        with open(file_path, 'w') as f:
            f.write(self.mol_svg)
    
    

    def get_drawer(self, mol_pic_width, mol_pic_height, tag_box_height):
        """Return an RDKit SVG drawer."""
        drawer = rdMolDraw2D.MolDraw2DSVG(width=mol_pic_width, 
                                    height=mol_pic_height + tag_box_height + 20 ,
                                    panelWidth=-1,
                                    panelHeight=mol_pic_height + tag_box_height
                                    )
        
        
        do = drawer.drawOptions()
        do.addStereoAnnotation = False
        do.explicitMethyl = True
        do.addStereoAnnotation = True


        do.annotationFontScale = 0.8
        do.baseFontSize = 0.6
        do.legendFontSize = 24


        do.setAtomPalette({
            6: (204/255, 153/255, 51/255, 1),
            8: (0/255, 153/255, 102/255, 1),
            7: (51/255, 102/255, 153/255, 1),
            15: (255/255, 153/255, 51/255, 1)
        })
        
        return drawer
    #endregion
    
    def add_mol_name_labele(self, mol_svg, mol_pic_width, mol_pic_height, mol_pic_size, tag_box_height):
        """Add the compound name and external link to the SVG."""
        box_width = int(mol_pic_width*0.8)
        box_height = 40

        
        x_box = int((mol_pic_width - box_width)/2)
        y_box = int(mol_pic_height + tag_box_height/2 + self.box_up_padding)
        
        x_text = int(x_box + box_width/2)
        y_text = int(y_box  +box_height/2 + 5)
        
        cpd_name_len = len(self.cpd_name)
        label_font_size = int(1.6*box_width/cpd_name_len)
        
        if label_font_size > 12:
            label_font_size = 12

        self.cpd_link = self.get_cpd_link()

        additional_svg = f'''
        <a href="{self.cpd_link}" target="_blank" style="cursor:pointer;">
        <g>
            <!-- label box -->
            <rect x="{x_box}" y="{y_box}" width="{box_width}" height="{box_height}" rx="10" ry="10" style="fill:#0066CC;" />
            <!-- compound name and external link -->
            <text x="{x_text }" y="{y_text}" font-size="{label_font_size}"  font-weight="bold" text-anchor="middle" fill="#FFFFFF">{self.cpd_name}</text>
        </g>
        </a>
        '''


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
        

        drawer = self.get_drawer(mol_pic_width=mol_pic_width, 
                                 mol_pic_height=mol_pic_height, 
                                 tag_box_height=tag_box_height
                                 
                                 )
        

        Chem.Draw.PrepareAndDrawMolecule(drawer, mol, kekulize=False)
        


        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
    
        
        if self.cpd_name != '':
            svg = self.add_mol_name_labele(mol_svg=svg,
                                        mol_pic_width=mol_pic_width,
                                        mol_pic_height=mol_pic_height,
                                        mol_pic_size=mol_pic_size, 
                                        tag_box_height= tag_box_height)
            

            

        if not os.path.exists(f'{cfg.DIR_CPD_SVG}'):
            os.makedirs(f'{cfg.DIR_CPD_SVG}')
        
        # Use a deterministic filename derived from CHEBI ID or canonical SMILES
        stable_hash = self._get_stable_filename_hash()
        file_name = f'{cfg.DIR_CPD_SVG}{stable_hash}.svg'
        
        # If already exists, skip re-writing to keep IO minimal
        if not os.path.exists(file_name):
            with open(f'{file_name}', "w", encoding="utf-8") as file:
                file.write(svg)
        
                
        return os.path.relpath(file_name, cfg.DIR_PROJECT_ROOT)
    
    def to_html(self):
        """Render this molecule as an HTML image block."""
        html_output = "<div style='display: flex; align-items: center;'>"
        html_output += f"<img src='{cfg.DIR_PROJECT_ROOT}/{self.mol_svg}' style='display:inline-block'/>"
        html_output += "</div>"

        return html_output
    
    def to_dict(self):
        return {
            'cpd_id'     : self.cpd_id,
            'cpd_smiles': self.cpd_smiles,
            'cpd_name': self.cpd_name,
            'cpd_ref_chebi': self.cpd_ref_chebi,
            'cpd_link':self.cpd_link,
            'mol_svg': self.mol_svg,
            'cpd_num': self.mol_num
        }
    

    def to_json(self):
        """Serialize this molecule to JSON."""
        return json.dumps(self.to_dict(), indent=4)
    
    
if __name__ == '__main__':
    smiles = 'C1=CC=C(C=C1)C(C)C(=O)O'
    name = 'ethanol'
    ref = 'CHEBI:63637'
    mol = Molecule(smiles, name, ref)
    print(mol.to_json())
