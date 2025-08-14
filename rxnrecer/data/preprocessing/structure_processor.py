"""
Structure processing utilities for RXNRECer
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import logging
import subprocess
import tempfile
import os

from ...utils.file_utils import save_dataframe, load_dataframe


class StructureProcessor:
    """Structure processing and analysis utilities."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def process_pdb_structures(self,
                             protein_data: pd.DataFrame,
                             pdb_dir: str,
                             output_path: str) -> pd.DataFrame:
        """
        Process PDB structures for proteins.
        
        Args:
            protein_data: DataFrame with protein information
            pdb_dir: Directory containing PDB files
            output_path: Path to save processed structure data
            
        Returns:
            DataFrame with structure information
        """
        self.logger.info("Processing PDB structures...")
        
        structure_data = []
        
        for _, row in protein_data.iterrows():
            uniprot_id = row['uniprot_id']
            pdb_files = self._find_pdb_files(uniprot_id, pdb_dir)
            
            if pdb_files:
                # Get best PDB structure
                best_pdb = self._select_best_pdb(pdb_files)
                structure_info = self._extract_structure_info(best_pdb)
                
                structure_data.append({
                    'uniprot_id': uniprot_id,
                    'pdb_id': structure_info['pdb_id'],
                    'resolution': structure_info['resolution'],
                    'method': structure_info['method'],
                    'chain': structure_info['chain'],
                    'pdb_file': best_pdb
                })
            else:
                structure_data.append({
                    'uniprot_id': uniprot_id,
                    'pdb_id': None,
                    'resolution': None,
                    'method': None,
                    'chain': None,
                    'pdb_file': None
                })
        
        structure_df = pd.DataFrame(structure_data)
        save_dataframe(structure_df, output_path, 'feather')
        
        self.logger.info(f"Structure data processed: {len(structure_df)} proteins")
        return structure_df
    
    def _find_pdb_files(self, uniprot_id: str, pdb_dir: str) -> List[str]:
        """Find PDB files for a given UniProt ID."""
        pdb_files = []
        pdb_path = Path(pdb_dir)
        
        # Look for PDB files with UniProt ID in filename
        for pdb_file in pdb_path.glob(f"*{uniprot_id}*.pdb"):
            pdb_files.append(str(pdb_file))
        
        return pdb_files
    
    def _select_best_pdb(self, pdb_files: List[str]) -> str:
        """Select the best PDB structure based on resolution and method."""
        best_pdb = None
        best_resolution = float('inf')
        
        for pdb_file in pdb_files:
            resolution = self._extract_resolution(pdb_file)
            if resolution and resolution < best_resolution:
                best_resolution = resolution
                best_pdb = pdb_file
        
        return best_pdb or pdb_files[0] if pdb_files else None
    
    def _extract_resolution(self, pdb_file: str) -> Optional[float]:
        """Extract resolution from PDB file."""
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('REMARK   2 RESOLUTION.'):
                        parts = line.split()
                        if len(parts) >= 4:
                            return float(parts[3])
        except Exception as e:
            self.logger.warning(f"Could not extract resolution from {pdb_file}: {e}")
        
        return None
    
    def _extract_structure_info(self, pdb_file: str) -> Dict:
        """Extract structure information from PDB file."""
        info = {
            'pdb_id': None,
            'resolution': None,
            'method': None,
            'chain': None
        }
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('HEADER'):
                        parts = line.split()
                        if len(parts) >= 2:
                            info['pdb_id'] = parts[1]
                    elif line.startswith('REMARK   2 RESOLUTION.'):
                        parts = line.split()
                        if len(parts) >= 4:
                            info['resolution'] = float(parts[3])
                    elif line.startswith('EXPDTA'):
                        info['method'] = line[10:].strip()
                    elif line.startswith('ATOM') and not info['chain']:
                        info['chain'] = line[21]
        except Exception as e:
            self.logger.warning(f"Could not extract info from {pdb_file}: {e}")
        
        return info
    
    def generate_3di_embeddings(self,
                               structure_data: pd.DataFrame,
                               output_path: str,
                               foldseek_path: str = "foldseek") -> pd.DataFrame:
        """
        Generate 3DI embeddings using FoldSeek.
        
        Args:
            structure_data: DataFrame with structure information
            output_path: Path to save 3DI embeddings
            foldseek_path: Path to FoldSeek executable
            
        Returns:
            DataFrame with 3DI embeddings
        """
        self.logger.info("Generating 3DI embeddings...")
        
        # Filter proteins with structures
        proteins_with_structures = structure_data[structure_data['pdb_file'].notna()]
        
        embeddings_data = []
        
        for _, row in proteins_with_structures.iterrows():
            pdb_file = row['pdb_file']
            uniprot_id = row['uniprot_id']
            
            try:
                # Generate 3DI embedding using FoldSeek
                embedding = self._generate_single_3di_embedding(pdb_file, foldseek_path)
                
                embeddings_data.append({
                    'uniprot_id': uniprot_id,
                    'pdb_id': row['pdb_id'],
                    '3di_embedding': embedding
                })
                
            except Exception as e:
                self.logger.warning(f"Failed to generate 3DI embedding for {uniprot_id}: {e}")
        
        embeddings_df = pd.DataFrame(embeddings_data)
        save_dataframe(embeddings_df, output_path, 'feather')
        
        self.logger.info(f"3DI embeddings generated: {len(embeddings_df)} proteins")
        return embeddings_df
    
    def _generate_single_3di_embedding(self, pdb_file: str, foldseek_path: str) -> np.ndarray:
        """Generate 3DI embedding for a single PDB file."""
        # This is a placeholder implementation
        # In practice, you would use FoldSeek to generate actual 3DI embeddings
        
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as temp_pdb:
            temp_pdb_path = temp_pdb.name
        
        try:
            # Copy PDB file to temporary location
            import shutil
            shutil.copy2(pdb_file, temp_pdb_path)
            
            # Run FoldSeek (placeholder command)
            # cmd = [foldseek_path, "createdb", temp_pdb_path, "temp_db"]
            # subprocess.run(cmd, check=True)
            
            # For now, return a dummy embedding
            # In practice, you would extract the actual 3DI embedding
            embedding = np.random.rand(1024)  # Placeholder
            
            return embedding
            
        finally:
            os.unlink(temp_pdb_path)
    
    def calculate_structure_similarity(self,
                                     structure_data: pd.DataFrame,
                                     output_path: str) -> pd.DataFrame:
        """
        Calculate structure similarity between proteins.
        
        Args:
            structure_data: DataFrame with structure information
            output_path: Path to save similarity matrix
            
        Returns:
            DataFrame with similarity scores
        """
        self.logger.info("Calculating structure similarities...")
        
        # Filter proteins with structures
        proteins_with_structures = structure_data[structure_data['pdb_file'].notna()]
        
        similarity_data = []
        
        # Calculate pairwise similarities
        for i, row1 in proteins_with_structures.iterrows():
            for j, row2 in proteins_with_structures.iterrows():
                if i < j:  # Avoid duplicate calculations
                    try:
                        similarity = self._calculate_pairwise_similarity(
                            row1['pdb_file'], row2['pdb_file']
                        )
                        
                        similarity_data.append({
                            'protein1': row1['uniprot_id'],
                            'protein2': row2['uniprot_id'],
                            'similarity_score': similarity
                        })
                        
                    except Exception as e:
                        self.logger.warning(f"Failed to calculate similarity: {e}")
        
        similarity_df = pd.DataFrame(similarity_data)
        save_dataframe(similarity_df, output_path, 'feather')
        
        self.logger.info(f"Structure similarities calculated: {len(similarity_df)} pairs")
        return similarity_df
    
    def _calculate_pairwise_similarity(self, pdb1: str, pdb2: str) -> float:
        """Calculate similarity between two PDB structures."""
        # This is a placeholder implementation
        # In practice, you would use TM-align or similar tools
        
        # For now, return a random similarity score
        return np.random.random()
    
    def validate_structure_quality(self,
                                 structure_data: pd.DataFrame,
                                 min_resolution: float = 3.0) -> pd.DataFrame:
        """
        Validate structure quality.
        
        Args:
            structure_data: DataFrame with structure information
            min_resolution: Minimum acceptable resolution
            
        Returns:
            DataFrame with quality validation results
        """
        self.logger.info("Validating structure quality...")
        
        # Add quality flags
        structure_data = structure_data.copy()
        structure_data['has_structure'] = structure_data['pdb_file'].notna()
        structure_data['good_resolution'] = (
            structure_data['resolution'].notna() & 
            (structure_data['resolution'] <= min_resolution)
        )
        structure_data['xray_method'] = structure_data['method'].str.contains('X-RAY', case=False, na=False)
        
        # Overall quality score
        structure_data['quality_score'] = (
            structure_data['has_structure'].astype(int) * 0.3 +
            structure_data['good_resolution'].astype(int) * 0.4 +
            structure_data['xray_method'].astype(int) * 0.3
        )
        
        self.logger.info(f"Quality validation completed: {len(structure_data)} structures")
        return structure_data
