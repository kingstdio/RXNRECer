"""
Basic tests for RXNRECer
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os

# Import RXNRECer modules
from rxnrecer.utils.file_utils import (
    fasta_to_dataframe, dataframe_to_fasta, save_dataframe, load_dataframe
)
from rxnrecer.utils.bio_utils import (
    validate_sequence, clean_sequence, get_protein_properties
)
from rxnrecer.config.settings import settings
from rxnrecer.config.model_config import ModelConfig, get_model_config


class TestFileUtils:
    """Test file utility functions."""
    
    def test_fasta_to_dataframe(self):
        """Test FASTA to DataFrame conversion."""
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">P12345\nMKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG\n")
            f.write(">P67890\nMKLIVWALVLAFLACQGAVLGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEKIPVDGIKIVGDMVEV\n")
            temp_fasta = f.name
        
        try:
            # Convert to DataFrame
            df = fasta_to_dataframe(temp_fasta)
            
            # Check structure
            assert len(df) == 2
            assert list(df.columns) == ['uniprot_id', 'seq']
            assert df.iloc[0]['uniprot_id'] == 'P12345'
            assert df.iloc[1]['uniprot_id'] == 'P67890'
            
        finally:
            os.unlink(temp_fasta)
    
    def test_dataframe_to_fasta(self):
        """Test DataFrame to FASTA conversion."""
        # Create test DataFrame
        df = pd.DataFrame({
            'uniprot_id': ['P12345', 'P67890'],
            'seq': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG',
                   'MKLIVWALVLAFLACQGAVLGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEKIPVDGIKIVGDMVEV']
        })
        
        # Convert to FASTA
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_fasta = f.name
        
        try:
            dataframe_to_fasta(df, temp_fasta)
            
            # Read back and verify
            with open(temp_fasta, 'r') as f:
                content = f.read()
            
            assert '>P12345' in content
            assert '>P67890' in content
            assert 'MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG' in content
            
        finally:
            os.unlink(temp_fasta)
    
    def test_save_load_dataframe(self):
        """Test DataFrame save and load functions."""
        # Create test DataFrame
        df = pd.DataFrame({
            'id': [1, 2, 3],
            'name': ['A', 'B', 'C'],
            'value': [1.1, 2.2, 3.3]
        })
        
        # Test TSV format
        with tempfile.NamedTemporaryFile(suffix='.tsv', delete=False) as f:
            temp_tsv = f.name
        
        try:
            save_dataframe(df, temp_tsv, 'tsv')
            loaded_df = load_dataframe(temp_tsv)
            pd.testing.assert_frame_equal(df, loaded_df)
        finally:
            os.unlink(temp_tsv)
        
        # Test CSV format
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
            temp_csv = f.name
        
        try:
            save_dataframe(df, temp_csv, 'csv')
            loaded_df = load_dataframe(temp_csv)
            pd.testing.assert_frame_equal(df, loaded_df)
        finally:
            os.unlink(temp_csv)


class TestBioUtils:
    """Test biological utility functions."""
    
    def test_validate_sequence(self):
        """Test sequence validation."""
        # Valid sequences
        assert validate_sequence("MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG")
        assert validate_sequence("ACDEFGHIKLMNPQRSTVWY")
        
        # Invalid sequences
        assert not validate_sequence("MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGGX")
        assert not validate_sequence("123456789")
        assert not validate_sequence("")
    
    def test_clean_sequence(self):
        """Test sequence cleaning."""
        # Test with valid characters
        seq = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        cleaned = clean_sequence(seq)
        assert cleaned == seq
        
        # Test with invalid characters
        seq_with_invalid = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGGX123"
        cleaned = clean_sequence(seq_with_invalid)
        assert cleaned == "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        
        # Test case conversion
        seq_lower = "mktvrqerlksivrilerskepvsgaqlaeelsvsrqvivqdiyalrslgynivatprgyvlagg"
        cleaned = clean_sequence(seq_lower)
        assert cleaned == "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
    
    def test_get_protein_properties(self):
        """Test protein property calculation."""
        seq = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        properties = get_protein_properties(seq)
        
        # Check that all properties are present
        expected_keys = ['length', 'molecular_weight', 'isoelectric_point', 'hydrophobicity', 'charge']
        for key in expected_keys:
            assert key in properties
        
        # Check specific values
        assert properties['length'] == len(seq)
        assert properties['molecular_weight'] > 0
        assert 0 <= properties['isoelectric_point'] <= 14
        assert isinstance(properties['hydrophobicity'], float)
        assert isinstance(properties['charge'], float)


class TestConfig:
    """Test configuration modules."""
    
    def test_settings(self):
        """Test settings configuration."""
        # Check that settings object exists
        assert hasattr(settings, 'PROJECT_ROOT')
        assert hasattr(settings, 'DATA_ROOT')
        assert hasattr(settings, 'RESULTS_ROOT')
        
        # Check that directories are created
        assert settings.RAW_DATA_DIR.exists()
        assert settings.PROCESSED_DATA_DIR.exists()
        assert settings.MODELS_DIR.exists()
    
    def test_model_config(self):
        """Test model configuration."""
        # Test default config
        config = ModelConfig()
        assert config.esm_out_dim == 1280
        assert config.gru_h_dim == 512
        assert config.att_dim == 32
        assert config.dropout_rate == 0.2
        
        # Test custom config
        custom_config = ModelConfig(
            esm_out_dim=640,
            gru_h_dim=256,
            att_dim=16,
            dropout_rate=0.1
        )
        assert custom_config.esm_out_dim == 640
        assert custom_config.gru_h_dim == 256
        
        # Test config presets
        small_config = get_model_config('small')
        assert small_config.gru_h_dim == 256
        assert small_config.att_dim == 16
        
        large_config = get_model_config('large')
        assert large_config.gru_h_dim == 1024
        assert large_config.att_dim == 64
    
    def test_config_serialization(self):
        """Test config serialization."""
        config = ModelConfig()
        
        # Test to_dict
        config_dict = config.to_dict()
        assert isinstance(config_dict, dict)
        assert 'esm_out_dim' in config_dict
        assert 'gru_h_dim' in config_dict
        
        # Test from_dict
        new_config = ModelConfig.from_dict(config_dict)
        assert new_config.esm_out_dim == config.esm_out_dim
        assert new_config.gru_h_dim == config.gru_h_dim


class TestIntegration:
    """Test integration scenarios."""
    
    def test_basic_workflow(self):
        """Test basic workflow with sample data."""
        # Create sample FASTA data
        sample_fasta = """
>P12345
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
>P67890
MKLIVWALVLAFLACQGAVLGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEKIPVDGIKIVGDMVEV
"""
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(sample_fasta)
            temp_fasta = f.name
        
        try:
            # Test complete workflow
            df = fasta_to_dataframe(temp_fasta)
            assert len(df) == 2
            
            # Validate sequences
            valid_sequences = []
            for _, row in df.iterrows():
                if validate_sequence(row['seq']):
                    valid_sequences.append(row)
            
            assert len(valid_sequences) == 2
            
            # Test saving results
            results_df = pd.DataFrame({
                'uniprot_id': [row['uniprot_id'] for row in valid_sequences],
                'prediction': ['RXN-12345', 'RXN-67890'],
                'confidence': [0.85, 0.92]
            })
            
            with tempfile.NamedTemporaryFile(suffix='.tsv', delete=False) as f:
                temp_results = f.name
            
            try:
                save_dataframe(results_df, temp_results, 'tsv')
                loaded_results = load_dataframe(temp_results)
                pd.testing.assert_frame_equal(results_df, loaded_results)
            finally:
                os.unlink(temp_results)
                
        finally:
            os.unlink(temp_fasta)


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
