"""
Basic usage example for RXNRECer
"""

import rxnrecer
from rxnrecer.utils.file_utils import fasta_to_dataframe, save_dataframe
from rxnrecer.utils.bio_utils import validate_sequence, get_protein_properties

def main():
    """Demonstrate basic RXNRECer usage."""
    
    print("🚀 RXNRECer Basic Usage Example")
    print("=" * 50)
    
    # 1. Basic information
    print(f"Version: {rxnrecer.get_version()}")
    print(f"Info: {rxnrecer.get_info()}")
    print()
    
    # 2. Create sample FASTA data
    sample_fasta = """
>P12345
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
>P67890
MKLIVWALVLAFLACQGAVLGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEKIPVDGIKIVGDMVEV
"""
    
    # Write to temporary file
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(sample_fasta)
        temp_fasta = f.name
    
    try:
        # 3. Load and process sequences
        print("📊 Loading sequences...")
        df = fasta_to_dataframe(temp_fasta)
        print(f"Loaded {len(df)} sequences")
        print()
        
        # 4. Validate sequences
        print("🔍 Validating sequences...")
        valid_sequences = []
        for _, row in df.iterrows():
            if validate_sequence(row['seq']):
                valid_sequences.append(row)
                print(f"✓ {row['uniprot_id']}: Valid sequence")
            else:
                print(f"✗ {row['uniprot_id']}: Invalid sequence")
        print()
        
        # 5. Calculate protein properties
        print("🧬 Calculating protein properties...")
        for seq_data in valid_sequences:
            properties = get_protein_properties(seq_data['seq'])
            print(f"Protein: {seq_data['uniprot_id']}")
            print(f"  Length: {properties['length']}")
            print(f"  Molecular Weight: {properties['molecular_weight']:.2f} Da")
            print(f"  Isoelectric Point: {properties['isoelectric_point']:.2f}")
            print(f"  Hydrophobicity: {properties['hydrophobicity']:.3f}")
            print(f"  Charge: {properties['charge']:.2f}")
            print()
        
        # 6. Create mock predictions
        print("🔮 Creating mock predictions...")
        predictions = []
        for seq_data in valid_sequences:
            pred = {
                'uniprot_id': seq_data['uniprot_id'],
                'top_predictions': 'RXN-12345;RXN-67890;RXN-11111',
                'top_probabilities': '0.8500;0.1200;0.0300',
                'confidence': 0.85
            }
            predictions.append(pred)
        
        predictions_df = rxnrecer.pd.DataFrame(predictions)
        print(f"Generated predictions for {len(predictions)} proteins")
        print()
        
        # 7. Save results
        print("💾 Saving results...")
        output_file = "example_predictions.tsv"
        save_dataframe(predictions_df, output_file, 'tsv')
        print(f"Results saved to: {output_file}")
        print()
        
        # 8. Display results
        print("📈 Prediction Results:")
        print(predictions_df.to_string(index=False))
        print()
        
        print("✅ Example completed successfully!")
        
    finally:
        import os
        os.unlink(temp_fasta)

if __name__ == "__main__":
    main()
