#!/usr/bin/env python3
"""
Example usage of RXNRECer for reaction prediction
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from rxnrecer import step_by_step_prediction

def main():
    """Example of using RXNRECer for prediction"""
    
    # Input and output files
    input_file = "data/sample/sample10.fasta"
    output_file = "results/sample_prediction.json"
    
    # Create results directory if it doesn't exist
    os.makedirs("results", exist_ok=True)
    
    print("RXNRECer Prediction Example")
    print("=" * 40)
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print()
    
    try:
        # Run prediction
        print("Running prediction...")
        step_by_step_prediction(
            input_data=input_file,
            output_file=output_file,
            output_format='json',
            getEquation=True,
            Ensemble=True,
            batch_size=5  # Small batch size for demo
        )
        
        print(f"✓ Prediction completed! Results saved to {output_file}")
        
        # Show first few results
        import json
        with open(output_file, 'r') as f:
            results = json.load(f)
        
        print(f"\nPredicted {len(results)} reactions:")
        for i, result in enumerate(results[:3]):  # Show first 3
            print(f"{i+1}. {result.get('uniprot_id', 'N/A')}: {result.get('predicted_reaction', 'N/A')}")
        
    except Exception as e:
        print(f"❌ Prediction failed: {e}")
        return False
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 