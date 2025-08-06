#!/usr/bin/env python3
"""
Quick prediction script for RXNRECer
Usage: python run_prediction.py <input_fasta> <output_file>
"""

import sys
import os
import argparse

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from rxnrecer import step_by_step_prediction

def main():
    parser = argparse.ArgumentParser(description='Run RXNRECer prediction')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('output', help='Output file path')
    parser.add_argument('--format', '-f', default='json', 
                       choices=['json', 'tsv', 'csv'], 
                       help='Output format (default: json)')
    parser.add_argument('--batch-size', '-b', type=int, default=100,
                       help='Batch size for processing (default: 100)')
    parser.add_argument('--ensemble', '-e', action='store_true',
                       help='Use ensemble prediction')
    parser.add_argument('--equations', action='store_true',
                       help='Include reaction equations')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"❌ Input file not found: {args.input}")
        return 1
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"RXNRECer Prediction")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Format: {args.format}")
    print(f"Batch size: {args.batch_size}")
    print(f"Ensemble: {args.ensemble}")
    print(f"Equations: {args.equations}")
    print("-" * 50)
    
    try:
        step_by_step_prediction(
            input_data=args.input,
            output_file=args.output,
            output_format=args.format,
            getEquation=args.equations,
            Ensemble=args.ensemble,
            batch_size=args.batch_size
        )
        
        print(f"✅ Prediction completed successfully!")
        print(f"Results saved to: {args.output}")
        return 0
        
    except Exception as e:
        print(f"❌ Prediction failed: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 