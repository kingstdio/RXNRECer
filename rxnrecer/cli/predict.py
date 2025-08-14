"""
Command-line interface for RXNRECer prediction
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Optional

from ..core.predictor import RXNRECerPredictor
from ..config.model_config import get_model_config
from ..config.settings import settings


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format=settings.LOG_FORMAT,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(settings.get_result_path('logs', 'rxnrecer.log'))
        ]
    )


def predict_command(args):
    """Main prediction command."""
    logger = logging.getLogger(__name__)
    
    try:
        # Validate input file
        input_path = Path(args.input)
        if not input_path.exists():
            logger.error(f"Input file not found: {args.input}")
            return 1
        
        # Validate output directory
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Get model configuration
        model_config = get_model_config(args.config)
        
        # Initialize predictor
        logger.info("Initializing RXNRECer predictor...")
        predictor = RXNRECerPredictor(
            model_config=model_config,
            model_path=args.model_path
        )
        
        # Perform prediction
        logger.info(f"Starting prediction for {input_path}")
        predictions = predictor.predict(
            input_data=str(input_path),
            batch_size=args.batch_size,
            top_k=args.top_k,
            include_equations=args.equations
        )
        
        # Save predictions
        logger.info(f"Saving predictions to {output_path}")
        predictor.save_predictions(
            predictions=predictions,
            output_file=str(output_path),
            output_format=args.format
        )
        
        # Print summary
        logger.info(f"Prediction completed successfully!")
        logger.info(f"Processed {len(predictions)} sequences")
        logger.info(f"Results saved to: {output_path}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Prediction failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def main():
    """Main entry point for prediction CLI."""
    parser = argparse.ArgumentParser(
        description='RXNRECer - Enzyme Reaction Prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic prediction
  rxnrecer predict input.fasta output.tsv
  
  # Prediction with custom parameters
  rxnrecer predict input.fasta output.json --batch-size 50 --top-k 10 --format json
  
  # Prediction with equations
  rxnrecer predict input.fasta output.tsv --equations
  
  # Verbose output
  rxnrecer predict input.fasta output.tsv --verbose
        """
    )
    
    # Required arguments
    parser.add_argument(
        'input',
        help='Input FASTA file with protein sequences'
    )
    parser.add_argument(
        'output',
        help='Output file path for predictions'
    )
    
    # Optional arguments
    parser.add_argument(
        '--model-path', '-m',
        help='Path to pre-trained model file'
    )
    parser.add_argument(
        '--config', '-c',
        choices=['default', 'small', 'large'],
        default='default',
        help='Model configuration preset (default: default)'
    )
    parser.add_argument(
        '--batch-size', '-b',
        type=int,
        help='Batch size for processing (default: from config)'
    )
    parser.add_argument(
        '--top-k', '-k',
        type=int,
        help='Number of top predictions to return (default: from config)'
    )
    parser.add_argument(
        '--format', '-f',
        choices=['tsv', 'csv', 'json'],
        default='tsv',
        help='Output format (default: tsv)'
    )
    parser.add_argument(
        '--equations', '-e',
        action='store_true',
        help='Include reaction equations in output'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Run prediction
    return predict_command(args)


if __name__ == '__main__':
    sys.exit(main())
