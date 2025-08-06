#!/usr/bin/env python3
"""
Basic test script for RXNRECer project
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

def test_imports():
    """Test if all essential modules can be imported"""
    print("Testing imports...")
    
    try:
        from config import conf as cfg
        print("✓ Config module imported")
    except Exception as e:
        print(f"✗ Config module import failed: {e}")
        return False
    
    try:
        from methods import Mactive
        print("✓ Mactive module imported")
    except Exception as e:
        print(f"✗ Mactive module import failed: {e}")
        return False
    
    try:
        from modules import commonfunction
        print("✓ Common function module imported")
    except Exception as e:
        print(f"✗ Common function module import failed: {e}")
        return False
    
    try:
        from modules.predict import predRXN
        print("✓ Prediction module imported")
    except Exception as e:
        print(f"✗ Prediction module import failed: {e}")
        return False
    
    return True

def test_config():
    """Test configuration settings"""
    print("\nTesting configuration...")
    
    from config import conf as cfg
    
    # Check essential paths
    essential_paths = [
        cfg.DIR_PROJECT_ROOT,
        cfg.DATA_ROOT,
        cfg.FILE_WEIGHT_PRODUCTION_BEST_MODEL,
        cfg.FILE_DS_DICT_ID2RXN
    ]
    
    for path in essential_paths:
        if os.path.exists(path):
            print(f"✓ Path exists: {path}")
        else:
            print(f"✗ Path missing: {path}")
            return False
    
    return True

def test_model_loading():
    """Test if the model can be loaded"""
    print("\nTesting model loading...")
    
    try:
        from rxnrecer import load_model
        model = load_model()
        print("✓ Model loaded successfully")
        return True
    except Exception as e:
        print(f"✗ Model loading failed: {e}")
        return False

def test_fasta_processing():
    """Test FASTA file processing"""
    print("\nTesting FASTA processing...")
    
    try:
        from rxnrecer import fasta_to_dataframe
        import pandas as pd
        
        # Test with sample data
        sample_fasta = "data/sample/sample10.fasta"
        if os.path.exists(sample_fasta):
            df = fasta_to_dataframe(sample_fasta)
            print(f"✓ FASTA processing successful, loaded {len(df)} sequences")
            return True
        else:
            print(f"✗ Sample FASTA file not found: {sample_fasta}")
            return False
    except Exception as e:
        print(f"✗ FASTA processing failed: {e}")
        return False

def main():
    """Run all tests"""
    print("RXNRECer Basic Test Suite")
    print("=" * 50)
    
    tests = [
        test_imports,
        test_config,
        test_fasta_processing,
        test_model_loading
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
        print()
    
    print("=" * 50)
    print(f"Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! The project is ready to use.")
        return True
    else:
        print("❌ Some tests failed. Please check the issues above.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 