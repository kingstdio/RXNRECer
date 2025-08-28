# RXNRECer Makefile
# Simple commands for building and managing the package

.PHONY: help clean build test install uninstall dist clean-all

help:  ## Show this help message
	@echo "RXNRECer Package Management"
	@echo "=========================="
	@echo ""
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

clean:  ## Clean build artifacts
	@echo "🧹 Cleaning build artifacts..."
	@rm -rf build/ dist/ *.egg-info/ __pycache__/ .pytest_cache/ .mypy_cache/ .ruff_cache/
	@find . -name "*.pyc" -delete 2>/dev/null || true
	@find . -name "*.pyo" -delete 2>/dev/null || true
	@find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	@echo "✅ Cleaned"

build: clean  ## Build the package
	@echo "🏗️  Building package..."
	@python -m build
	@echo "✅ Built successfully"

test: build  ## Test the built package
	@echo "🧪 Testing package..."
	@python build/build.py
	@echo "✅ Tests passed"

install: build  ## Install the package in development mode
	@echo "📦 Installing package..."
	@pip install -e .
	@echo "✅ Installed successfully"

uninstall:  ## Uninstall the package
	@echo "🗑️  Uninstalling package..."
	@pip uninstall rxnrecer -y
	@echo "✅ Uninstalled successfully"

dist: build  ## Create distribution packages
	@echo "📦 Creating distribution packages..."
	@mkdir -p build/dist
	@mv dist/* build/dist/ 2>/dev/null || true
	@mv *.egg-info build/ 2>/dev/null || true
	@echo "✅ Distribution packages created in build/dist/"

clean-all: clean  ## Clean everything including build directory
	@echo "🧹 Cleaning everything..."
	@rm -rf build/
	@echo "✅ Everything cleaned"

version:  ## Show current version
	@python -c "import rxnrecer; print(f'RXNRECer version: {rxnrecer.__version__}')"

check:  ## Check package structure
	@echo "🔍 Checking package structure..."
	@python -c "import rxnrecer; print('✅ Package imports successfully')"
	@rxnrecer --help > /dev/null && echo "✅ CLI works" || echo "❌ CLI failed"
	@rxnrecer-download-data --help > /dev/null && echo "✅ Download command works" || echo "❌ Download command failed"
	@rxnrecer-cache --help > /dev/null && echo "✅ Cache command works" || echo "❌ Cache command failed"
