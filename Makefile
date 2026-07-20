# RXNRECer Makefile
# Provides convenient commands for version management and development

.PHONY: help version current bump-major bump-minor bump-patch set-version install test clean build release

PYTHON ?= /hpcfs/fhome/shizhenkun/miniconda3/envs/rxnrecer-release/bin/python
PIP ?= /hpcfs/fhome/shizhenkun/miniconda3/envs/rxnrecer-release/bin/pip

# Default target
help:
	@echo "RXNRECer Development Commands"
	@echo "=============================="
	@echo ""
	@echo "Version Management:"
	@echo "  version          Show current version"
	@echo "  bump-major       Bump major version (1.0.0 -> 2.0.0)"
	@echo "  bump-minor       Bump minor version (1.0.0 -> 1.1.0)"
	@echo "  bump-patch       Bump patch version (1.0.0 -> 1.0.1)"
	@echo "  set-version VER  Set specific version (e.g., make set-version VER=1.4.0)"
	@echo "  update-docs      Update all documentation version references"
	@echo ""
	@echo "Development:"
	@echo "  install          Install package in development mode"
	@echo "  test             Run tests"
	@echo "  clean            Clean build artifacts"
	@echo "  build            Build package"
	@echo "  release          Build and release to PyPI"
	@echo ""

# Version management
version current:
	@$(PYTHON) scripts/version_manager.py current

bump-major:
	@$(PYTHON) scripts/version_manager.py bump major

bump-minor:
	@$(PYTHON) scripts/version_manager.py bump minor

bump-patch:
	@$(PYTHON) scripts/version_manager.py bump patch

set-version:
	@if [ -z "$(VER)" ]; then \
		echo "❌ Error: Please specify version with VER=1.4.0"; \
		exit 1; \
	fi
	@$(PYTHON) scripts/version_manager.py set $(VER)

update-docs:
	@echo "📚 Updating documentation version references..."
	@$(PYTHON) scripts/update_docs.py

# Development commands
install:
	@echo "📦 Installing RXNRECer in development mode..."
	$(PIP) install -e .

test:
	@echo "🧪 Running tests..."
	$(PYTHON) -m unittest discover -s tests -v

clean:
	@echo "🧹 Cleaning build artifacts..."
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf __pycache__/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} +

build:
	@echo "🔨 Building package..."
	$(PYTHON) -m build

release: clean build
	@echo "🚀 Releasing to PyPI..."
	$(PYTHON) -m twine upload dist/*

# Git integration
commit-version:
	@echo "📝 Committing version changes..."
	git add .
	git commit -m "Bump version to $(shell $(PYTHON) scripts/version_manager.py current)"

tag-version:
	@echo "🏷️  Creating version tag..."
	git tag v$(shell $(PYTHON) scripts/version_manager.py current)

# Version synchronization
sync-versions:
	@echo "🔄 Synchronizing all version references..."
	@python scripts/sync_versions.py

# Complete version bump workflow
bump-and-commit: bump-patch sync-versions commit-version tag-version
	@echo "✅ Version bumped, synchronized, committed, and tagged!"
	@echo "💡 Run 'git push origin release --tags' to push changes"

# Pre-commit workflow
pre-commit: sync-versions
	@echo "✅ Pre-commit version sync completed"

# Build workflow with version sync
build-with-sync: sync-versions clean build
	@echo "✅ Build completed with version synchronization"

# Release workflow with version sync
release-with-sync: sync-versions clean build release
	@echo "✅ Release completed with version synchronization"
