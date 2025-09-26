# 使用 setuptools_scm 的示例配置
# 这种方法会自动从 git 标签生成版本号

# pyproject.toml 配置示例:
"""
[build-system]
requires = ["setuptools>=45", "wheel", "setuptools_scm[toml]>=6.2"]
build-backend = "setuptools.build_meta"

[project]
name = "rxnrecer"
dynamic = ["version"]

[tool.setuptools_scm]
write_to = "rxnrecer/_version.py"
version_scheme = "python-simplified-semver"
local_scheme = "dirty-tag"

[tool.setuptools.dynamic]
version = {attr = "rxnrecer._version:version"}
"""

# 使用方式:
# 1. 创建 git 标签: git tag v1.3.0
# 2. 版本号自动从标签生成
# 3. 开发版本自动添加 .dev0+hash
