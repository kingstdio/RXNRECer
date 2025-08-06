# RXNRECer - Cursor 设置指南

## 环境设置

### 1. Conda 环境
项目使用 conda 环境 `rxnrecer`，已经配置好所有必要的依赖包：

```bash
# 激活环境
conda activate rxnrecer

# 检查环境
python --version  # Python 3.10.13
pip list | grep torch  # PyTorch 2.5.1+cu121
```

### 2. Cursor 配置
项目已配置 `.vscode/settings.json`，自动设置：
- Python 解释器路径
- 代码分析路径
- 文件排除规则
- 代码格式化设置

## 快速开始

### 1. 测试项目状态
```bash
python test_basic.py
```

### 2. 运行预测
```bash
# 基本预测
python run_prediction.py data/sample/sample10.fasta results/prediction.json

# 使用集成预测
python run_prediction.py data/sample/sample10.fasta results/prediction.json --ensemble

# 包含反应方程式
python run_prediction.py data/sample/sample10.fasta results/prediction.json --equations

# 自定义批次大小
python run_prediction.py data/sample/sample10.fasta results/prediction.json --batch-size 50
```

### 3. 查看帮助
```bash
python run_prediction.py --help
```

## 项目结构

```
RXNRECer/
├── rxnrecer.py              # 主入口文件
├── run_prediction.py        # 快速预测脚本
├── test_basic.py           # 基础测试脚本
├── example_usage.py        # 使用示例
├── config/                 # 配置文件
├── modules/                # 核心模块
├── methods/                # 方法实现
├── tools/                  # 工具库
├── data/                   # 数据目录
└── results/                # 结果目录
```

## 主要功能

### 1. 反应预测
- 基于蛋白质序列预测酶催化的反应
- 支持多种嵌入模型（ESM、T5、UniRep）
- 集成预测提高准确性

### 2. 模型特性
- 使用 ESM-2 模型进行蛋白质序列嵌入
- BGRU + 注意力机制
- 支持 GPU 加速

### 3. 输出格式
- JSON：结构化数据
- TSV：表格格式
- CSV：逗号分隔格式

## 开发建议

### 1. 代码组织
- 使用模块化设计
- 保持函数职责单一
- 添加适当的错误处理

### 2. 测试
- 运行 `test_basic.py` 确保环境正常
- 使用小数据集进行测试
- 验证输出格式正确性

### 3. 性能优化
- 调整批次大小
- 使用 GPU 加速
- 考虑内存使用

## 常见问题

### 1. 环境问题
```bash
# 如果遇到导入错误
conda activate rxnrecer
python -c "import torch; print(torch.cuda.is_available())"
```

### 2. 路径问题
```bash
# 检查配置文件路径
python -c "from config import conf as cfg; print(cfg.DIR_PROJECT_ROOT)"
```

### 3. 模型加载问题
```bash
# 检查模型文件
ls -la data/model/production_185846best.pth
```

## 下一步

1. 运行基础测试确保环境正常
2. 尝试小规模预测
3. 根据需要调整参数
4. 开始重构项目结构

## 联系信息

如有问题，请检查：
1. Conda 环境是否正确激活
2. 依赖包是否完整安装
3. 数据文件是否存在
4. GPU 是否可用 