# RXNRECer 项目重构计划

## 1. 项目概述
RXNRECer 是一个基于深度学习的酶反应预测系统，通过蛋白质序列预测酶催化的化学反应。

## 2. 重构后的项目结构

```
rxnrecer/
├── README.md                    # 项目说明文档
├── setup.py                     # 安装脚本
├── requirements.txt             # 依赖包列表
├── requirements-dev.txt         # 开发环境依赖
├── pyproject.toml              # 项目配置文件
├── .gitignore                  # Git忽略文件
├── LICENSE                     # 许可证文件
│
├── rxnrecer/                   # 主包目录
│   ├── __init__.py
│   ├── core/                   # 核心功能模块
│   │   ├── __init__.py
│   │   ├── predictor.py        # 预测器主类
│   │   ├── model.py           # 模型定义
│   │   ├── data_loader.py     # 数据加载器
│   │   └── trainer.py         # 训练器
│   │
│   ├── models/                 # 模型模块
│   │   ├── __init__.py
│   │   ├── esm_embedding.py   # ESM嵌入模型
│   │   ├── bgru_model.py      # BGRU模型
│   │   ├── attention.py       # 注意力机制
│   │   └── ensemble.py        # 集成模型
│   │
│   ├── data/                   # 数据处理模块
│   │   ├── __init__.py
│   │   ├── fasta_parser.py    # FASTA文件解析
│   │   ├── dataset.py         # 数据集类
│   │   ├── preprocessor.py    # 数据预处理
│   │   └── augmentation.py    # 数据增强
│   │
│   ├── utils/                  # 工具模块
│   │   ├── __init__.py
│   │   ├── bio_utils.py       # 生物信息学工具
│   │   ├── file_utils.py      # 文件操作工具
│   │   ├── sequence_utils.py  # 序列处理工具
│   │   ├── structure_utils.py # 结构分析工具
│   │   └── evaluation_utils.py # 评估工具
│   │
│   ├── config/                 # 配置模块
│   │   ├── __init__.py
│   │   ├── settings.py        # 基础配置
│   │   ├── model_config.py    # 模型配置
│   │   └── paths.py           # 路径配置
│   │
│   └── cli/                    # 命令行接口
│       ├── __init__.py
│       ├── predict.py         # 预测命令
│       ├── train.py           # 训练命令
│       └── evaluate.py        # 评估命令
│
├── scripts/                    # 脚本目录
│   ├── setup_environment.sh   # 环境设置脚本
│   ├── download_data.py       # 数据下载脚本
│   ├── preprocess_data.py     # 数据预处理脚本
│   └── benchmark.py           # 基准测试脚本
│
├── tests/                      # 测试目录
│   ├── __init__.py
│   ├── test_predictor.py
│   ├── test_models.py
│   ├── test_data.py
│   └── test_utils.py
│
├── docs/                       # 文档目录
│   ├── api.md                 # API文档
│   ├── user_guide.md          # 用户指南
│   ├── developer_guide.md     # 开发者指南
│   └── examples/              # 示例代码
│
├── data/                       # 数据目录
│   ├── raw/                   # 原始数据
│   ├── processed/             # 处理后数据
│   ├── models/                # 预训练模型
│   └── samples/               # 示例数据
│
├── results/                    # 结果目录
│   ├── predictions/           # 预测结果
│   ├── evaluations/           # 评估结果
│   └── logs/                  # 日志文件
│
└── notebooks/                  # Jupyter笔记本
    ├── exploration.ipynb      # 数据探索
    ├── model_analysis.ipynb   # 模型分析
    └── results_analysis.ipynb # 结果分析
```

## 3. 重构步骤

### 第一阶段：项目结构重组
1. 创建新的目录结构
2. 移动和重命名文件
3. 创建必要的__init__.py文件
4. 更新导入路径

### 第二阶段：代码重构
1. 提取公共方法到utils模块
2. 重构核心预测逻辑
3. 统一配置管理
4. 优化模型架构

### 第三阶段：环境配置
1. 创建虚拟环境
2. 更新依赖管理
3. 创建安装脚本
4. 配置开发环境

### 第四阶段：测试和文档
1. 编写单元测试
2. 创建API文档
3. 编写用户指南
4. 创建示例代码

## 4. 重复代码提取计划

### 4.1 文件操作相关
- `fasta_to_dataframe()` - 在多个文件中重复
- `table2fasta()` - 在bioFunctionLib.py中
- JSON文件读写操作

### 4.2 数据处理相关
- 数据集分割方法（ta_ds_spliter, tb_ds_spliter等）
- 评估指标计算
- 序列处理工具

### 4.3 模型相关
- ESM嵌入提取
- 注意力机制实现
- 模型加载和保存

### 4.4 生物信息学工具
- BLAST比对
- 结构分析
- 序列比对

## 5. 依赖管理

### 5.1 核心依赖
- torch >= 2.0.0
- fair-esm >= 2.0.0
- transformers >= 4.30.0
- pandas >= 1.5.0
- numpy >= 1.23.0
- biopython >= 1.80

### 5.2 开发依赖
- pytest >= 7.0.0
- black >= 22.0.0
- flake8 >= 5.0.0
- mypy >= 1.0.0
- jupyter >= 1.0.0

## 6. 配置管理

### 6.1 环境配置
- 开发/生产环境切换
- GPU/CPU设备配置
- 路径配置管理

### 6.2 模型配置
- 模型参数配置
- 训练参数配置
- 预测参数配置

## 7. 测试策略

### 7.1 单元测试
- 核心功能测试
- 工具函数测试
- 模型组件测试

### 7.2 集成测试
- 端到端预测测试
- 数据流程测试
- 性能基准测试

## 8. 文档计划

### 8.1 用户文档
- 安装指南
- 快速开始
- API参考
- 故障排除

### 8.2 开发者文档
- 架构设计
- 贡献指南
- 代码规范
- 发布流程
