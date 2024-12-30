import os
import sys
import numpy as np

# 动态添加路径
current_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, current_dir)
sys.path.insert(1, os.path.join(current_dir, '../../'))

# 从配置文件导入配置
from config import conf as cfg


class Tdi:
    def __init__(self, name: str = '', seq: str = '', token_3di: str = '', matrix_3di: np.ndarray = None):
        """
        初始化 Tdi 类的属性
        :param name: 名称
        :param seq: 序列
        :param token_3di: 氨基酸特性标记
        :param matrix_3di: 3Di 嵌入矩阵
        """
        self.name = name
        self.seq = seq
        self.token_3di = token_3di
        self.matrix_3di = matrix_3di if matrix_3di is not None else np.array([])

    def read_3di_file(self, path: str) -> None:
        """
        从 3di 文件中读取数据并解析为类属性
        :param path: 文件路径
        """
        if not os.path.exists(path):
            print(f"文件不存在: {path}")
            return

        try:
            with open(path, 'r', encoding='utf-8') as file:
                content = file.read()

            # 去除换行符并分割字段
            content = content.strip().split('\t')

            # 检查字段数量是否足够
            if len(content) < 4:
                raise IndexError(f"文件内容少于预期的 4 个字段:{content}")
                
            # 解析字段
            self.name = content[0]
            self.seq = content[1]
            self.token_3di = content[2]
            
            # 将最后一个字段解析为浮点数数组，并重塑为矩阵
            matrix_3di = np.array(content[3].split(','), dtype=float)
            self.matrix_3di = matrix_3di.reshape(-1, 10)

        except FileNotFoundError:
            print(f"文件未找到: {path}")
        except ValueError as e:
            print(f"数据解析错误: {e}")
        except IndexError as e:
            print(f"文件格式错误: {e}")

    def show(self) -> None:
        """
        显示当前对象的所有属性
        """
        print(f'name: {self.name}')
        print(f'seq: {self.seq}')
        print(f'token_3di: {self.token_3di}')
        print(f'matrix_3di:\n{self.matrix_3di}')

    def get_name(self) -> str:
        """
        返回对象的名称
        :return: 名称字符串
        """
        return self.name

    def get_seq(self) -> str:
        """
        返回对象的序列
        :return: 序列字符串
        """
        return self.seq

    def get_token_3di(self) -> str:
        """
        返回对象的氨基酸特性标记
        :return: 氨基酸特性字符串
        """
        return self.token_3di

    def get_matrix_3di(self) -> np.ndarray:
        """
        返回对象的 3Di 嵌入矩阵
        :return: NumPy 数组
        """
        return self.matrix_3di
