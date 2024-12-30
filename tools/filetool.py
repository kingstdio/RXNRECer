'''
Author: Zhenkun Shi
Date: 2023-06-21 14:43:20
LastEditors: Zhenkun Shi
LastEditTime: 2023-09-22 18:36:56
FilePath: /preaction/tools/filetool.py
Description: 

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
from retrying import retry
from matplotlib.pyplot import axis
import urllib3
import pandas as pd

from shutil import copyfile
from sys import exit
import zipfile
import tarfile
import shutil
import subprocess


#region 下载文件
@retry(stop_max_attempt_number=10, wait_random_min=10, wait_random_max=20)
def download(download_url, save_file):
    """[从网址下载文件]
    Args:
        download_url ([Stirng]): [下载路径]
        save_file ([String]): [保存文件路径]
    """
    http = urllib3.PoolManager()
    response = http.request('GET', download_url)
    with open(save_file, 'wb') as f:
        f.write(response.data)
    response.release_conn()
#endregion


def wget(download_url, save_file, verbos=False):
    process = os.popen('which wget') # return file
    output = process.read()
    if output =='':
        print('wget not installed')
    else:
        if verbos:
            cmd = 'wget ' + download_url + ' -O ' + save_file
        else:
            cmd = 'wget -q ' + download_url + ' -O ' + save_file
        print (cmd)
        process = os.popen(cmd)
        output = process.read()
    process.close()
    
def download_with_wget(url, dst):
    """
    Download a file using wget and handle errors.

    Parameters:
    url (str): The URL of the file to download.
    dst (str): The destination path to save the file.

    Returns:
    bool: True if download is successful, False if not.
    """
    try:
        # 调用 wget 命令并捕获输出
        result = subprocess.run(
            ["wget", "-O", dst, url],
            stdout=subprocess.PIPE,  # 捕获标准输出
            stderr=subprocess.PIPE,  # 捕获错误输出
            text=True  # 输出为文本而非字节
        )
        
        # 检查返回码和输出内容
        if result.returncode == 0:
            print(f"File downloaded successfully to {dst}")
            return True
        else:
            return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

def convert_DF_dateTime(inputdf):
    """[Covert unisprot csv records datatime]

    Args:
        inputdf ([DataFrame]): [input dataFrame]

    Returns:
        [DataFrame]: [converted DataFrame]
    """
    inputdf.date_integraged = pd.to_datetime(inputdf['date_integraged'])
    inputdf.date_sequence_update = pd.to_datetime(inputdf['date_sequence_update'])
    inputdf.date_annotation_update = pd.to_datetime(inputdf['date_annotation_update'])
    inputdf = inputdf.sort_values(by='date_integraged', ascending=True)
    inputdf.reset_index(drop=True, inplace=True)
    return inputdf

def get_file_names_in_dir(dataroot, filetype='all'):
    """返回某个文件夹下的文件列表，可以指定文件类型

    Args:
        dataroot (string): 文件夹目录
        filetype (str, optional): 文件类型. Defaults to ''.

    Returns:
        DataFrame: columns=['filename','filetype','filename_no_suffix']
    """
    exist_file_df = pd.DataFrame(os.listdir(dataroot), columns=['filename'])

    if len(exist_file_df)!=0:
        exist_file_df['filetype'] = exist_file_df.filename.apply(lambda x: x.split('.')[-1])
        exist_file_df['filename_no_suffix'] = exist_file_df.apply(lambda x: x['filename'].replace(('.'+str(x['filetype']).strip()), ''), axis=1)
        if filetype !='all':
            exist_file_df = exist_file_df[exist_file_df.filetype==filetype]
        return exist_file_df
    else:
        return pd.DataFrame(columns=['filename','filetype','filename_no_suffix'])


def copy(source, target):
    """ 拷贝文件

    Args:
        source (string): source
        target (string): target
    """
    try:
        copyfile(src=source, dst=target)
    except IOError as e:
        print("Unable to copy file. %s" % e)
        exit(1)
    except:
        print("Unexpected error:", sys.exc_info())
        exit(1)

def delete(filepath):
    """删除文件

    Args:
        filepath (string): 文件全路径
    """
    try:
        os.remove(filepath)
    except IOError as e:
        print("Unable to delete file. %s" % e)
        exit(1)
    except:
        print("Unexpected error:", sys.exc_info())
        exit(1)

#region unzip file
def unzipfile(filename, target_dir):
    """uzip file
    Args:
        zipfile (string): zip file full path
        target_dir (string): target dir
    """
    with zipfile.ZipFile(filename,"r") as zip_ref:
        zip_ref.extractall(target_dir)
#endregion




#region 解压 tar 文件到指定目录
def extract_tar(tar_file_path, target_directory):
    """
    解压 tar 文件到指定目录
    
    参数:
    tar_file_path (str): 要解压的 tar 文件路径
    target_directory (str): 解压到的目标目录
    """
    # 创建目标目录（如果不存在）
    os.makedirs(target_directory, exist_ok=True)
    
    # 打开 tar 文件并解压到目标目录
    with tarfile.open(tar_file_path, 'r') as tar:
        tar.extractall(path=target_directory)
    print("解压完成")
#endregion


def isfileExists(filepath):
   return os.path.exists(filepath)


def checkFileExists_with_dir_make(filepath):

    if os.path.exists(filepath) == False:
        file_dir= os.path.dirname(filepath)
        if not os.path.exists(file_dir):
            os.makedirs(file_dir)
        res =  False
    else:
        res = True

    return res


def cp_pdb(src, dst):
    """
    Copy a PDB file from source to destination, ensuring the destination directory exists.
    Returns:
        1: if the copy is successful.
        0: if the copy fails or the source file does not exist.
    """

    # if not os.path.exists(os.path.dirname(dst)):
    #     print(os.path.dirname(dst))
    #     os.makedirs(os.path.dirname(dst), exist_ok=True)
    
    # 复制文件
    if os.path.exists(src):
        
        if not os.path.exists(os.path.dirname(dst)):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
        
        shutil.copy(src, dst)
        return True
    else:
        return False
    
    
def get_pdb_files(directory):
    """
    Get a list of all .pdb files in the specified directory and its subdirectories.

    Parameters:
    directory (str): The root directory to search for .pdb files.

    Returns:
    list: A list of paths to .pdb files.
    """
    pdb_files = []
    for root, dirs, files in os.walk(directory):  # 遍历目录及子目录
        for file in files:
            if file.endswith(".pdb"):  # 筛选 .pdb 文件
                pdb_files.append(os.path.join(root, file))
    return pdb_files



def is_file_empty(file_path):
    """检查文件是否为空"""
    if os.path.isfile(file_path):
        size = os.path.getsize(file_path)
        if size == 0:
            return True
        else:
            return False
    else:
        raise FileNotFoundError(f"文件不存在: {file_path}")