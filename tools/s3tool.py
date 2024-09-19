import sys,os

sys.path.append(os.path.dirname(os.path.realpath('__file__')))
sys.path.append('../../')
from config import s3conf as s3cfg
import os
import boto3
from boto3.s3.transfer import TransferConfig
import tools.filetool as fileTool
import tempfile
from botocore.exceptions import ClientError
import pandas as pd
import logging

# 配置日志记录器
logging.basicConfig(level=logging.INFO)  # 设置日志级别为 INFO
# 创建日志记录器
logger = logging.getLogger(__name__)


#region  1. 保存 pandas DataFrame 到 S3 存储空间
def save_df_to_s3(df, file_key, format='feather', bucket_name=s3cfg.BUCKET_NAME, endpoint_url=s3cfg.ENDPOINT_URL, aws_access_key_id=s3cfg.ACCESS_KEY_ID, aws_secret_access_key=s3cfg.SECRET_ACCESS_KEY):
    """
    保存 pandas DataFrame 到 S3 存储空间

    Args:
        df (DataFrame): pandas DataFrame
        file_key (str): 保存文件名及保存路径
        format (str, optional): 文件格式，默认为'feather'。
        bucket_name (str, optional): 存储桶的名称。默认为 s3cfg.BUCKET_NAME。
        endpoint_url (str, optional): 终端节点 URL。默认为 s3cfg.ENDPOINT_URL。
        aws_access_key_id (str, optional): AWS 访问密钥 ID。默认为 s3cfg.AWS_ACCESS_KEY_ID。
        aws_secret_access_key (str, optional): AWS 秘密访问密钥。默认为 s3cfg.AWS_SECRET_ACCESS_KEY。
    """

    s3 = boto3.client('s3', endpoint_url=endpoint_url,
                      aws_access_key_id=aws_access_key_id,
                      aws_secret_access_key=aws_secret_access_key)
    GB = 1024 ** 3
    config = TransferConfig(multipart_threshold=5*GB)

    temp_file = None

    try:
        if format == 'feather':
            temp_file = tempfile.NamedTemporaryFile(suffix='.feather', delete=False)
            df.to_feather(temp_file.name)
        elif format == 'hdf':
            temp_file = tempfile.NamedTemporaryFile(suffix='.hdf', delete=False)
            df.to_hdf(temp_file.name, key='data')
        else:  # 默认为 CSV 格式
            temp_file = tempfile.NamedTemporaryFile(suffix='.csv', delete=False)
            df.to_csv(temp_file.name, index=False)

        s3.upload_file(temp_file.name, bucket_name, file_key, Config=config)

        logger.info(f"文件上传成功。{temp_file.name}")

    except Exception as e:
        print(f'上传文件失败：{e}')

    finally:
        if temp_file:
            temp_file.close()  # 关闭并删除临时文件
#endregion              
            
#region 2. 将字符串以文件的形式保存到S3
def save_string_to_s3(text_string, file_key, format='html', bucket_name=s3cfg.BUCKET_NAME, endpoint_url=s3cfg.ENDPOINT_URL, aws_access_key_id=s3cfg.ACCESS_KEY_ID, aws_secret_access_key=s3cfg.SECRET_ACCESS_KEY):
    """将字符串以文件的形式保存到S3

    Args:
        text_string (string): 字符串
        file_key (string): 保存文件名及保存路径
        format (str, optional): _description_. Defaults to 'html'.
        bucket_name (_type_, optional): _description_. Defaults to s3cfg.BUCKET_NAME.
        endpoint_url (_type_, optional): _description_. Defaults to s3cfg.ENDPOINT_URL.
        aws_access_key_id (_type_, optional): _description_. Defaults to s3cfg.AWS_ACCESS_KEY_ID.
        aws_secret_access_key (_type_, optional): _description_. Defaults to s3cfg.AWS_SECRET_ACCESS_KEY.
    """

    s3 = boto3.client('s3', endpoint_url=endpoint_url,
                      aws_access_key_id=aws_access_key_id,
                      aws_secret_access_key=aws_secret_access_key)
    GB = 1024 ** 3
    config = TransferConfig(multipart_threshold=5*GB)

    temp_file = None

    metadata = {'Content-Type': f'txt'}
    metadata ={"ContentType":"txt"}
    try:
        if format == 'html':
            temp_file = tempfile.NamedTemporaryFile(suffix='.html', delete=False)
            metadata =  {"ContentType":"application/html"}
        elif format == 'jpg':
            temp_file = tempfile.NamedTemporaryFile(suffix='.jpg', delete=False)
            metadata = {"ContentType":"image/jpg"}
        elif format == 'svg':
            metadata = {"ContentType":"image/svg+xml"}
            temp_file = tempfile.NamedTemporaryFile(suffix='.svg', delete=False)
        else:  # 默认为 txt 格式
            temp_file = tempfile.NamedTemporaryFile(suffix='.txt', delete=False)

            
        with open(temp_file.name, 'w') as file:
            file.write(text_string)

        s3.upload_file(temp_file.name, bucket_name, file_key, Config=config, ExtraArgs=metadata)
        
        logger.info(f"文件上传成功。{temp_file.name}")

    except Exception as e:
        print(f'上传文件失败：{e}')

    finally:
        if temp_file:
            temp_file.close()  # 关闭并删除临时文件        
#endregion            

#region 3. 从S3存储空间读取文件到 pandas DataFrame
def read_df_from_s3(file_key, format='feather', bucket_name=s3cfg.BUCKET_NAME, endpoint_url=s3cfg.ENDPOINT_URL, aws_access_key_id=s3cfg.ACCESS_KEY_ID, aws_secret_access_key=s3cfg.SECRET_ACCESS_KEY):
    """
    从 S3 存储空间读取文件到 pandas DataFrame

    Args:
        file_key (str): 文件在 S3 中的键值
        format (str, optional): 文件格式，默认为'feather'。
        bucket_name (str, optional): 存储桶的名称。默认为 s3cfg.BUCKET_NAME。
        endpoint_url (str, optional): 终端节点 URL。默认为 s3cfg.ENDPOINT_URL。
        aws_access_key_id (str, optional): AWS 访问密钥 ID。默认为 s3cfg.AWS_ACCESS_KEY_ID。
        aws_secret_access_key (str, optional): AWS 秘密访问密钥。默认为 s3cfg.AWS_SECRET_ACCESS_KEY。

    Returns:
        DataFrame or None: 读取的 DataFrame 数据或者为 None
    """

    try:
        s3 = boto3.client('s3', endpoint_url=endpoint_url,
                          aws_access_key_id=aws_access_key_id,
                          aws_secret_access_key=aws_secret_access_key)

        with tempfile.NamedTemporaryFile(suffix=f'.{format}', delete=False) as temp_file:
            temp_file_path = temp_file.name

            s3.download_file(bucket_name, file_key, temp_file_path)
            # print(f"文件 '{file_key}' 下载到 '{temp_file_path}'")

            if format == 'feather':
                df = pd.read_feather(temp_file_path)
            elif format == 'hdf':
                df = pd.read_hdf(temp_file_path)
            else:  # 默认为 CSV 格式
                df = pd.read_csv(temp_file_path)

            os.remove(temp_file_path) # 删除临时文件
            return df

    except ClientError as e:
        if e.response['Error']['Code'] == '404':
            print(f"文件 '{file_key}' 不存在。")
        else:
            print(f"下载文件出错：{e}")
        return None
    except Exception as ex:
        print(f"错误：{ex}")
        return None
#endregion

#region 4. 上传本地文件到S3            
def upload_file_to_s3(file, file_key, format='', bucket_name=s3cfg.BUCKET_NAME, endpoint_url=s3cfg.ENDPOINT_URL, aws_access_key_id=s3cfg.ACCESS_KEY_ID, aws_secret_access_key=s3cfg.SECRET_ACCESS_KEY):
    s3 = boto3.client('s3', endpoint_url=endpoint_url,
                      aws_access_key_id=aws_access_key_id,
                      aws_secret_access_key=aws_secret_access_key)
    GB = 1024 ** 3
    config = TransferConfig(multipart_threshold=5*GB)
    
   
    try:
        if format == 'html':
            metadata =  {"ContentType":"application/html"}
        elif format == 'jpg':
            metadata = {"ContentType":"image/jpg"}
        elif format == 'svg':
            metadata = {"ContentType":"image/svg+xml"}
        elif format == 'txt':
            metadata = {"ContentType":"text/plain"}
        elif format == 'pdf':
            metadata = {"ContentType":"application/pdf"}
        elif format == 'json':
            metadata = {"ContentType":"application/json"}
        else:  # 默认为空
            metadata ={}
            
    
        s3.upload_file(file, bucket_name, file_key, Config=config,  ExtraArgs=metadata)
        # logger.info(f"文件上传成功。{file}")
    except Exception as e:
        print(f'上传文件失败：{e}')
#endregion

#region 5. 从S3下载文件到本地    
def download_file_from_s3(save_path, 
                          file_key, 
                          bucket_name=s3cfg.BUCKET_NAME, 
                          endpoint_url=s3cfg.ENDPOINT_URL, 
                          aws_access_key_id=s3cfg.ACCESS_KEY_ID, 
                          aws_secret_access_key=s3cfg.SECRET_ACCESS_KEY):
    s3 = boto3.client('s3', endpoint_url=endpoint_url,
                      aws_access_key_id=aws_access_key_id,
                      aws_secret_access_key=aws_secret_access_key)
    GB = 1024 ** 3
    config = TransferConfig(multipart_threshold=5*GB)
    s3.download_file(bucket_name, file_key, save_path, Config=config)
    logger.info(f'文件下载成功。保存路径: {save_path}')                
#endregion



def save_files_and_upload_to_s3(file_content:list,
                                local_path:list,
                                s3_file_keys:list,
                                format='',
                                bucket_name=s3cfg.BUCKET_NAME, 
                                endpoint_url=s3cfg.ENDPOINT_URL, 
                                aws_access_key_id=s3cfg.ACCESS_KEY_ID, 
                                aws_secret_access_key=s3cfg.SECRET_ACCESS_KEY
                                ):
    
    # print(local_path)
    for filecontent, lpath, s3path in zip(file_content, local_path, s3_file_keys):
        if fileTool.checkFileExists_with_dir_make(lpath) == True:
            logger.info(f'File exists: {local_path}')
        
        with open(lpath, mode='w') as f:
            f.write(filecontent)
            
        upload_file_to_s3(file=lpath,
                        file_key=s3path,
                        format=format,
                        bucket_name=bucket_name, 
                        endpoint_url=endpoint_url,
                        aws_access_key_id=aws_access_key_id, 
                        aws_secret_access_key=aws_secret_access_key)

        # logger.info('save success')
        
        
        

if __name__ == '__main__':
    print('success')