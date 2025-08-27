'''
Date: 2025-06-24 13:56:30
LastEditors: Zhenkun Shi kingstdio@gmail.com
LastEditTime: 2025-06-24 13:59:52
FilePath: /RXNRECer/modules/llm/chat.py
'''
import openai
import json
from typing import Any, Dict, List, Optional


import openai
import httpx
from typing import Any, Optional


class Chat:
    def __init__(
        self,
        name: str,
        url: str,
        api_key: str,
        proxy: Optional[str] = None  # 例如: "http://127.0.0.1:7890"
    ):
        self.name = name
        self.url = url
        self.api_key = api_key

        # 设置代理客户端（如果提供）
        client_options = {
            "base_url": self.url,
            "api_key": self.api_key,
        }

        if proxy:
            self.proxy_client = httpx.Client(proxies=proxy)
            client_options["http_client"] = self.proxy_client

        self.client = openai.Client(**client_options)

    def chat(
        self,
        message: str,
        system_prompt: Optional[str] = None,
        temperature: float = 0,
        max_tokens: int = 10240,
        response_format: Optional[dict] = None,
    ) -> Any:
        # 可选参数封装
        extra_kwargs = {}
        if response_format:
            extra_kwargs["response_format"] = response_format

        # 拼接消息
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": message})

        # 发起请求
        resp = self.client.chat.completions.create(
            model=self.name,
            messages=messages,
            temperature=temperature,
            max_tokens=max_tokens,
            **extra_kwargs
        )
        return resp
