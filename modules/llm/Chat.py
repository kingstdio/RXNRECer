import openai


class Chat:
    def __init__(self, name, url, api_key):
        self.name = name
        self.url = url 
        self.api_key = api_key
        self.client =  openai.Client(base_url = self.url, api_key= self.api_key)
        
    def chat(self, message, temperature=0, max_tokens=2048):
        
        response = self.client.chat.completions.create(
            model= self.name,  
            messages=[
                {
                    "role": "user", 
                    "content": message
                },
            ],
            temperature=temperature,
            max_tokens=max_tokens,
        )



        return response
        