# RXNRECer-S3 LLM API Setup

RXNRECer-S3 is the optional explanation layer of RXNRECer. It takes candidate reactions from RXNRECer-S1/S2 and uses a large language model (LLM) to generate structured, human-readable rationales for expert review.

S3 is not required for standard S1/S2 prediction. Use it when you want natural-language evidence summaries for candidate enzyme-reaction annotations.

## Quickstart

1. Install RXNRECer and download the required data and model files.

```bash
pip install rxnrecer
rxnrecer-download-data
```

2. Configure an OpenAI-compatible chat-completions endpoint.

```bash
export LLM_API_URL="https://api.example.com/v1"
export LLM_API_KEY="your_api_key_here"
export LLM_MODEL="provider/model-version"
```

3. Run S3 mode.

```bash
rxnrecer -i input.fasta -o output.json -m s3 -f json
```

For reproducible runs, use a fixed provider, a fixed model name/version, and deterministic decoding settings when the provider supports them.

## Configuration Variables

RXNRECer currently reads the following environment variables for S3:

| Variable | Required | Description |
| --- | --- | --- |
| `LLM_API_URL` | Yes | Base URL for an OpenAI-compatible API, usually ending in `/v1`. |
| `LLM_API_KEY` | Yes | API key or local-server token. For unauthenticated local servers, use a placeholder such as `EMPTY`. |
| `LLM_MODEL` | Recommended | Exact model name or version exposed by the provider. Defaults to `openai/gpt-4.1` if unset. |

The underlying client uses the OpenAI chat-completions format. Provider-specific SDKs are not needed when the endpoint implements that protocol.

## Hosted API Example

Many hosted gateways expose an OpenAI-compatible endpoint. The exact model name depends on the provider.

```bash
export LLM_API_URL="https://api.provider.example/v1"
export LLM_API_KEY="sk-..."
export LLM_MODEL="provider/model-version"

rxnrecer -i input.fasta -o output.json -m s3 -f json
```

OpenRouter-style endpoints usually follow the same pattern:

```bash
export LLM_API_URL="https://openrouter.ai/api/v1"
export LLM_API_KEY="sk-or-v1-..."
export LLM_MODEL="openai/gpt-4.1"
```

For beginners, general provider-side setup guides such as <https://docs.aimlapi.com/> can be useful for learning where to find API keys, base URLs, and model identifiers. Always use the values provided by your actual provider.

## Local vLLM Server

RXNRECer-S3 can also use a local vLLM server when it exposes the OpenAI-compatible API.

Start a local server, for example:

```bash
python -m vllm.entrypoints.openai.api_server \
  --model /path/to/local/model \
  --served-model-name local-rxnrecer-llm \
  --host 127.0.0.1 \
  --port 8000
```

Then point RXNRECer to that server:

```bash
export LLM_API_URL="http://127.0.0.1:8000/v1"
export LLM_API_KEY="EMPTY"
export LLM_MODEL="local-rxnrecer-llm"

rxnrecer -i input.fasta -o output.json -m s3 -f json
```

## Local SGLang Server

SGLang can be used in the same way if it is launched with an OpenAI-compatible frontend.

Example:

```bash
python -m sglang.launch_server \
  --model-path /path/to/local/model \
  --served-model-name local-rxnrecer-llm \
  --host 127.0.0.1 \
  --port 30000
```

Configure RXNRECer:

```bash
export LLM_API_URL="http://127.0.0.1:30000/v1"
export LLM_API_KEY="EMPTY"
export LLM_MODEL="local-rxnrecer-llm"

rxnrecer -i input.fasta -o output.json -m s3 -f json
```

Use the exact startup options supported by your installed vLLM or SGLang version.

## Python API

For scripted workflows, set the configuration values before running S3:

```python
from rxnrecer.config import config as cfg

cfg.LLM_API_URL = "http://127.0.0.1:8000/v1"
cfg.LLM_API_KEY = "EMPTY"
cfg.LLM_MODEL = "local-rxnrecer-llm"
```

Advanced users can also pass a model name directly to the lower-level S3 helper:

```python
from rxnrecer.lib.llm import qa

result = qa.batch_chat(
    res_rxnrecer=s3_input_df,
    api_key="EMPTY",
    api_url="http://127.0.0.1:8000/v1",
    llm_model="local-rxnrecer-llm",
)
```

## Reproducibility Notes

LLM outputs can change across providers, model releases, decoding settings, and prompt versions. For reproducible S3 analyses, record:

- API provider or local inference engine.
- Exact model name/version.
- API base URL or local deployment identifier.
- RXNRECer version.
- Prompt template version, if customized.
- Inference date.
- Input FASTA and generated S3 input JSON.
- Raw LLM response and parsed RXNRECer output.
- Decoding settings if exposed by the provider, such as temperature, top-p, seed, max tokens, or JSON mode.

When possible, use deterministic settings such as `temperature=0` and a fixed model version.

## Interpreting S3 Output

S3 explanations are intended to help users understand candidate reactions returned by RXNRECer-S1/S2. They should be treated as evidence summaries for expert interpretation, not as experimental validation.

Recommended practice:

- Verify important claims against curated databases, deterministic sequence-analysis tools, structural analyses, or biochemical experiments.
- Keep raw API responses for auditability.
- Reject outputs that cite evidence not present in the input context or known external sources.
- Do not use S3 explanations as ground-truth labels for benchmark evaluation.

## Troubleshooting

### Authentication error

Check that `LLM_API_KEY` is set and accepted by the provider. For local servers without authentication, use the placeholder expected by the server or client, such as `EMPTY`.

### Model not found

Check that `LLM_MODEL` exactly matches the model name exposed by the hosted provider, vLLM server, or SGLang server.

### Connection refused

For local APIs, confirm that the server is running and that `LLM_API_URL` includes the correct host, port, and `/v1` suffix.

### Output is not valid JSON

Retry with the same model and deterministic decoding settings. If the provider supports JSON mode or schema-constrained output, enable it at the server or gateway level.

### Different models give different rationales

This is expected for generative models. For studies that need reproducibility, choose one provider/model version and keep it fixed for all S3 runs.
