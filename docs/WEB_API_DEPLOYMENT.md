# Public Webserver and Local Usage

RXNRECer can be used through the public webserver or installed locally as a Python package.

## Public Webserver

The public RXNRECer webserver is available at:

<https://rxnrecer.biodesign.ac.cn/>

The webserver provides a browser-based interface for testing RXNRECer without local installation. It is intended for interactive use, quick checks, demonstrations, and users who prefer not to configure a local Python environment.

Typical webserver workflow:

1. Open <https://rxnrecer.biodesign.ac.cn/>.
2. Submit protein sequences through the web interface.
3. Run RXNRECer prediction.
4. Review reaction-level predictions and supporting information in the browser.
5. Download or save results for downstream review when needed.

## Local CLI/Python Workflow

For reproducible local analysis, large-scale batch inference, or integration into custom pipelines, install RXNRECer and run it locally:

```bash
pip install rxnrecer
rxnrecer-download-data

rxnrecer -i input.fasta -o output.tsv -m s1
rxnrecer -i input.fasta -o output.tsv -m s2
rxnrecer -i input.fasta -o output.json -m s3 -f json
```

The local workflow is the recommended path when users need full control over input files, package versions, downloaded data/model bundles, LLM settings, and generated outputs.

## How the Two Modes Fit Together

The webserver and the local package are complementary:

- Use the public webserver for convenient interactive testing and public access.
- Use the CLI for reproducible batch jobs and command-line workflows.
- Use the Python API when integrating RXNRECer into scripts, notebooks, or institutional pipelines.
- Use S3 mode when structured LLM-based explanations are needed; see [LLM API Setup](LLM_API_SETUP.md) for local S3 configuration.

## Reproducibility Notes

For local runs, record:

- RXNRECer version.
- Data and model bundle version.
- Input FASTA file.
- Selected mode (`s1`, `s2`, or `s3`).
- Hardware and runtime configuration for large jobs.
- Output TSV/JSON files.
- LLM provider, model name/version, and decoding settings if S3 is enabled.

For webserver jobs, users should save downloaded outputs and record the submission date, selected analysis options, and input sequences. When exact reruns are required, the local CLI/Python workflow provides the most controlled execution path.
