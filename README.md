# CROssBAR-BioCypher-Migration
This is a repo for CROssBARv2 data to the Neo4j database via [BioCypher](https://github.com/saezlab/BioCypher). CROssBARv2 is, an extended and improved version of our previous work (for v1 please check [CROssBAR](https://github.com/cansyl/CROssBAR)), a heterogeneous general purpose biomedical knowledge graph (KG) based system.

### This repo is currently under development. Therefore, you may encounter some problems while replicating it. Feel free to open issue if you encounter any problems.

## Installation
The project uses [Poetry](https://python-poetry.org). You can install like this:

```
git clone https://github.com/HUBioDataLab/CROssBAR-BioCypher-Migration.git
cd CROssBAR-BioCypher-Migration
poetry install
```

Poetry will create a virtual environment according to your configuration (either centrally or in the project folder). You can activate it by running `poetry shell` inside the project directory.

### Note about pycurl
You may encounter an error when executing the UniProt adapter about the SSL
backend in pycurl: `ImportError: pycurl: libcurl link-time ssl backend (openssl)
is different from compile-time ssl backend (none/other)`

Should this happen, it can be fixed as described here:
https://stackoverflow.com/questions/68167426/how-to-install-a-package-with-poetry-that-requires-cli-args
by running `poetry shell` followed by `pip list`, noting the version of pycurl,
and then running `pip install --compile --install-option="--with-openssl"
--upgrade --force-reinstall pycurl==<version>` to provide the correct SSL
backend.
