# Chunk ZIP

## Installation

```shell
pip install git+http://github.com/DingWB/czip

python setup.py install
```

## Implementation

|                                  | allcools | ballcools | bmzip |
| -------------------------------- | -------- | --------- | ----- |
| Format                           | .tsv.gz  | .ballc    | .mz   |
| Compression algorithm            | bgzip    | bgzip     | bmzip |
| Support Random Access ?          | Yes      | Yes       | Yes   |
| Need extra index file for query? | Yes      | yes       | No    |
| Quickly Merge?                   | No       | No        | Yes   |

![img.png](docs/images/tab1.png)

![docs/images/img.png](docs/images/img.png)

## Usage

[Documentation](https://dingwb.github.io/czip)