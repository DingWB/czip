import fire
from .cz import (
    Reader,
    Writer,
    extract
)
from .allc import (AllC, generate_ssi, bed2cz,
                   merge_cz, extractCG,
                   merge_cell_type,
                   combp, prepare_methylpy, agg_beta)

__version__ = "0.2.1"

def main():
    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire(
        {
            "Writer": Writer,
            'Reader': Reader,
            'AllC': AllC,
            'bed2cz': bed2cz,
            'generate_ssi': generate_ssi,
            'merge_cz': merge_cz,
            'merge_cell_type': merge_cell_type,
            'extract': extract,
            'extractCG': extractCG,
            "combp": combp,
            'prepare_methylpy': prepare_methylpy,
            'agg_beta': agg_beta,
        }
	)

if __name__=="_main__":
    main()
