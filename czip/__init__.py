import fire
from .cz import (
    Reader,
    Writer,
    extract
)
from .allc import (AllC, generate_ssi1, bed2cz, generate_ssi2,
                   merge_cz, extractCG, aggregate,
                   merge_cell_type,
                   combp, prepare_methylpy, intersect, annot_dmr,
                   methylpy_heatmap)

__version__ = "0.2.1"

def main():
    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire(
        {
            "Writer": Writer,
            'Reader': Reader,
            'AllC': AllC,
            'bed2cz': bed2cz,
            'generate_ssi1': generate_ssi1,
            'generate_ssi2': generate_ssi2,
            'merge_cz': merge_cz,
            'merge_cell_type': merge_cell_type,
            'extract': extract,
            'extractCG': extractCG,
            'aggregate': aggregate,
            "combp": combp,
            'prepare_methylpy': prepare_methylpy,
            'intersect': intersect,
            'annot_dmr': annot_dmr,
            'methylpy_heatmap': methylpy_heatmap,
        }
	)

if __name__=="_main__":
    main()
