import fire
from .bmz import (
    Reader,
    Writer,
    test_difference,
    extract
)
from .allc import (AllC, generate_ssi, allc2mz, prepare_sky,
                   merge_mz, extractCG, merge_cell_type)

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(
        {
            "Writer": Writer,
            'Reader': Reader,
            'AllC': AllC,
            'allc2mz': allc2mz,
            'test_diff': test_difference,
            'generate_ssi': generate_ssi,
            'prepare_sky': prepare_sky,
            'merge_mz': merge_mz,
            'merge_cell_type': merge_cell_type,
            'extract': extract,
            'extractCG': extractCG,
        }
	)

if __name__=="_main__":
	main()