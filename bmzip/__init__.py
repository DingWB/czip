import fire
from .bmz import (
    Reader,
    Writer,
    test_difference
)
from .allc import AllC, generate_context_ssi, allc2mz

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(
        {
            "Writer": Writer,
            'Reader': Reader,
            'AllC': AllC,
            'allc2mz': allc2mz,
            'test_diff': test_difference,
            'generate_context_ssi': generate_context_ssi,
        }
	)

if __name__=="_main__":
	main()