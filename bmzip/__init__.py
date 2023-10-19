import fire
from .bmz import (
    Reader,
    Writer,
    test_difference
)
from .allc import AllC, allc2mz

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(
        {
            "Writer": Writer,
            'Reader': Reader,
            'AllC': AllC,
            'allc2mz': allc2mz,
            'test_diff': test_difference,
        }
	)

if __name__=="_main__":
	main()