import fire
from .bmz import (
	Reader,
	Writer
)
from .allc import AllC,allc2mz

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(
		{
			"Writer": Writer,
			'Reader': Reader,
			'AllC':AllC,
            'allc2mz':allc2mz,
		}
	)

if __name__=="_main__":
	main()