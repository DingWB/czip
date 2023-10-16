import fire
from .bmz import (
	Reader,
	Writer
)
from .ballc import AllC

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(
		{
			"Writer": Writer,
			'Reader': Reader,
			'AllC':AllC,
		}
	)

if __name__=="_main__":
	main()