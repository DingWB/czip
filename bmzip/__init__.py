import fire
from .bmz import (
	Reader,
	Writer
)

def main():
	fire.core.Display = lambda lines, out: print(*lines, file=out)
	fire.Fire(
		{
			"Writer": Writer,
			'Reader': Reader,
		}
	)

if __name__=="_main__":
	main()