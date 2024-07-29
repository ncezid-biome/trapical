from bin.Clock import Clock
from bin.Config import Config
from bin.createMsa import _createMsa
from bin.parseAlleleCalls import _parseAlleleCalls

__author__ = "Joseph S. Wirth"
__version__ = "0.0.1"


def _main() -> None:
    clock = Clock()
    
    config = Config(__version__)
    
    seqs,core = _parseAlleleCalls(config)
    
    _createMsa(config, seqs, core)
    
    print(f'\ntotal runtime: {clock.getTimeString()}')
