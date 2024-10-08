from bin.Clock import Clock
from bin.Config import Config
from bin.createMsa import _createMsa
from bin.parseAlleleCalls import _parseAlleleCalls


__author__ = "Joseph S. Wirth"
__version__ = "0.1.0"


def _main() -> None:
    """main runner function
    """
    # start the clock
    clock = Clock()
    
    try:
        # parse command line arguments and build the Config object
        config = Config(__version__, __author__)  
        
        # only do work if help wasn't requested
        if not config.helpRequested:
            # get the core loci allele calls
            core = _parseAlleleCalls(config)
            
            # create the multiple sequence alignment
            _createMsa(config, core)
            
            # print the status
            print(f'\ntotal runtime: {clock.getTimeString()}')
    
    except Exception as e:
        Clock._killWheel()
        raise Exception(e)
