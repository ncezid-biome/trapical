from __future__ import annotations
import getopt, os, sys

class Config():
    FORMAT = "fasta"
    __DEF_FNA_DIR = "loci_fastas"
    __DEF_ALN_DIR = "loci_alignments"
    __DEF_CPUS = 1
    __DEF_OUTFN = "msa.aln"
    
    def __init__(self, version:str) -> Config:
        # type hint variables
        self.pulsenetDir:str
        self.fnaDir:str
        self.fnaFiles:list[str]
        self.alnDir:str
        self.alnFiles:list[str]
        self.outFn:str
        self.cpus:int
        
        # import version
        self.version:str = version
        
        # parse the arguments
        self.__parseArgs()
    
    def __parseArgs(self) -> None:
        IN_FLAGS = ('-i', '--in')
        FNA_FLAGS = ('-f', '--fna_dir')
        ALN_FLAGS = ('-a', '--aln_dir')
        OUT_FLAGS = ('-o', '--out')
        CPUS_FLAGS = ('-n', '--num_threads')
        VERS_FLAGS = ('-v', '--version')
        HELP_FLAGS = ('-h', '--help')
        SHORT_OPTS = IN_FLAGS[0][-1] + ':' + \
                     FNA_FLAGS[0][-1] + ':' + \
                     ALN_FLAGS[0][-1] + ':' + \
                     OUT_FLAGS[0][-1] + ':' + \
                     CPUS_FLAGS[0][-1] + ':' + \
                     VERS_FLAGS[0][-1] + \
                     HELP_FLAGS[0][-1]
        LONG_OPTS = (IN_FLAGS[1][2:] + '=',
                     FNA_FLAGS[1][2:] + '=',
                     ALN_FLAGS[1][2:] + '=',
                     OUT_FLAGS[1][2:] + '=',
                     CPUS_FLAGS[1][2:] + '=',
                     VERS_FLAGS[1][2:],
                     HELP_FLAGS[1][2:])
        
        # set default values
        self.pulsenetDir = None
        self.fnaDir = os.path.join(os.curdir, Config.__DEF_FNA_DIR)
        self.alnDir = os.path.join(os.curdir, Config.__DEF_ALN_DIR)
        self.cpus = Config.__DEF_CPUS
        self.outFn = os.path.join(os.curdir, Config.__DEF_OUTFN)
        
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        for opt,arg in opts:
            if opt in IN_FLAGS:
                self.pulsenetDir = arg
            
            elif opt in FNA_FLAGS:
                self.fnaDir = arg
            
            elif opt in ALN_FLAGS:
                self.alnDir = arg
            
            elif opt in OUT_FLAGS:
                self.outFn = arg
            
            elif opt in CPUS_FLAGS:
                self.cpus = int(arg)
        
        if self.pulsenetDir is None:
            raise BaseException('must specify an input directory')
    
        
    
    
    