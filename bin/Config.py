from __future__ import annotations
import getopt, os, shutil, sys

class Config():
    FORMAT = "fasta"
    __DEF_FNA_DIR = "loci_fastas"
    __DEF_ALN_DIR = "loci_alignments"
    __DEF_CPUS = 1
    __DEF_OUTFN = "msa.aln"
    __DEF_FORCE = False
    __DEF_HELP = False
    
    def __init__(self, version:str, author:str) -> Config:
        # type hint variables that will be populated by parseArgs
        self.pulsenetDir:str
        self.fnaDir:str
        self.alnDir:str
        self.outFn:str
        self.cpus:int
        self.helpRequested:bool
        self.__force:bool
        
        # type hint variables that will be populated during runtime
        self.fnaFiles:list[str]
        self.alnFiles:list[str]
        
        # import version and author
        self.__version:str = version
        self.__author:str = author
        
        # parse the commandline arguments
        self.__parseArgs()
    
    def __checkForClustalo() -> None:
        """checks that clustalo is installed and in the path

        Raises:
            BaseException: clustalo is not installed properly
        """
        # constant
        CLUSTALO = "clustalo"
        ERR_MSG = ' is not installed or not in the PATH'
        
        # make sure that isPcr is in the path
        if not shutil.which(CLUSTALO):
            raise BaseException(f"'{CLUSTALO}'{ERR_MSG}")
    
    def __checkExistingFile(self, path:str, forceFlag:str) -> None:
        """checks if a directory exists and removes if forcing

        Args:
            path (str): the path to check
            forceFlag (str): the string of the force flag

        Raises:
            FileExistsError: the directory exists and not forcing
        """
        ERR_MSG = f" file already exists; use '{forceFlag}' to overwrite"
        
        # if it exists
        if os.path.exists(path):
            # remove if forcing otherwise error
            if self.__force:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)
            else:
                raise FileExistsError(f"'{path}'{ERR_MSG}")
    
    def __parseArgs(self) -> None:
        """parses command line arguments

        Raises:
            NotADirectoryError: the input directory does not exist
            PermissionError: cannot write to the output file
            ValueError: cpus specified is not an int
            BaseException: must specify all required arguments
        """
        # flags
        IN_FLAGS = ('-i', '--in')
        FNA_FLAGS = ('-f', '--fna_dir')
        ALN_FLAGS = ('-a', '--aln_dir')
        OUT_FLAGS = ('-o', '--out')
        CPUS_FLAGS = ('-n', '--num_threads')
        FORCE_FLAGS = ('--force',)
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
                     FORCE_FLAGS[0][2:],
                     VERS_FLAGS[1][2:],
                     HELP_FLAGS[1][2:])
        
        # error messages
        ERR_MSG_1 = "input directory does not exist"
        ERR_MSG_2 = " is not an integer"
        ERR_MSG_3 = "must specify a value for '--in'"
        ERR_MSG_4 = "cannot write to "
        
        def printHelp() -> None:
            GAP = " "*4
            EOL = "\n"
            SEP = ', '
            DEF_OPEN = ' (default: '
            CLOSE = ')'
            WIDTH = max(map(len, LONG_OPTS)) + len(GAP) + 4
            HELP_MSG = f"{EOL}trapical rapidly aligns pulsenet isolate core alleles{EOL}" + \
                       f"{GAP}{self.__author}, 2024{EOL}" + \
                       f"{GAP}trapical v{self.__version}{EOL*2}" + \
                       f"usage:{EOL}" + \
                       f"{GAP}python3 trapical.py [-{SHORT_OPTS.replace(':','')}]{EOL*2}" + \
                       f"required arguments:{EOL}" + \
                       f"{GAP}{SEP.join(IN_FLAGS):<{WIDTH}}[directory] the directory containing the Pulsenet2.0 WGMLST results{EOL*2}" + \
                       f"optional arguments:{EOL}" + \
                       f"{GAP}{SEP.join(FNA_FLAGS):<{WIDTH}}[directory] the directory where locus fasta files should be written{DEF_OPEN}{self.fnaDir}{CLOSE}{EOL}" + \
                       f"{GAP}{SEP.join(ALN_FLAGS):<{WIDTH}}[directory] the directory where aligned locus files should be written{DEF_OPEN}{self.alnDir}{CLOSE}{EOL}" + \
                       f"{GAP}{SEP.join(OUT_FLAGS):<{WIDTH}}[file] the name of the output file{DEF_OPEN}{self.outFn}{CLOSE}{EOL}" + \
                       f"{GAP}{SEP.join(CPUS_FLAGS):<{WIDTH}}[int] the number of processors for parallel processing{DEF_OPEN}{Config.__DEF_CPUS}{CLOSE}{EOL}" + \
                       f"{GAP}{SEP.join(VERS_FLAGS):<{WIDTH}}prints the version{EOL}" + \
                       f"{GAP}{SEP.join(HELP_FLAGS):<{WIDTH}}prints this meessage{EOL}" + \
                       f"{GAP}{SEP.join(FORCE_FLAGS):<{WIDTH}}if specified, overwrite existing files and directories{EOL}"
            
            print(HELP_MSG)
        
        # set default values
        self.pulsenetDir = None
        self.helpRequested = Config.__DEF_HELP
        self.fnaDir = os.path.join(os.curdir, Config.__DEF_FNA_DIR)
        self.alnDir = os.path.join(os.curdir, Config.__DEF_ALN_DIR)
        self.cpus = Config.__DEF_CPUS
        self.outFn = os.path.join(os.curdir, Config.__DEF_OUTFN)
        self.__force = Config.__DEF_FORCE
        
        # print help if requested
        if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or sys.argv[1:] == []:
            self.helpRequested = True
            printHelp()
        
        # print version if requested
        elif VERS_FLAGS[0] in sys.argv or VERS_FLAGS[1] in sys.argv:
            self.helpRequested = True
            print(f'trapical v{self.__version}')
        
        # parse command line arguments
        else:
            # make sure clustalo is installed
            Config.__checkForClustalo()
            
            # determine if we are forcing
            if FORCE_FLAGS[0] in sys.argv:
                self.__force = True
            
            opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
            for opt,arg in opts:
                # parse input directory
                if opt in IN_FLAGS:
                    # make sure the directory exists
                    if not os.path.isdir(arg):
                        raise NotADirectoryError(ERR_MSG_1)
                    
                    # save the value
                    self.pulsenetDir = arg
                
                # parse fasta directory
                elif opt in FNA_FLAGS:
                    self.fnaDir = os.path.abspath(arg)
                
                # parse the alignment directory
                elif opt in ALN_FLAGS:
                    self.alnDir = os.path.abspath(arg)
                
                # parse the output file
                elif opt in OUT_FLAGS:
                    # save the value
                    self.outFn = os.path.abspath(arg)
                
                # parse the number of cpus
                elif opt in CPUS_FLAGS:
                    try:
                        self.cpus = int(arg)
                    except ValueError:
                        raise ValueError(f"'{arg}'{ERR_MSG_2}")

            # make sure that an input directory was specified
            if self.pulsenetDir is None:
                raise BaseException(ERR_MSG_3)
            
            # handle existing directories and files
            self.__checkExistingFile(self.fnaDir, FORCE_FLAGS[0])
            self.__checkExistingFile(self.alnDir, FORCE_FLAGS[0])
            self.__checkExistingFile(self.outFn, FORCE_FLAGS[0])
            
            # make sure output directory has write permissions
            if not os.access(os.path.dirname(self.outFn), os.W_OK):
                raise PermissionError(f"{ERR_MSG_4}'{self.outFn}'")
            
            # make the fasta directory
            os.mkdir(self.fnaDir)
            
            # make the alignment directory (if it isn't the same as the fasta directory)
            if os.path.abspath(self.alnDir) != os.path.abspath(self.fnaDir):
                os.mkdir(self.alnDir)
