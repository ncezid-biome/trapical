from Bio import SeqIO
from Bio.Seq import Seq
from bin.Clock import Clock
from typing import Generator
from bin.Config import Config
from Bio.SeqRecord import SeqRecord
import multiprocessing, os, subprocess


def __alignOneSequence(fna:str, alnDir:str) -> str:
    """uses clustalo to align one sequence

    Args:
        fna (str): the fasta file to
        alnDir (str): the directory where the fasta file will be saved

    Raises:
        subprocess.CalledProcessError: error if clustalo fails

    Returns:
        str: the filename of the aligned file
    """
    # constants
    CMD = 'clustalo'
    IN  = '--in'
    OUT = '--out'
    EXT = '.aln'
    SUFFIX = ['--threads', '1', '--force']
    ERR_MSG = 'clustalo failed'
    
    # create the alignment filename
    aln = os.path.join(alnDir, os.path.splitext(os.path.basename(fna))[0] + EXT)
    
    # build the command
    cmd = [CMD, IN, fna, OUT, aln]
    cmd.extend(SUFFIX)
    
    # run clustalo
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        raise subprocess.CalledProcessError(ERR_MSG)
    
    return aln


def __alignAllSequences(fnas:list[str], alnDir:str, cpus:int) -> list[str]:
    """aligns all the sequences in parallel

    Args:
        fnas (list[str]): a list of fasta files to align
        alnDir (str): the directory where alignments will be saved
        cpus (int): the number of cpus to use for parallel processing

    Raises:
        ChildProcessError: a call to clustalo failed

    Returns:
        list[str]: a list of alignment files
    """
    # message
    ERR_MSG = 'one or more calls to clustalo failed'
    
    # generator to get arguments
    def genArgs() -> Generator[tuple[str,str],None,None]:
        for fn in fnas:
            yield (fn, alnDir)
    
    # align sequences in parallel
    pool = multiprocessing.Pool(cpus)
    out = pool.starmap(__alignOneSequence, genArgs())
    pool.close()
    pool.join()
    
    # make sure everything worked
    if len(out) != len(fnas):
        raise ChildProcessError(ERR_MSG)
    
    return out


def __getVariableSitesForOneLocus(fn:str, frmt:str) -> tuple[str,list[int]]:
    """gets a list of the sites in an aligned file that are variable

    Args:
        fn (str): the alignment file to evaluate
        frmt (str): the sequence format

    Returns:
        tuple[str,list[int]]: the filename, a list of indices of the string
    """
    # initialize variables
    out = list()
    
    # load all the sequence data into memory
    seqs = [rec.seq for rec in SeqIO.parse(fn, frmt)]
    
    # for each position in the alignment
    for position in range(len(seqs[0])):
        # get the character for the first sequence
        char = seqs[0][position]
        
        # for each other sequence
        for idx in range(1, len(seqs)):
            # save the position if the character is different; stop looping
            if seqs[idx][position] != char:
                out.append(position)
                break

    return fn, out
    

def __getAllVariableSites(alns:list[str], frmt:str, cpus:int) -> dict[int,Seq]:
    """gets all the variable sites for each locus in parallel

    Args:
        alnFiles (list[str]): a list of alignment files
        cpus (int): the number of cpus for parallel processing

    Returns:
        dict[int,Seq]: key=allele hash code; val=sequence
    """
    # initialize variables
    out = dict()
    rec:SeqRecord
    
    # create argument list
    args = [(f, frmt) for f in alns]
    
    # get the variable positions in parallel
    pool = multiprocessing.Pool(cpus)
    results = pool.starmap(__getVariableSitesForOneLocus, args)
    pool.close()
    pool.join()
    
    # use the variable positions to extract the sequences
    for fn,sites in results:
        for rec in SeqIO.parse(fn, frmt):
            out[int(rec.id)] = Seq('').join([rec.seq[x] for x in sites])
    
    return out


def __restructureCore(core:dict[str,dict[str,int]]) -> dict[str,list[int]]:
    """restructures the core alleles to remove loci names as a key

    Args:
        core (dict[str,dict[str,int]]): key=locus; val=dict: key=genome name; val=allele hash code

    Returns:
        dict[str,list[int]]: key=genome name; val=list of allele hashcodes
    """
    # initialize output
    out = dict()
    
    # for each locus
    for locus in core.keys():
        # for each genome name
        for name in core[locus].keys():
            # add the hash code to the list
            try:
                out[name].append(core[locus][name])
            
            # create a new list if one does not exist yet
            except KeyError:
                out[name] = [core[locus][name]]
    
    return out


def __writeVariableSites(fn:str, core:dict[str,list[int]], sites:dict[int,Seq], frmt:str) -> None:
    """writes concatenated variable sites to file

    Args:
        fn (str): the file to write to
        core (dict[str,list[int]]): the dictionary produced by _restructureCore
        sites (dict[int,Seq]): the dictionary produced by _getAllVariableSites
        frmt (str): the output sequence format
    """
    # open the file in append mode
    with open(fn, 'a') as fh:
        # for each genome name
        for name in core.keys():
            # extract the variable sequences and concatenate them
            seq = Seq('').join([sites[x] for x in core[name] if x in sites.keys()])
            
            # write the concatenated sequence to file
            SeqIO.write(SeqRecord(seq, name, '', ''), fh, frmt)


def _createMsa(config:Config, core:dict[str,dict[str,int]]) -> None:
    """creates a multiple sequence alignment

    Args:
        config (Config): a Config object
        core (dict[str,dict[str,int]]): key=locus; val=dict: key=genome name; val=hash
    """
    # messages
    MSG_1 = 'aligning sequences'
    MSG_2 = 'determining variable sites'
    MSG_3A = 'writing alignment to file ('
    MSG_3B = ' characters per isolate)'
    
    # initialize clock
    clock = Clock()
    
    # align sequences
    clock.printStart(MSG_1)
    config.alnFiles = __alignAllSequences(config.fnaFiles, config.alnDir, config.cpus)
    clock.printDone()
    
    # get the variable sites for each locus
    clock.printStart(MSG_2)
    sites = __getAllVariableSites(config.alnFiles, config.FORMAT, config.cpus)
    clock.printDone()
    
    # write results to file
    clock.printStart(f'{MSG_3A}{sum(map(len, sites.values()))}{MSG_3B}')
    core = __restructureCore(core)
    __writeVariableSites(config.outFn, core, sites, config.FORMAT)
    clock.printDone()
