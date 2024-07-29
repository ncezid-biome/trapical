from Bio import SeqIO
from Bio.Seq import Seq
from bin.Clock import Clock
from bin.Config import Config
from Bio.SeqRecord import SeqRecord
import glob, gzip, json, multiprocessing, os


def __parseOneJson(fn:str) -> dict[int,Seq]:
    """parses one allele json file

    Args:
        fn (str): a gzipped json file

    Returns:
        dict[int,Seq]: key=allele id; val=sequence
    """
    # initialize output
    out = dict()
    
    # parse the json file
    with gzip.open(fn, 'r') as fh:
        data:dict[str,dict] = json.load(fh)
    
    # for each called locus
    for val in data['calls'].values():
        # for each called allele
        for call in val:
            # get the id and sequence
            id = call['id']
            seq = Seq(call['seq'])

            # save them in the dictionary
            out[id] = seq
    
    return out


def __importSequences(directory:str, cpus:int) -> dict[int,Seq]:
    """gets sequences from parsed json files in parallel

    Args:
        directory (str): the directory containing all the results
        cpus (int): the number of parallel processes

    Returns:
        dict[int,Seq]: key=allele id; val-sequence
    """
    # constant
    JSON_FILE = "allele_calls.json.gz"
    
    # initialize variables
    out = dict()
    files = glob.glob(os.path.join(directory, "*", JSON_FILE))
    
    # parse files in parallel
    pool = multiprocessing.Pool(cpus)
    results = pool.map(__parseOneJson, files)
    pool.close()
    pool.join()
    
    # combine the results
    while results != []:
        out.update(results.pop())
    
    return out


def __parseOneAlleleFile(fn:str) -> tuple[str,dict[str,int]]:
    """parses one core gene allele calls file

    Args:
        fn (str): the file to parse

    Returns:
        tuple[str,dict[str,int]]: name, key=locus; val=dict: key=name; val=allele id
    """
    # initialize variables
    firstLine = True
    indices = dict()
    out = dict()
    
    # go through each line in the csv
    with gzip.open(fn,'rt') as fh:
        for line in fh:
            # split the line into a row
            row = line.rstrip().split(',')
            
            # link the headers to the indices of this file
            if firstLine:
                for idx in range(1,len(row)):
                    indices[idx] = row[idx]
                firstLine = False
            
            # extract the name and allele for each locus
            else:
                name = row[0]
                
                for idx in range(1,len(row)):
                    # if it cannot be coerced to an int, then it is '?'; skip
                    try:
                        out[indices[idx]] = int(row[idx])
                    except ValueError:
                        pass
    
    return name,out


def __importCoreAlleles(directory:str, cpus:int) -> dict[str,dict[str,int]]:
    """gets the allele ids for all core genes in the input genomes in parallel

    Args:
        directory (str): directory containing the results
        cpus (int): number of parallel processes

    Returns:
        dict[str,dict[str,int]]: key=locus; val=dict: key=name; val=allele id
    """
    # constant
    CORE_FILE = "calls_core_standard.csv.gz"
    
    # initialize variables
    out = dict()
    files = glob.glob(os.path.join(directory, "*", CORE_FILE))
    
    # open the pool and parse files in parallel
    pool = multiprocessing.Pool(cpus)
    results = pool.map(__parseOneAlleleFile, files)
    pool.close()
    pool.join()
    
    # combine the results
    while results != []:
        name,core = results.pop()
        
        # first level key should be locus; then name
        for locus,id in core.items():
            out[locus] = out.get(locus, dict())
            out[locus].update({name: id})
    
    # only keep loci that are present (not '?') in every genome
    return {k:v for k,v in out.items() if len(v) == len(files)}


def __writeFastas(outdir:str, core:dict[str,dict[str,int]], seqs:dict[int,Seq], frmt:str) -> list[str]:
    """writes a fasta file for each locus in the core genes

    Args:
        outdir (str): the output directory where fasta files will be written
        core (dict[str,dict[str,int]]): the dictionary produced by __importCoreAlleles
        seqs (dict[int,Seq]): the dictionary produced by __importSequences
        frmt
    """
    # constants
    EXT = ".fna"
    
    # initialize output
    out = list()
    
    # for each locus
    for locus in core.keys():
        # create the new filename and save it
        fn = os.path.join(outdir, locus + EXT)
        out.append(fn)
        
        # remove the file if it exists
        if os.path.exists(fn):
            os.remove(fn)
        
        # open the file
        with open(fn, 'a') as fh:
            # write each unique allele to the file
            for uid in {x for x in core[locus].values()}:
                SeqIO.write(SeqRecord(seqs[uid], str(uid), '', ''), fh, frmt)


def _parseAlleleCalls(config:Config) -> dict[str,dict[str,int]]:
    """parses pulsenet2.0 allele caller results

    Args:
        config (Config): a Config object

    Returns:
        dict[str,dict[str,int]]: key=locus; val=dict: key=genome name; val=hash
    """
    # messages
    MSG_1 = "importing sequence data"
    MSG_2 = "importing core alleles"
    MSG_3 = "writing a fasta file for each core locus"
    
    # initialize clock
    clock = Clock()
    
    # import all the sequences
    clock.printStart(MSG_1)
    seqs = __importSequences(config.pulsenetDir, config.cpus)
    clock.printDone()
    
    # import the core alleles from the csv files
    clock.printStart(MSG_2)
    core = __importCoreAlleles(config.pulsenetDir, config.cpus)
    clock.printDone()
    
    # write all the fastas to file
    clock.printStart(MSG_3)
    config.fnaFiles = __writeFastas(config.fnaDir, core, seqs, config.FORMAT)
    clock.printDone()
    
    return core
