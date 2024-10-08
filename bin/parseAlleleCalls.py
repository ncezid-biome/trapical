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
            seq = Seq(call['seq']).upper()

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


def __writeFastas(outdir:str, core:dict[str,dict[str,int]], seqs:dict[int,Seq], frmt:str) -> tuple[list[str],dict[str,dict[str,int]]]:
    """writes a fasta file for each locus in the core genes

    Args:
        outdir (str): the output directory where fasta files will be written
        core (dict[str,dict[str,int]]): the dictionary produced by __importCoreAlleles
        seqs (dict[int,Seq]): the dictionary produced by __importSequences
        frmt (str): the sequence file format to write
    
    Returns:
        tuple[
            list[str]:                a list of fasta files for each locus with >1 allele
            dict[str,dict[str,int]]:  key=filename; val=dict: key=character; val=count
        ]
    """
    # constants
    EXT = ".fna"
    
    # initialize variables
    fastasToAlign = list()
    conservedCounts = dict()
    
    # for each locus
    for locus in core.keys():
        # get the hashes
        hashes = {x for x in core[locus].values()}
        
        # create the new filename
        fn = os.path.join(outdir, locus + EXT)
        
        # only add a fasta to the list if there are multiple alleles to align
        if len(hashes) > 1:
            fastasToAlign.append(fn)
        
        # count the conserved characters for loci with only one allele
        else:
            # initialize the entry for this file
            conservedCounts[fn] = dict()
            
            # for each allele
            for uid in hashes:
                # get all the unique characters
                chars = {str(x) for x in seqs[uid]}
                
                # save the counts for this file
                conservedCounts[fn] = {c:seqs[uid].count(c) for c in chars}
        
        # write the fasta for this locus
        with open(fn, 'a') as fh:
            # write each unique allele to the file
            for uid in hashes:
                SeqIO.write(SeqRecord(seqs[uid], str(uid), '', ''), fh, frmt)

    return fastasToAlign, conservedCounts


def _parseAlleleCalls(config:Config) -> tuple[dict[str,dict[str,int]],dict[str,dict[str,int]]]:
    """parses pulsenet2.0 allele caller results

    Args:
        config (Config): a Config object

    Returns:
        dict[str,dict[str,int]]: key=locus; val=dict: key=genome name; val=hash
        dict[str,dict[str,int]]: key=filename; val=dict: key=character; val=count
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
    config.fnaFiles, config.conservedCounts = __writeFastas(config.fnaDir, core, seqs, config.FORMAT)
    clock.printDone()
    
    return core
