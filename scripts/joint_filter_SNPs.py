'''
Imports 
'''
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import sys
from os import listdir
from os.path import isfile, join
import logging
import vcf

logging.basicConfig(format='%(asctime)s %(levelname)s:  %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
log = logging.getLogger(__name__)

def read_in_files(tumorfiles, normalfiles):
    '''
    Read in all files
    
    returns:
        acclist: list of unique identifiers for the samples
        scfiles: dataframes 
    '''
    log.debug("read_in_files START")
    log.info("Reading in files")

    acclist = [f.split("/")[-1][-2:] for f in tumorfiles]
    nacclist = [f.split("/")[-1][-2:] for f in normalfiles]
    scfiles = [pd.read_table(f).set_index('loc') for f in tumorfiles]
    
    nfiles = [pd.read_table(f).set_index('loc') for f in normalfiles]



    log.info("Read %d tumor files", len(scfiles))
    log.info("Read %d normal files", len(nfiles))
    log.debug("read_in_files DONE")
    
    return acclist, scfiles, nacclist, nfiles

def write_out_files(filename, df):

    log.debug("write_out_files START")
    log.info("Writing out results to %s", filename)
    df.to_csv(filename)

    log.debug("write_out_files DONE")
    

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-t",
        "--tumorfiles",
        metavar="TUMORFILE",
        nargs='*', required=True
    )

    parser.add_argument(
        "-n",
        "--normalfiles",
        metavar="NORMALFILE",
        nargs="*", required=True
    )

    parser.add_argument(
        "-o",
        "--outputfile",
        metavar="OUTPUT",
        help="output file",
        required=True
    )

    parser.add_argument(
        "-v", "--verbose",
        action = 'store_true',
        default=False,
        help="Log full debug output"
    )

    parser.add_argument(
       "--VCF", metavar = "VCFFILE",
        required=True
    )
    args = parser.parse_args()


    if args.verbose:
        log.setLevel(logging.DEBUG)
    else: 
        log.setLevel(logging.INFO)
    #if args.logtofile:
    #    log.basicConfig(filename=join(args.output_dir, "log.txt"))
    
    log.info("Running %s", sys.argv[1])
    log.info("Args: Tumor files: %s", args.tumorfiles)
    log.info("Args: Normal files: %s", args.normalfiles)
    log.info("Args: Output file: %s", args.outputfile)
    
    return args.tumorfiles, args.normalfiles, args.outputfile, args.VCF


from scipy.stats import norm, beta
gamma = 0.001
def beta_posterior_test(n_a, n_b, offset):
    """
    Determines if a locus should be considered heterozygous.
    Arguments:
        n_a (int): number of a alleles counted. Used as alpha parameter for the beta distribution
        n_b (int): number of b alleles counted. Used as beta parameter for the beta distribution
        gamma (float): parameter used for deciding heterozygosity; determined via a beta distribution
            with 1 - gamma confidence
    Returns:
        A boolean indicating whether or not the allele should be considered heterozygous.
    """

    if n_a == -1 or n_b == -1: return False

    p_lower = gamma / 2.0
    p_upper = 1 - p_lower

    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], n_a + 1, n_b + 1)
    o = offset
    return c_lower <= 0.5+o and c_upper >= 0.5-o


def beta_posterior_test_joint(x):
    As = x['a_count']
    Bs = x['b_count']

    return beta_posterior_test(As, Bs,0.1)




def get_heterozygous_locs(scfiles):
    ''' 
    Get a list of heterozygous loci and their alleles
    Returns:
        alleles           Dataframe containing the a and b alleles for each loc
    '''
    
    log.debug("get_heterozygous_locs START")
    log.info("Finding heterozygous loci")
    ### Filter to the heterozygous loci
    alleles = reduce(lambda x, y: pd.merge(x, y, left_on=['chr', 'a_allele', 'b_allele'], right_on=['chr', 'a_allele', 'b_allele'], left_index=True, right_index=True), scfiles)
    #alleles = all_merged[['loc']]
    alleles['a'] = (alleles[['a_x', 'a_y']] > 0).sum(axis=1)
    alleles['b'] = (alleles[['b_x', 'b_y']] > 0).sum(axis=1)
    is_heterozygous = (alleles[['a', 'b']] > 0).sum(axis=1) == 2
    alleles = alleles[is_heterozygous]
    

    log.info("Found %d heterozygous loci", len(alleles))
    log.debug("Assigning alleles")

    alleles = alleles[['a_allele', 'b_allele']]

    log.debug("get_heterozygous_locs DONE")
    return alleles
    

def filter_df(scfiles, acclist, nfiles, nacclist):
    '''
    Filters the dataframe to just the rows that are heterozygous and gets the 
    count of each allele
    '''

    scfiles = [s[~s.index.duplicated()] for s in scfiles]
    nfiles = [s[~s.index.duplicated()] for s in nfiles]
    files = scfiles + nfiles

    alleles = files[0][['chr', 'a_allele', 'b_allele']]
    alleles = alleles[~alleles.index.duplicated()]
    alleles['a_count']=0
    alleles['b_count']=0
    alleles['a_sample_count']=0
    alleles['b_sample_count']=0
    log.debug("filter_df START")
    #log.debug("Alleles head: %s", str(alleles.head()))
    Acols = []
    Bcols = []
    
    for acc, df in enumerate(files):
    
        alleles['a_'+str(acc)] = df['a']
        alleles['b_'+str(acc)] = df['b']
        Acols.append('a_'+str(acc))
        Bcols.append('b_'+str(acc))
        
        
    alleles = alleles.fillna(0)    
    alleles['a_count'] += alleles.apply(lambda x: sum(x[Acols]), axis=1) 
    alleles['b_count'] += alleles.apply(lambda x: sum(x[Bcols]), axis=1)     
        

    #alleles2 = alleles[(alleles['a_count'] > 0) & (alleles['b_count'] > 0) ]
    alleles2 = alleles

    log.debug("Merged dataframe #%d shape: %s", acc, str(alleles2.shape))
    log.debug("filter_df END")
    
    return alleles2


import sys
import argparse

def get_a_and_b_allele(vcffile, scfiles, nfiles):

    vcf_reader = vcf.Reader(open(vcffile))
    pos_allele_map = {}
    for i,record in enumerate(vcf_reader):
        pos_allele_map[str(record.POS)] = [str(record.REF), map(str,record.ALT)]

    def allele_a(x):
            return pos_allele_map[str(x.name)][0]


    def allele_b(x):
        return pos_allele_map[str(x.name)][1][0]
     
        return None
    
    def get_a_allele(x):
        return x[x['a_allele']]
    

    locs = None
    for f1 in scfiles + nfiles:

        f1['a_allele'] = f1.apply(allele_a, axis=1)
        f1['b_allele'] = f1.apply(allele_b, axis=1)
  
        f1['a'] = f1.apply(lambda x: x[x['a_allele']], axis=1)
        f1['b'] = f1.apply(lambda x: x[x['b_allele']], axis=1)

    
def main():
    tumorfiles, normalfiles, outfile, vcffile = parse_arguments()
    if len(tumorfiles) + len(normalfiles) == 0: 
        raise Error("No tumor or normal supplied")
    acclist, scfiles, nacclist, nfiles = read_in_files(tumorfiles, normalfiles)
    get_a_and_b_allele(vcffile, scfiles, nfiles) 
    
    log.info("Filtering single cell files for heterozygous SNPs")
    merged = filter_df(scfiles, acclist, nfiles, nacclist)
    merged = merged[merged.apply(beta_posterior_test_joint, axis=1)]
    write_out_files(outfile, merged)

if __name__ == '__main__':
    main()

