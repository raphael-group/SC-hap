

# Requirements

    - Snakemake: https://snakemake.readthedocs.io/en/stable/
    - Wext: https://github.com/raphael-group/wext

# EXAMPLE

Run example as:
    `snakemake --configfile example/config --cores 8`

# USAGE

## Input file
    
    1. Data file. Comma-separated file, with one line per variant with the following columns:
        - `loc`: The genomic position of the variant
        - `allele_a, allele_b`
        - `a_i, b_i`: For each sample i, the count of the a allele and b allele respectively. 

        Example:
            
```
loc,allele_a,allele_b,a_0,b_0,a_1,b_1,a_2,b_2,a_3,b_3,a_4,b_4,a_5,b_5,a_6,b_6,a_7,b_7,a_8,b_8
1897,G,A,5,0,0,0,0,0,0,0,0,0,4,0,1,0,2,0
11754,C,T,1,0,0,0,0,0,0,3,0,0,4,0,8,1,0,0
13656,A,G,2,0,4,0,2,2,0,0,0,0,0,0,4,0,0,0
```

    2. Config file. YAML formatted file with the following entries:
        - `outdir`: output directory
        - `cutoff`: cutoff parameter(s)


## Output file
    1. `hapcut_input_f/mda_fragments.{CUTOFF}.txt`
        Output fragments in hapcut2 format
