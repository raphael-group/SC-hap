SCRIPT_DIR  = './scripts/'

###
#   Config
###

### Set of chromosomes to run on
CUTOFF = config['cutoff']
OUTDIR = config['outdir']

rule all:
    input:
        expand(OUTDIR+"/hapcut_input_f/mda_fragments.{CUTOFF}.txt", CUTOFF = CUTOFF)

rule create_adj_list:
    input:
        "{OUTDIR}/input/input.csv"
    output:
        "{OUTDIR}/wext_input/adj_list.txt"        
    shell:
        "python scripts/create_wext_input.py {input} {output}"
   
rule run_wext_preprocessing:
    input:
        "{OUTDIR}/wext_input/adj_list.txt"
    output:
        "{OUTDIR}/wext_input/adj_list.txt.data.json",
        "{OUTDIR}/wext_input/adj_list.txt.weights.npy",
    threads: 30
    shell: 
        "bash scripts/run_wext_preprocessing.bash {input} {threads}"

rule enumerate_cases:
    input:
        data="{OUTDIR}/wext_input/adj_list.txt.data.json",
        weights="{OUTDIR}/wext_input/adj_list.txt.weights.npy",
    output:
        cases="{OUTDIR}/caselist/cases.txt"
    run:
        # Create list of pairs
        # run parallel with available threads
        import numpy as np
        import json         
        
        data=input.data
        with open(data[0]) as data_file:    
            data = json.load(data_file)
        num_genes = len(data[u'genes'])
        print("Number of genes:{}".format(num_genes))
        with open(output.cases[0], 'w') as temp: 
            for v in range(1, (num_genes + 100)//100):
                w = v*100
                temp.write(str(w)+'\n')

rule run_wext:
    input:
        data="{OUTDIR}/wext_input/adj_list.txt.data.json",
        weights="{OUTDIR}/wext_input/adj_list.txt.weights.npy",
        cases="{OUTDIR}/caselist/cases.txt"
    output:
        file="{OUTDIR}/wext_output/total.txt"
    params:
        outdir="{OUTDIR}/wext_output/"
    threads: 30
    log:
        "log/wext.txt"
    shell:
        "parallel -j {threads} --joblog {log} python scripts/wext_calc.py {{1}} 100 {input.data} {input.weights} {params.outdir}  ::: `cat {input.cases}`; "
        "cat {params.outdir}/WextResult_*.txt > {output.file}"

rule create_mda_fragments:
    input:
        wext_output="{OUTDIR}/wext_output/total.txt",
        input = "{OUTDIR}/input/input.csv"
    output:
        mda_fragments="{OUTDIR}/hapcut_input_f/mda_fragments.{CUTOFF}.txt",
    shell:
        "python scripts/create_hapcut_input_fishers.py {input.input} {input.wext_output} {output.mda_fragments} {wildcards.CUTOFF};"




