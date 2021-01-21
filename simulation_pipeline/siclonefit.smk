configfile: "sim.yaml"

#requires python 
seeds = [ i+1 for i in range(config["nseeds"])]


rule all:
    input:
        expand("sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_ternary.tsv",
                  seed=seeds, 
                  nmuts=config["nmuts"], 
                  nreads=config["nreads"],
                  tdoublet=config["tdoublet"], 
                  alpha_fp=config["fp"], 
                  alpha_fn = config["fn"],
                  ado=config["ado"], 
                  bprec=config['bprec'], 
                  loss= config["loss"], 
                  cna =config["cna"]),
#        expand("prediction/siclonefit/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_prediction",
#                  seed=seeds, 
#                  nmuts=config["nmuts"], 
#                  nreads=config["nreads"],
#                  tdoublet=config["tdoublet"], 
#                  alpha_fp=config["fp"], 
#                  alpha_fn = config["fn"],
#                  ado=config["ado"], 
#                  bprec=config['bprec'], 
#                  loss= config["loss"], 
#                  cna=config["cna"])


rule ternarize:
    input:
        dis_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_discrete.tsv",
    output:
        outfile = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_ternary.tsv",
    shell:
        'python siclonefit_input.py -i {input.dis_file} -o {output.outfile}'

rule siclonefit:
    input:
        inputfile = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_ternary.tsv",
    output: directory("prediction/siclonefit/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_prediction")
    benchmark: "prediction/siclonefit/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_benchmark.log" 
    params:
        ncells = config["ncells"],
    log: 
        std = "prediction/siclonefit/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}.log",
        err = "prediction/siclonefit/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}.err.log" 
    shell:    
        "java -jar SiCloneFiTComplete.jar -n {wilcards.nmuts} -m {params.ncells} --ipMat {input.inputfile} -outDir {output} -df 1  "
        " > {log.std} 2> {log.err}" 
