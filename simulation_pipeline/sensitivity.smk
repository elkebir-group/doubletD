configfile: "sensitivity.yaml"

seeds = [ i+1 for i in range(config["nseeds"])]

def my_function(beta_prec, mult):
    return mult*beta_prec

rule all:
    input:
        expand( "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alph_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_ad.csv",
                  seed=seeds, 
                  nmuts=config["nmuts"], 
                  nreads=config["nreads"],
                  tdoublet=config["tdoublet"], 
                  alpha_fp=config["fp"], 
                  alph_fn = config["fn"],
                  ado=config["ado"], 
                  bprec=config['bprec'], 
                  loss= config["loss"], 
                  cna=config["cna"]),
        expand("prediction/doubletD/s{seed}/i{pdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_prediction.out",
                  seed=seeds, 
                  nmuts=config["nmuts"], 
                  nreads=config["nreads"],
                  tdoublet=config["tdoublet"], 
                  pdoublet=config["pdoublet"], 
                  alpha_fp=config["fp"], 
                  alpha_fn = config["fn"],
                  ado=config["ado"], 
                  bprec=config['bprec'], 
                  loss= config["loss"], 
                  cna=config["cna"]),
        expand("prediction/doubletD/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_b{pprec}_prediction.out",
                  seed=seeds, 
                  nmuts=config["nmuts"], 
                  nreads=config["nreads"],
                  tdoublet=config["tdoublet"], 
                  alpha_fp=config["fp"], 
                  alpha_fn = config["fn"],
                  ado=config["ado"], 
                  bprec=config['bprec'], 
                  pprec =  config['pprec'],
                  loss= config["loss"], 
                  cna=config["cna"]),
        
# rule simulate_tree:
#     output:
#         ad_file = "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_ad.csv",
#         dp_file = "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_dp.csv",
#         doublet_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_doublets.csv",
#         mutalleles_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_mutAlleles.csv",
#         refalleles_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_refAlleles.csv",
#         cellgenotype_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_cellgenotypes.csv",
#         binary_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_binary.csv",
#         tree_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_tree.csv",
#         param_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_param.csv",
#         prev_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_prev.csv",
#         loh_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_loh.csv",
#         cnv_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_cna.csv",
#         clusters_file = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_clusters.csv",
#     params:
#         nnodes = config["nnodes"],
#         ncells = config["ncells"],
#         hetero = config["hetero"],
#         rare = config['rare'],
#         nbdisp = config['nbdisp'],
#     log:
#         std = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}.log", 
#         err = "sim_aux/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}.err.log" 
#     shell:
#         "./scsim -tree {output.tree_file} -params {output.param_file}  "
#         " -mut-copies {output.mutalleles_file} -ref-copies {output.refalleles_file} -cell_genotypes {output.cellgenotype_file} -loh-file {output.loh_file} "
#         " -ad {output.ad_file} -dp {output.dp_file} -doublets {output.doublet_file} -scs {output.binary_file} -cnv-file {output.cnv_file} "
#         " -clusters {output.clusters_file} -prev {output.prev_file} -s {wildcards.seed} -delta {wildcards.tdoublet} -depth {wildcards.nreads} "
#         " -alpha-fp {wildcards.alpha_fp} -alpha-fn {wildcards.alpha_fn} -ado {wildcards.ado} -t {wildcards.nmuts} -n {params.nnodes} -exp {params.ncells} -het {params.hetero} -beta {wildcards.bprec} "
#         " -cna-rate {wildcards.cna} -loss-rate {wildcards.loss} -copy-error "
#         " -rare-thresh {params.rare} -nbdispersion {params.nbdisp} "
#         " > {log.std} 2> {log.err}"


rule detectDoublets:
    input:
        ad_file = "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_ad.csv",
        dp_file = "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_dp.csv",
    output: 
        outfile ="prediction/doubletD/s{seed}/i{pdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_prediction.out",
    benchmark: "prediction/doubletD/s{seed}/i{pdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_benchmark.log" 
    log: 
        std = "prediction/doubletD/s{seed}/i{pdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}.log",
        err = "prediction/doubletD/s{seed}/i{pdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}.err.log" 
    shell:
        "python ../scripts/doubletD.py --inputTotal {input.dp_file} --inputAlternate {input.ad_file}  "
        " --delta {wildcards.pdoublet} --beta {wildcards.ado} -o {output.outfile} "
        " --prec {wildcards.bprec} --alpha_fp {wildcards.alpha_fp} --alpha_fn {wildcards.alpha_fn} "
        " --estimate --cellcoal --missing > {log.std} 2> {log.err}"


rule testBetaPrec:
    input:
        ad_file = "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_ad.csv",
        dp_file = "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_dp.csv",
    output: 
        outfile ="prediction/doubletD/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_b{pprec}_prediction.out",
    benchmark: "prediction/doubletD/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_b{pprec}_benchmark.log" 
    log: 
        std = "prediction/doubletD/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_b{pprec}.log",
        err = "prediction/doubletD/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_b{pprec}.err.log" 
    shell:
        "python ../scripts/doubletD.py --inputTotal {input.dp_file} --inputAlternate {input.ad_file}  "
        " --delta {wildcards.tdoublet} --beta {wildcards.ado} -o {output.outfile} "
        " --prec {wildcards.pprec} --alpha_fp {wildcards.alpha_fp} --alpha_fn {wildcards.alpha_fn} "
        " --estimate --cellcoal --missing > {log.std} 2> {log.err}"

# rule scrublet:
#     input:
#         ad_file = "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_ad.csv",
#         dp_file = "read_counts/s{seed}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_dp.csv",
#     output: 
#         outfile ="prediction/scrublet/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_prediction.out",
#     benchmark: "prediction/doubletD/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}_benchmark.log" 
#     log: 
#         std = "prediction/scrublet/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}.log",
#         err = "prediction/scrublet/s{seed}/i{tdoublet}/m{nmuts}_r{nreads}_d{tdoublet}_fp{alpha_fp}_fn{alpha_fn}_a{ado}_p{bprec}_c{cna}_l{loss}.err.log" 
#     shell:    
#         "python scrublet_benchmark.py --inputAlternate {input.ad_file} --inputTotal {input.dp_file} "
#         " --delta {wildcards.tdoublet} -o {output.outfile} "
#         " > {log.std} 2> {log.err}"

