from pathlib import Path

rule dna2protein:
    input:
        fasta=get_fasta_genes,
    output:
        faa="results/{date}/translate/{sample}.faa",
    log:
        "logs/{date}/analysis/dna2prot/{sample}.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/dna2protein.py"


if config["card"]["data"]["use_local"]:

    rule copy_local_CARD_data:
        output:
            json=get_card_db_file(),
        params:
            local_folder=lambda wildcards, output: Path(config["card"]["data"]["local_path"]).parent,
            root=get_root(),
            folder=lambda wildcards, output: Path(output.json).parent,
            tar_name=get_card_tar_file(),
        log:
            "logs/CARD_data_local_copy.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.folder}/ && "
            "(cd {params.local_folder}/ && "
            "tar cpfz - {params.tar_name} | (cd {params.root}/{params.folder}/ ; tar xpfz -)) && "
            "(cd {params.folder}/ && tar -xvf {params.tar_name})) > {log} 2>&1"
'''
else:

    rule download_CARD_data:
        output:
            model=get_card_db_folder(),
        params:
            download=config["card"]["data"]["url"],
            folder=get_card_db_folder(),
        log:
            "logs/CARD_data_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"
'''
#rgi load -i resources/card_db/card.json --local
#rgi main -i results/test_115/analysis/prokka/ab_115/ab_115.faa -o results/test_115/analysis/rgi_115.out -t protein -a DIAMOND -d wgs --clean --local

rule CARD_load_DB:
    input:
        get_card_db_file(),
    output:
        touch("logs/CARD_load_DB.done"),
    log:
        "logs/CARD_load_DB.log",
    conda:
        "../envs/card.yaml"
    shell:
        "rgi load -i {input} > {log} 2>&1"# --local


rule CARD_run:
    input:
        faa=rules.dna2protein.output.faa,
        db=rules.CARD_load_DB.output,
    output:
        txt="results/{date}/analysis/card/{sample}/{sample}.txt",
    params:
        path_wo_ext=lambda wildcards, output: Path(output.txt).with_suffix('')
    log:
        "logs/{date}/analysis/card/{sample}.log",
    threads: 12
    conda:
        "../envs/card.yaml"
    shell:
        "rgi main -i {input.faa} -o {params.path_wo_ext} -t protein "
        "-a DIAMOND -d wgs --clean > {log} 2>&1"# --local

"""
rule clone_plm_arg:
    output:
        main=get_plm_arg_main(),
        dummy=touch("logs/clone_plm_arg.done"),
    log:
        "logs/plm_arg_clone.log",
    params:
        git_url=config["plm_arg"]["git"],
        folder=lambda wildcards, output: Path(output.main).parent,
    conda:
        "../envs/unix.yaml"
    shell:
        "(rm -r -f {params.folder} && "
        "git clone --recursive {params.git_url} {params.folder}) > {log} 2>&1"


if config["plm_arg"]["model"]["use_local"]:

    rule copy_local_plm_arg_model:
        input:
            rules.clone_plm_arg.output.dummy,
        output:
            model=get_plm_arg_model_file(),
        params:
            local=config["plm_arg"]["model"]["local_path"],
            folder=get_plm_arg_model_path(),
        log:
            "logs/plm_arg_model_local_copy.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "cp {params.local} {params.folder}/ > {log} 2>&1"

else:

    rule download_plm_arg_model:
        input:
            rules.clone_plm_arg.output.dummy,
        output:
            model=get_plm_arg_model_file(),
        params:
            download=config["plm_arg"]["model"]["url"],
            folder=get_plm_arg_model_path(),
        log:
            "logs/plm_arg_model_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


if config["plm_arg"]["regression"]["use_local"]:

    rule copy_local_plm_arg_regression:
        input:
            rules.clone_plm_arg.output.dummy,
        output:
            reg=get_plm_arg_regression_file(),
        params:
            local=config["plm_arg"]["regression"]["local_path"],
            folder=get_plm_arg_model_path(),
        log:
            "logs/plm_arg_regression_local_copy.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "cp {params.local} {params.folder}/ > {log} 2>&1"

else:

    rule download_plm_arg_regression:
        input:
            rules.clone_plm_arg.output.dummy,
        output:
            reg=get_plm_arg_regression_file(),
        params:
            download=config["plm_arg"]["regression"]["url"],
            folder=get_plm_arg_model_path(),
        log:
            "logs/plm_arg_regression_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


rule run_plm_arg:
    input:
        model=get_plm_arg_model_file(),
        reg=get_plm_arg_regression_file(),
        main=get_plm_arg_main(),
        fasta=rules.prokka.output.faa[0],
    output:
        tsv="results/{date}/analysis/plm_arg/{sample}/{sample}.tsv",
    params:
        folder=lambda wildcards, input: Path(input.main).parent,
        main= lambda wildcards, input: Path(input.main).name,
        root=get_root(),
    log:
        "logs/{date}/analysis/plm_arg/{sample}.log",
    threads: 10
    conda:
        "../envs/plm_arg.yaml"
    shell:
        "(cd {params.folder}/ && "
        "python {params.main} predict -i {params.root}/{input.fasta} "
        "-o {params.root}/{output.tsv}) > {log} 2>&1"


rule extract_plm_arg_results:
    input:
        tsv=rules.run_plm_arg.output.tsv,
    output:
        arg="results/{date}/analysis/plm_arg/{sample}/{sample}_arg.csv",
        non_arg="results/{date}/analysis/plm_arg/{sample}/{sample}_non_arg.csv",
    log:
        "logs/{date}/analysis/plm_arg/{sample}_extract.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_plm_arg.py"


rule plot_plm_arg_results:
    input:
        csv=expand("results/{{date}}/analysis/plm_arg/{sample}/{sample}_arg.csv", sample=get_samples()),
    output:
        html="results/{date}/report/arg_plot.html",
    log:
        "logs/{date}/report/arg_plot.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_resistance.py"


"""


