import os


configfile: "config/config.yaml"

## helper functions
def get_root():
    return os.getcwd()


def get_resource_path():
    return config["data_handling"]["resources"]


def get_run_date():
    return config["run_date"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fasta_genes(wildcards):
    return (pep.sample_table.loc[wildcards.sample]["fasta_genes"],)


def get_card_db_file():
    name = config["card"]["data"]["dbfile"]
    path = "{}CARD_db/{}".format(get_resource_path(), name)
    return path


def get_card_tar_file():
    if config["card"]["data"]["use_local"]:
        name = Path(config["card"]["data"]["local_path"]).name
    else:
        name = "card-data.tar.bz2"
    #path = "{}CARD_db/{}".format(get_resource_path(), name)
    return name


def get_plm_arg_main():
    script = config["plm_arg"]["main"]
    script_path = "{}{}".format(get_resource_path(), script)
    return script_path


def get_plm_arg_model_path():
    path = "{}PLM-ARG/models/".format(get_resource_path())
    return path


def get_plm_arg_model_file():
    name = (config["plm_arg"]["model"]["url"]).split("/")[-1]
    path = "{0}{1}".format(get_plm_arg_model_path(), name)
    return path


def get_plm_arg_regression_file():
    name = (config["plm_arg"]["regression"]["url"]).split("/")[-1]
    path = "{0}{1}".format(get_plm_arg_model_path(), name)
    return path

