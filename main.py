#from configparser import ConfigParser
import argparse
import yaml

#import models
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier

#import dataloaders
from data_loader.HGTDB_data_loader import HGTDBDataLoader
from data_loader.NCBI_data_loader import NCBIDataLoader

#import bases
from base.binary_classifier import BinaryClassifier

ROOT_CONFIGURATION_FOLDER = 'configuration/'
ROOT_PARTITION_FOLDER = 'partition_file/'

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c',
        '--config',
        type = str,
        default = 'test',
        help = 'Choose which config file'
    )
    parser.add_argument(
        '-m',
        '--mode',
        type = str,
        default = 'train',
        help = 'Choose mode: train or eval'
    )
    return parser.parse_args()



def load_algorithm(algorithm):
    if algorithm == 'LGBM':
        return LGBMClassifier()
    elif algorithm == 'XGBOOST':
        return XGBClassifier()
    elif algorithm == 'HGB':
        raise NotImplementedError('Not yet implemented')

def load_dataloader(dataloader, partition_file, data_type):
    if dataloader == 'HGTDB':
        return HGTDBDataLoader(data_type, partition_file)
    elif dataloader == 'NCBI':
        return NCBIDataLoader(partition_file, data_type)
        
    
def load_type( model_type, name, algorithm_type,  algorithm, dataloader, mode='train'):
    if model_type == 'binary_classifier':
        return BinaryClassifier(name, algorithm_type, algorithm, dataloader, mode)
    
    
    
    
    
    

def assert_config(args):
    list_of_algorithm = [
        'LGBM',
        'XGBOOST',
        'HGB']
    list_of_algorithm_types = [
        'c',
        'd'
    ]
    list_of_model_types = [
        'binary_classifier'
    ]
    list_of_dataset = [
        'HGTDB',
        'NCBI'
    ]
    list_of_data_type = ['A','B','C','D','E','F']
    list_of_mode = ['train', 'eval']
    
    print('-*-'*20)
    ## check args
    if args.mode not in list_of_mode:
        raise ValueError('Unknown mode chosen! choose train or load!')
    else:
        print(f'Selected mode: \n\t{args.mode}')
    
    # load yaml file
    configuration_file = ROOT_CONFIGURATION_FOLDER+args.config+'.yaml'
    print(f'Configuration file: \n\t{configuration_file}')
    try:
        with open(configuration_file, "r") as yamlfile:
            configuration = yaml.load(yamlfile, Loader=yaml.FullLoader)
    except FileNotFoundError:
        print(f"File {configuration_file} not found!")
        
    ## check configs
    # check dataset configs
    if configuration['Dataset']['data_loader'] not in list_of_dataset:
        raise ValueError('Unknwon Dataset!')
    else:
        print(f"\tDataset loaded: \n\t\t{configuration['Dataset']['data_loader']} ")
        
    if configuration['Dataset']['data_type'] not in list_of_data_type:
        raise ValueError('Unknwon Datatype!')
    else:
        print(f"\tData type selected \n\t\t{configuration['Dataset']['data_type']} ")
    
    # check model configs
    if configuration['Model']['type'] not in list_of_model_types:
        raise ValueError('Unknown model type!')
    else:
        print(f"\tModel type: \n\t\t{configuration['Model']['type']}")

    if configuration['Model']['algorithm_type'] not in list_of_algorithm_types:
        raise ValueError('Unknown algorithm type! chooses either classsical (c) or deep (d)')
    else:
        print(f"\talgorithm type: \n\t\t{configuration['Model']['algorithm_type']}")
        
    if configuration['Model']['algorithm'] not in list_of_algorithm:
        raise ValueError(f"algorithm: {configuration['Model']['algorithm']} not found! \nList of supported algo {list_of_algorithm}")
    
    
    print('-*-'*20)
        
        
    
    return configuration, args
    
    
def main():
    configuration, args = assert_config(parse_args())
    
    algorithm = load_algorithm(configuration['Model']['algorithm'])
    dataloader = load_dataloader(configuration['Dataset']['data_loader'], configuration['Dataset']['partition_file'], configuration['Dataset']['data_type'])
    
    
    if args.mode == 'train':
        test_model = load_type( configuration['Model']['type'], args.config, configuration['Model']['algorithm_type'],algorithm, dataloader)
        test_model.model_train()
        test_model.model_eval()
        test_model.save_model()
    elif args.mode == 'eval':
        test_model = load_type( configuration['Model']['type'], args.config, configuration['Model']['algorithm_type'],algorithm, dataloader, 'eval')
        test_model.model_eval()
    
    


if __name__ == '__main__':
    main()

    
    
    # print out filtered table as partition file
    print('Done')