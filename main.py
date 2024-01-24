#from configparser import ConfigParser
import argparse
import yaml
import logging

#import models
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier
from sklearn.ensemble import HistGradientBoostingClassifier
from dl_algo.long_short_term_memory import *

#import dataloaders
from data_loader.HGTDB_data_loader import HGTDBDataLoader, HGTDBDatasetSequential
from data_loader.NCBI_data_loader import NCBIDataLoader

#import bases
from base.binary_classifier import BinaryClassifier
from base.binary_classifier_dl import BinaryClassifierDLSequential



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



def load_algorithm(algorithm_type, algorithm):
    if algorithm_type == 'c':
        # for classical algorithm
        
        if algorithm == 'LGBM':
            return LGBMClassifier()
        elif algorithm == 'XGBOOST':
            return XGBClassifier()
        elif algorithm == 'HGBC':
            raise HistGradientBoostingClassifier()
        else:
            raise ValueError(f'no such thing as {algorithm}')
    elif algorithm_type == 'd':
        # for deep learning
        # TODO: Params could be usefule here...
        # Models should be configurable easier way!
        
        if algorithm == 'LSTM_v3':
            # need to make this configurable for now leave it
            return LSTMHGTTagger_v3(4,100,1)
        elif algorithm == 'LSTM_v4':
            # need to make this configurable for now leave it
            return LSTMHGTTagger_v4(4,100,1)
        elif algorithm == 'LSTM_v5':
            # need to make this configurable for now leave it
            return LSTMHGTTagger_v5(4,100,1)
        else:
            raise ValueError(f'no such thing as {algorithm}')

def load_dataloader(dataloader, partition_file, data_type):
    if dataloader == 'HGTDB':
        return HGTDBDataLoader(data_type, partition_file)
    elif dataloader == 'NCBI':
        return NCBIDataLoader(partition_file, data_type)
    elif dataloader == 'sequential':
        # currently working only with partition_file/HGTDB_firmicutes_trisplit.csv
        # and data type only A!
        hgtdb_train = HGTDBDatasetSequential(data_type,partition_file, 'train')
        hgtdb_valid = HGTDBDatasetSequential(data_type,partition_file, 'valid')
        hgtdb_test = HGTDBDatasetSequential(data_type,partition_file, 'test')
        return hgtdb_train, hgtdb_valid, hgtdb_test
        
    
def load_type( model_type, name, algorithm_type,  algorithm, dataloader,params, mode='train'):
    if algorithm_type == 'c':
        # for classical algorithm
        
        return BinaryClassifier(name, algorithm_type, algorithm, dataloader,params, mode)
    elif algorithm_type == 'd':
        # for deep learning
        train,valid,test = dataloader
        # return BinaryClassifierDLSequential(name, algorithm_type, algorithm, params['loss'], params['optimizer'],params['learning_rate'], train,valid,test, mode)
        return BinaryClassifierDLSequential(name, algorithm_type, algorithm, params, train, valid, test, mode)
    
    
    
    
    
    

def assert_config(args):
    list_of_algorithm = [
        'LGBM',
        'XGBOOST',
        'HGBC',
        'LSTM_v3',
        'LSTM_v4',
        'LSTM_v5']
    list_of_algorithm_types = [
        'c',
        'd'
    ]
    list_of_model_types = [
        'binary_classifier'
    ]
    list_of_dataset = [
        'HGTDB',
        'NCBI',
        'sequential'
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
        print(f"\tData type selected: \n\t\t{configuration['Dataset']['data_type']} ")
        print(f"\tPartition file: \n\t\t{configuration['Dataset']['partition_file'].replace('partition_file/','')} ")
    
    # check model configs
    if configuration['Model']['type'] not in list_of_model_types:
        raise ValueError('Unknown model type!')
    else:
        print(f"\tModel type: \n\t\t{configuration['Model']['type']}")

    # check algo type
    if configuration['Model']['algorithm_type'] not in list_of_algorithm_types:
        raise ValueError('Unknown algorithm type! chooses either classsical (c) or deep (d)')
    else:
        print(f"\talgorithm type: \n\t\t{configuration['Model']['algorithm_type']}")
        
    # check algo availability
    if configuration['Model']['algorithm'] not in list_of_algorithm:
        raise ValueError(f"algorithm: {configuration['Model']['algorithm']} not found! \nList of supported algo {list_of_algorithm}")
    else:
        print(f"\talgorithm chosen: \n\t\t{configuration['Model']['algorithm']}")
    
    
    print('-*-'*20)
        
        
    
    return configuration, args
    
    
def main():
    configuration, args = assert_config(parse_args())
    
    algorithm = load_algorithm(configuration['Model']['algorithm_type'], configuration['Model']['algorithm'])
    dataloader = load_dataloader(configuration['Dataset']['data_loader'], configuration['Dataset']['partition_file'], configuration['Dataset']['data_type'])
    
    
    print(configuration['Model']['params'])
    
    if args.mode == 'train':
        test_model = load_type( configuration['Model']['type'], args.config, configuration['Model']['algorithm_type'], algorithm, dataloader, configuration['Model']['params'])
        if configuration['Model']['algorithm_type'] == 'c':
            test_model.model_train()
        elif configuration['Model']['algorithm_type'] == 'd':
            test_model.model_train(configuration['Model']['epochs'])
        test_model.model_eval()
        test_model.save_model()
    elif args.mode == 'eval':
        test_model = load_type( configuration['Model']['type'], args.config, configuration['Model']['algorithm_type'], algorithm, dataloader, configuration['Model']['params'], 'eval')
        test_model.model_eval()
    
    


if __name__ == '__main__':
    main()

    
    
    # print out filtered table as partition file
    print('Done')