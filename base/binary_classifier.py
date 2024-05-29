import pandas as pd
import numpy as np
from data.utils.HGTDB.ml_data_loader import ml_load_species
from sklearn import metrics
import joblib
import os
import time
from datetime import datetime
import logging
import torch
import torch.nn.utils.rnn as rnn
from data_loader.utils.prep_genome import prep_genome

# this import need not be here for future
import torch.nn as nn
import torch.optim as optim

'''
TODO: Model loading
    - Think about when you have a model and only want to test on new data
    model loading and then testing meant that func eval needs to be 
    able to take in new partition file only for the test!
'''


class BinaryClassifier:
    def __init__(self, name, algotrithm_type, algorithm, dataloader,params=None, mode='train'):
        
        self.logger = None
        now = datetime.now()
        self.dt_string = now.strftime("%Y_%m_%d_%H_%M_%S")
        self.name = name
        self.mode = mode
        
        # create folder first to save
        self.save_folder_root = 'models/binary_classifier/'
        # create root save folder
        if not os.path.isdir(self.save_folder_root):
            os.makedirs(self.save_folder_root)
            print("creating folder : ", self.save_folder_root)
            print(f'Saving in {self.save_folder_root}')
        else:
            print(f'Saving in {self.save_folder_root}')
        
        self._init_logger()
        
        if self.mode == 'train':
            self.algorithm = algorithm
            # set params by kwargs
            if params is not None:
                self.algorithm.set_params(**params)
                print(self.algorithm)
        elif self.mode == 'eval':
            self.algorithm = self._load_model()
        elif self.mode == 'annotate':
            self.algorithm = self._load_model()
            self.annotate_folder_root = 'annotated_file_HGT/'
            if not os.path.isdir(self.annotate_folder_root):
                os.makedirs(self.annotate_folder_root)
                print("creating folder : ", self.annotate_folder_root)
                print(f'Saving in {self.annotate_folder_root}')
            else:
                print(f'Saving in {self.annotate_folder_root}')

        
        self.X_train = None
        self.Y_train = None
        self.X_valid = None
        self.Y_valid = None
        self.X_test = None
        self.Y_test = None
        self.accuracy = None
        self.precision = None
        self.recall = None
        self.f_one = None
        self.roc_auc = None
        
        
        if algotrithm_type not in ['c','d']:
            raise ValueError(f'Unknown model type {algotrithm_type}')
        else:
            self.algotrithm_type = algotrithm_type
            
        
        
        

        self.dataloader = dataloader
        self._set_data()
        if self.mode != 'annotate':
            # dataloader is non when annotating
            self._print_data_statistics()
          
        #self._dataset_prep('E')
        
    def _init_logger(self):
        # create log folder to save logging
        folder_name = os.path.join(self.save_folder_root,self.name)
        if not os.path.isdir(folder_name):
            os.makedirs(folder_name)
            print("creating folder : ", folder_name)
            
            
        # init logger
        #logging.basicConfig(level=logging.INFO, filename=self.save_folder_root+"logs/"+self.name,filemode="w")
        logging.basicConfig(level=logging.INFO, filename=os.path.join(folder_name,'log_'+'BinaryClassifier'),filemode="w")
        self.logger=logging.getLogger()
        self.logger.info('ACCESS:')
        self.logger.info(self.dt_string)
    
    def _print_data_statistics(self):
        print('\tDataset statistics:')
        print(f'\t\tLength of Training set: {len(self.X_train)}')
        print(f'\t\tLength of Valid set: {len(self.X_valid)}')
        print(f'\t\tLength of Test set: {len(self.X_test)}')
        self.logger.info('\tDataset statistics:')
        self.logger.info(f'\t\tLength of Training set: {len(self.X_train)}')
        self.logger.info(f'\t\tLength of Valid set: {len(self.X_valid)}')
        self.logger.info(f'\t\tLength of Test set: {len(self.X_test)}')
        
        
        
        
    def _set_data(self):
        # for HGTDB
        if self.mode == 'annotate':
            # Use to annotate
            print('Annotation mode')
        else:
            self.X_train, self.Y_train, self.X_valid, self.Y_valid, self.X_test, self.Y_test  = self.dataloader.dataset_prep()
        
    
    def model_train(self):
        self.logger.info('Training Starts')
        self.algorithm.fit(self.X_train,self.Y_train)
        print('Valid set')
        self._get_accuracy(self.X_valid, self.Y_valid)
        self._get_precision(self.X_valid, self.Y_valid)
        self._get_recall(self.X_valid, self.Y_valid)
        self._get_f1(self.X_valid, self.Y_valid)
        self._get_roc_auc(self.X_valid, self.Y_valid)
        
        print('*'*30)
        print('EVAL')
        print('ACC\tPREC\tREC\tf1\tROCAUC')
        print(f'{self.accuracy}\t{self.precision}\t{self.recall}\t{self.f_one}\t{self.roc_auc}')
        self.logger.info('ACC\tPREC\tREC\tf1\tROCAUC')
        self.logger.info(f'{self.accuracy}\t{self.precision}\t{self.recall}\t{self.f_one}\t{self.roc_auc}')
        

    def model_eval(self):
        self.logger.info('Evaluation Starts')
        self._get_accuracy(self.X_test, self.Y_test)
        self._get_precision(self.X_test, self.Y_test)
        self._get_recall(self.X_test, self.Y_test)
        self._get_f1(self.X_test, self.Y_test)
        self._get_roc_auc(self.X_test, self.Y_test)
        print(f'acc: {self.accuracy}')
        self.logger.info(f'acc: {self.accuracy}')
        print(f'prec: {self.precision}')
        self.logger.info(f'prec: {self.precision}')
        print(f'recall: {self.recall}')
        self.logger.info(f'recall: {self.recall}')
        print(f'f1: {self.f_one}')
        self.logger.info(f'f1: {self.f_one}')
        print(f'roc_auc: {self.roc_auc}')
        self.logger.info(f'roc_auc: {self.roc_auc}')
        
    def model_annotate(self, annotated_dataset):
        # TODO
        _,length_of_genomes,to_annotate,ori_path = annotated_dataset.dataset_prep()
        #predictions = self.algorithm.predict(to_annotate)
        
        
        prev_index = 0
        current_index = 0
        print(f'Length of predicitons {len(to_annotate)}')
        for idx, length in enumerate(length_of_genomes):
            # get the indices
            
            
            print(f'hgts for genome: {ori_path[idx]} cur_idx {prev_index} end index {prev_index+length} of length {length}')
            # slices them
            to_predict = to_annotate[int(prev_index):int(prev_index+length)]
            hgts = self.algorithm.predict(to_predict)
            #assert length == len(hgts), 'ERROR'
            self._output_annotation(ori_path[idx],hgts)
            
            
            #end iteration and store starting index for nex iteration
            prev_index= prev_index + length
        
        #raise NotImplementedError('Not Yet Implemented')
        
    def model_single_annotate(self, single_genome):
        '''
        this is used for the snakemake pipeline
        '''
        # TODO
        raise NotImplementedError('Not yet implemented!')
        
        
        
    def _get_accuracy(self,x_set,y_set): 
        predictions = self.algorithm.predict(x_set)
        self.accuracy = metrics.accuracy_score(y_set, predictions)
        #print(f'acc: {self.accuracy}')
        #self.logger.info(f'acc: {self.accuracy}')
            
    def _get_precision(self,x_set,y_set):
        predictions = self.algorithm.predict(x_set)
        self.precision = metrics.precision_score(y_set, predictions)
        #print(f'prec: {self.precision}')
        #self.logger.info(f'prec: {self.precision}')

            
    def _get_recall(self, x_set, y_set):
        predictions = self.algorithm.predict(x_set)
        self.recall = metrics.recall_score(y_set, predictions)
        #print(f'recall: {self.recall}')
        #self.logger.info(f'recall: {self.recall}')
        
    def _get_f1(self, x_set, y_set):
        predictions = self.algorithm.predict(x_set)
        self.f_one = metrics.f1_score(y_set, predictions)
        #print(f'f1: {self.f_one}')
        #self.logger.info(f'f1: {self.f_one}')

        
        
    def _get_roc_auc(self, x_set, y_set):
        #predictions = self.model.predict(X_test)
        #self.roc_auc = metrics.roc_auc_score(Y_test, self.model.decision_function(X_test))
        
        y_proba = self.algorithm.predict_proba(x_set)[:, 1]
        self.roc_auc = metrics.roc_auc_score(y_set, y_proba)
        #print(f'roc_auc: {self.roc_auc}')
        #self.logger.info(f'roc_auc: {self.roc_auc}')

            
    def save_model(self):
        '''
        creates folder 
        and save models and perhaps report?
        '''
        '''
        # create root save folder
        if not os.path.isdir(self.save_folder_root):
            os.makedirs(self.save_folder_root)
            print("creating folder : ", self.save_folder_root)
            
        ''' 
        # create root save folder
        if not os.path.isdir(self.save_folder_root):
            os.makedirs(self.save_folder_root)
            print("creating folder : ", self.save_folder_root)
            
        # create folder for current experiment
        folder_name = os.path.join(self.save_folder_root,self.name)
        if not os.path.isdir(folder_name):
            os.makedirs(folder_name)
            print("creating folder : ", folder_name)
            
        
        joblib.dump(self.algorithm, os.path.join(folder_name,'model_'+'BinaryClassifier'))
        
    def _load_model(self):
        folder_name = os.path.join(self.save_folder_root,self.name)
        print(f'loading from {folder_name}')
        # check if desired model exists
        if not os.path.isdir(folder_name):
            print(f'folder {folder_name} does not exist!')
            
        return joblib.load(os.path.join(folder_name,'model_'+'BinaryClassifier'))
    
    def _output_annotation(self, ori_path, hgts):
        folder_name = os.path.join(self.annotate_folder_root,self.name)
        if not os.path.isdir(folder_name):
            os.makedirs(folder_name)
            print("creating annotation folder : ", folder_name)
            
        genome_dict = prep_genome(ori_path)
        output_name = ori_path.split('/')[-1]
        output_name = output_name.replace(output_name.split('.')[-1],'fasta')
        output_name = folder_name+'/'+output_name
        f = open(output_name, 'w')
        for i,gene in enumerate(genome_dict):
            f.write(">" + gene+ " [HGT:" + str(int(hgts[i]))+"]"+ "\n" + genome_dict[gene]['sequence'] + "\n")
        f.close()

        
        