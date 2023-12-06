import pandas as pd
import numpy as np
from data.utils.HGTDB.ml_data_loader import ml_load_species
from sklearn import metrics
import joblib
import os
import time
from datetime import datetime



'''
TODO: Model loading
    - Think about when you have a model and only want to test on new data
    model loading and then testing meant that func eval needs to be 
    able to take in new partition file only for the test!
'''


class BinaryClassifier:
    def __init__(self, name, algotrithm_type, algorithm, dataloader, mode='train'):
        
        now = datetime.now()
        self.dt_string = now.strftime("%Y_%m_%d_%H_%M_%S")
        self.name = name
        self.save_folder_root = 'models/binary_classifier/'
        
        if mode == 'train':
            self.algorithm = algorithm
        elif mode  == 'eval':
            self.algorithm = self._load_model()

        
        self.X_train = None
        self.Y_train = None
        self.X_test = None
        self.Y_test = None
        self.accuracy = None
        self.precision = None
        self.roc_auc = None
        
        
        if algotrithm_type not in ['c','d' ]:
            raise ValueError(f'Unknown model type {algotrithm_type}')
        else:
            self.algotrithm_type = algotrithm_type
            
        
        
        

        self.dataloader = dataloader
        self._set_data()
        self._print_data_statistics()
          
        #self._dataset_prep('E')
        
    
    def _print_data_statistics(self):
        print('\tDataset statistics:')
        print(f'\t\tLength of Training set: {len(self.X_train)}')
        print(f'\t\tLength of Test set: {len(self.X_test)}')
        
        
        
        
    def _set_data(self):
        # for HGTDB
        self.X_train,self.Y_train, self.X_test, self.Y_test  = self.dataloader.dataset_prep()
        
    
    def model_train(self):
        self.algorithm.fit(self.X_train,self.Y_train)

    def model_eval(self):
        self._get_accuracy()
        self._get_precision()
        self._get_roc_auc()
        
        
    def _get_accuracy(self): 
        if self.accuracy is None:
            predictions = self.algorithm.predict(self.X_test)
            self.accuracy = metrics.accuracy_score(self.Y_test, predictions)
            print(f'acc: {self.accuracy}')
        else:
            print(f'acc: {self.accuracy}')
            
    def _get_precision(self):
        if self.precision is None:
            predictions = self.algorithm.predict(self.X_test)
            self.precision = metrics.precision_score(self.Y_test, predictions)
            print(f'prec: {self.precision}')
        else:
            print(f'prec: {self.precision}')
        
        
    def _get_roc_auc(self):
        #predictions = self.model.predict(X_test)
        #self.roc_auc = metrics.roc_auc_score(Y_test, self.model.decision_function(X_test))
        
        if self.roc_auc is None:
            y_proba = self.algorithm.predict_proba(self.X_test)[:, 1]
            self.roc_auc = metrics.roc_auc_score(self.Y_test, y_proba)
            print(f'roc_auc: {self.roc_auc}')
        else:
            print(f'roc_auc: {self.roc_auc}')
            
    def save_model(self):
        '''
        creates folder 
        and save models and perhaps report?
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
            

    
    