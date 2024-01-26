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

# this import need not be here for future
import torch.nn as nn
import torch.optim as optim






class BinaryClassifierDLBase:
    def __init__(self, name, algotrithm_type, algorithm, params, train_dataset, valid_dataset, test_dataset, mode='train'):
        # params is a dict! containg
        # - loss
        # - optimizer and learning rate
        
        # does not need algorithm_type
        # algorithm is basically model!
        now = datetime.now()
        self.dt_string = now.strftime("%Y_%m_%d_%H_%M_%S")
        self.name = name
        
        
        # create folder first to save
        self.save_folder_root = 'models/binary_classifier/'
        # create root save folder
        if not os.path.isdir(self.save_folder_root):
            os.makedirs(self.save_folder_root)
            print("creating folder : ", self.save_folder_root)
            print(f'Saving in {self.save_folder_root}')
        else:
            print(f'Saving in {self.save_folder_root}')
        
        self.epoch = 0
        self.algorithm = algorithm
        self.loss = None
        self.optimizer = None
        
        self._set_params(**params)
        # self.learning_rate = learning_rate    
        #self._set_dl_loss(loss)
        #self._set_dl_optimizer(optimizer)
        
        # load model and optim state
        if mode  == 'eval':
            self._load_model()
            # self.algorithm = self._load_model()
            
        

        self.train_dataset = train_dataset
        self.valid_dataset = valid_dataset
        self.test_dataset = test_dataset
        
        self.train_loader  = None
        self.valid_loader  = None
        self.test_loader  = None
        
        self.accuracy = None
        self.precision = None
        self.roc_auc = None
        
        if algotrithm_type not in ['c','d']:
            raise ValueError(f'Unknown model type {algotrithm_type}')
        else:
            self.algotrithm_type = algotrithm_type

        #self._init_dataloader()
        self._set_data()
        #self._print_data_statistics()
        
    def _set_data(self,  batch_size = 3):
        print(f'Batch size: {batch_size}')
        self.train_loader = torch.utils.data.DataLoader(dataset=self.train_dataset,batch_size=batch_size,shuffle=False)
        self.valid_loader = torch.utils.data.DataLoader(dataset=self.valid_dataset,batch_size=batch_size,shuffle=False)
        self.test_loader = torch.utils.data.DataLoader(dataset=self.test_dataset,batch_size=batch_size,shuffle=False)
        
    def _set_params(self, loss, optimizer, learning_rate):
        # params here are defined as 
        # - loss
        # - optimizer and its learning rate
        
        self._set_dl_loss(loss)
        self._set_dl_optimizer(optimizer, learning_rate)
    
    
    def _set_dl_loss(self, loss='BCE'):
        if loss == 'BCE':
            self.loss = nn.BCELoss()
        else:
            raise ValueError(f'Unknown loss: {loss}')
        
    def _set_dl_optimizer(self, optimizer='Adam', learning_rate=0.01):
        if optimizer == 'Adam':
            self.optimizer = optim.Adam(self.algorithm.parameters(), lr=learning_rate)
        else:
            raise ValueError(f'Unknown optimizer {optimizer}')
        
        
    
        
    def save_model(self):
        # how do is
        # create folder for current experiment
        folder_name = os.path.join(self.save_folder_root,self.name)
        if not os.path.isdir(folder_name):
            os.makedirs(folder_name)
            print("creating folder : ", folder_name)
            
        torch.save({
            'epoch': self.epoch,
            'model_state_dict': self.algorithm.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'loss': self.loss
            }, os.path.join(folder_name,'model_'+'BinaryClassifier'))
        
    def _load_model(self):
        folder_name = os.path.join(self.save_folder_root,self.name)
        print(f'loading from {folder_name}')
        # check if desired model exists
        if not os.path.isdir(folder_name):
            print(f'folder {folder_name} does not exist!')
            
        checkpoint = torch.load(os.path.join(folder_name,'model_'+'BinaryClassifier'))
        
        self.epoch = checkpoint['epoch']
        self.algorithm.load_state_dict(checkpoint['model_state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        
        
        
        
class BinaryClassifierDLSequential(BinaryClassifierDLBase):
    def __init__(self, name, algotrithm_type, algorithm, params, train_dataset, valid_dataset, test_dataset, mode='train'):
        print('INTIALIZING DL')
        super().__init__(name, algotrithm_type, algorithm, params, train_dataset, valid_dataset, test_dataset, mode)
        
    def model_train(self, epochs = 15):
        
        for epoch in range(epochs):
            self.epoch += 1
            print(f'training epoch:{self.epoch}')
            acc_loss = 0.
            acc_vl_loss = 0.
            y_true_list = []
            y_pred_list = []
            
            # train !
            for (idx, data) in enumerate(self.train_loader):
                # print(idx)
                self.algorithm.train()
                self.algorithm.zero_grad()
                
                
                inputs = (data['datum'], data['seq_length'])
                output,_ = self.algorithm(inputs)
                
                # bad code but it works for now
                # it basically packs padded and unpad them....
                targets = rnn.pack_padded_sequence(data['label'], lengths=data['seq_length'], batch_first=True, enforce_sorted=False)
                targets,_ = rnn.pad_packed_sequence(targets, batch_first=True)
                
                
                # calc loss
                loss = self.loss(output, targets)
                loss.backward()

                # optimizer step
                self.optimizer.step()
                acc_loss += loss.item()
                
            # eval
            for (idx, data) in enumerate(self.valid_loader):
                self.algorithm.eval()
                
                with torch.no_grad():
                    inputs = (data['datum'], data['seq_length'])
                    output, input_sizes = self.algorithm(inputs)
                    
                    targets = rnn.pack_padded_sequence(data['label'], lengths=data['seq_length'], batch_first=True, enforce_sorted=False)
                    targets,_ = rnn.pad_packed_sequence(targets, batch_first=True)
                    # calc loss
                    loss = self.loss(output, targets)

                    acc_vl_loss += loss.item()


                
                #targets = torch.squeeze(targets,2)
                # the squeeze here is needed and fine
                #output = torch.squeeze(output,2)
                
                # this one not sure! it takes away the batch!
                


                #for i in range(targets.size(0)):
                #    y_true_list = [*y_true_list, *targets[i].detach().numpy()]
                #    y_pred_list = [*y_pred_list, *output[i].detach().round().numpy()]

            #acc = metrics.accuracy_score(y_true_list, y_pred_list)
            #print(f'val acc: {acc}')
                
            print(f'acculumative loss: {acc_loss/len(self.train_loader)}')
            print(f'acculumative vl loss: {acc_vl_loss/len(self.valid_loader)}')

    
    def model_eval(self):
        #y_true_list = []
        #y_pred_list = []
        y_pred_true = []
        list_input_sizes = []
        
        
        print(len(self.test_loader))
        for (idx, data) in enumerate(self.test_loader):
            self.algorithm.eval()
            with torch.no_grad():
                inputs = (data['datum'], data['seq_length'])
                output, input_sizes = self.algorithm(inputs)
                
                targets = rnn.pack_padded_sequence(data['label'], lengths=data['seq_length'], batch_first=True, enforce_sorted=False)
                targets,_ = rnn.pad_packed_sequence(targets, batch_first=True)
                
                #accuracy = (output.round() == targets).float().mean()
                #print(accuracy)
                
            targets = torch.squeeze(targets,2)
            # the squeeze here is needed and fine
            output = torch.squeeze(output,2)
            
            # this one not sure! it takes away the batch!
            

            for i in range(targets.size(0)):
                #y_true_list = [*y_true_list, *targets[i].detach().numpy()]
                #print(f'{i} {output}')
                #y_pred_list = [*y_pred_list, *output[i].detach().round().numpy()]
                list_input_sizes.append(input_sizes[i])
                # collected are all padded data!
                y_pred_true.append((targets[i].detach().numpy().tolist(), output[i].detach().round().numpy().tolist()))
                
        for idx in range(len(y_pred_true)):
            print('***'*20)
            print(f'Test genome no {idx} of length: {list_input_sizes[idx]} ')
            slice_index = list_input_sizes[idx]
            targs, preds = y_pred_true[idx]
            # slice the padding!
            targs, preds = targs[:slice_index], preds[:slice_index]
            #print(len(targs))
            acc = metrics.accuracy_score(targs, preds)
            prec = metrics.precision_score(targs, preds)
            reca = metrics.recall_score(targs, preds)
            f1_score = metrics.f1_score(targs, preds)
            conf_mat = metrics.confusion_matrix(targs, preds)
            print(f'test acc: {acc}')
            print(f'test prec: {prec}')
            print(f'test recall: {reca}')
            print(f'test f1: {f1_score}')
            print('Confusion matrix:')
            print(conf_mat)
            
        
        #print('***'*20)
        #print(f'Trained for {self.epoch} epochs')
        #print(f'Support: {len(y_pred_list)}')
        #print(f'HGTs predicted: {y_pred_list.count(1)}')
        #print(f'HGTs actual: {y_true_list.count(1)}')
        #print('***'*20)
        # print('-*-'*20)
        #acc = metrics.accuracy_score(y_true_list, y_pred_list)
        #prec = metrics.precision_score(y_true_list, y_pred_list)
        #reca = metrics.recall_score(y_true_list, y_pred_list)
        #f1_score = metrics.f1_score(y_true_list, y_pred_list)
        #conf_mat = metrics.confusion_matrix(y_true_list, y_pred_list)
        #print(f'test acc: {acc}')
        #print(f'test prec: {prec}')
        #print(f'test recall: {reca}')
        #print(f'test f1: {f1_score}')
        #print('Confusion matrix:')
        #print(conf_mat)
