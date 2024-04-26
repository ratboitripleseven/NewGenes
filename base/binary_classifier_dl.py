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

# for annotation
from data_loader.utils.prep_genome import prep_genome
SEQUENCES_FOLDER = 'data/NCBI/sequence_files'




class BinaryClassifierDLBase:
    def __init__(self, name, algotrithm_type, algorithm, params, train_dataset, valid_dataset, test_dataset, mode='train'):
        # params is a dict! containg
        # - loss
        # - optimizer and learning rate
        
        # does not need algorithm_type
        # algorithm is basically model!
        self.logger = None
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
        
        self._init_logger()
        
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
            self.logger.info("Model Loaded") 
            self._load_model()
        elif mode  == 'annotate':
            self.logger.info("Model Loaded") 
            self._load_model()
            self.annotate_folder_root = 'annotated_file_HGT/'
            if not os.path.isdir(self.annotate_folder_root):
                os.makedirs(self.annotate_folder_root)
                print("creating folder : ", self.annotate_folder_root)
                print(f'Saving in {self.annotate_folder_root}')
            else:
                print(f'Saving in {self.annotate_folder_root}')
            
        

        self.train_dataset = train_dataset
        self.valid_dataset = valid_dataset
        self.test_dataset = test_dataset
        
        self.train_loader  = None
        self.valid_loader  = None
        self.test_loader  = None
        
        self.accuracy = None
        self.precision = None
        self.roc_auc = None
        
        # TODO: check this I think it is redundant
        if algotrithm_type not in ['c','d','e']:
            raise ValueError(f'Unknown model type {algotrithm_type}')
        else:
            self.algotrithm_type = algotrithm_type

        #self._init_dataloader()
        if mode != 'annotate':
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
    
    def output_annotation(self, list_of_genomes, list_of_hgts):
        '''
        Called when annotating
        Output fast file
        
        This need to be redone! this is so badly coded
        '''
        folder_name = os.path.join(self.annotate_folder_root,self.name)
        if not os.path.isdir(folder_name):
            os.makedirs(folder_name)
            print("creating annotation folder : ", folder_name)
            
        print(list_of_genomes)
        for idx, genome in enumerate(list_of_genomes):
            genome_dict = prep_genome(genome)
            
            output_name = genome.split('/')[-1]
            output_name = output_name.replace(output_name.split('.')[-1],'fasta')
            output_name = folder_name+'/'+output_name
            f = open(output_name, 'w')
            #print(genome_dict)
            for i,gene in enumerate(genome_dict):
                f.write(">" + gene+ " [HGT:" + str(int(list_of_hgts[idx][i]))+"]"+ "\n" + genome_dict[gene]['sequence'] + "\n")
            f.close()
            
            
        
        
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
        print('INTIALIZING DL SEQUENTIAL')
        super().__init__(name, algotrithm_type, algorithm, params, train_dataset, valid_dataset, test_dataset, mode)
        
    def model_train(self, epochs = 15):
        self.logger.info('START OF TRAINING')
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
            self.logger.info(f'epoch: {self.epoch}')
            self.logger.info(f'\tacculumative loss: {acc_loss/len(self.train_loader)}')
            self.logger.info(f'\tacculumative vl loss: {acc_vl_loss/len(self.valid_loader)}')

    
    def model_eval(self):
        #y_true_list = []
        #y_pred_list = []
        y_pred_true = []
        list_input_sizes = []
        y_pred_not_rounded = []
        
        self.logger.info('START OF EVAL')
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
                y_pred_not_rounded.append( output[i].detach().numpy().tolist())
                
        all_targs = []
        all_preds = []        
        
        for idx in range(len(y_pred_true)):
            print('***'*20)
            print(f'Test genome no {idx} of length: {list_input_sizes[idx]} ')
            slice_index = list_input_sizes[idx]
            targs, preds = y_pred_true[idx]
            preds_not_rounded =y_pred_not_rounded[idx]
            # slice the padding!
            targs, preds = targs[:slice_index], preds[:slice_index]
            preds_not_rounded = preds_not_rounded[:slice_index]
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
            print(f'preds:\n\t{preds_not_rounded}')
            print('Confusion matrix:')
            print(conf_mat)
            
            # aggregate all
            all_targs.extend(targs)
            all_preds.extend(preds)
            
            self.logger.info(f'test acc: {acc}')
            self.logger.info(f'test prec: {prec}')
            self.logger.info(f'test recall: {reca}')
            self.logger.info(f'test f1: {f1_score}')
            self.logger.info(f'preds:\n\t{preds_not_rounded}')
            self.logger.info('Confusion matrix:')
            self.logger.info(conf_mat)
            
        print('***'*20)
        print('Complete Evaluation')
        print('***'*20)
        total_acc = metrics.accuracy_score(all_targs, all_preds)
        total_prec = metrics.precision_score(all_targs, all_preds)
        total_reca = metrics.recall_score(all_targs, all_preds)
        total_f1_score = metrics.f1_score(all_targs, all_preds)
        total_conf_mat = metrics.confusion_matrix(all_targs, all_preds)
        
        print(f'test acc: {total_acc}')
        print(f'test prec: {total_prec}')
        print(f'test recall: {total_reca}')
        print(f'test f1: {total_f1_score}')
        print('Confusion matrix:')
        print(total_conf_mat)
        
        self.logger.info('***'*20)
        self.logger.info('Complete Evaluation')
        self.logger.info('***'*20)
        self.logger.info(f'test acc: {total_acc}')
        self.logger.info(f'test prec: {total_prec}')
        self.logger.info(f'test recall: {total_reca}')
        self.logger.info(f'test f1: {total_f1_score}')
        self.logger.info('Confusion matrix:')
        self.logger.info(conf_mat)
        
        
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
    
    def model_annotate(self, annotated_dataset):
        '''
        Similar to eval but no need to model_eval
        Instead of getting evaluation metrics, It annotates the genomes
        '''
        #init data loader for annotation
        annotated_loader = torch.utils.data.DataLoader(dataset=annotated_dataset,batch_size=1,shuffle=False)
        
        y_pred_rounded = []
        list_input_sizes = []
        genome_id = []
        
        y_pred_not_rounded = []
        
        self.logger.info('START OF ANNOTATION')
        print(len(annotated_loader))
        for (idx, data) in enumerate(annotated_loader):
            self.algorithm.eval()
            with torch.no_grad():
                inputs = (data['datum'], data['seq_length'])
                output, input_sizes = self.algorithm(inputs)
                genome_ids = data['data_id']
                
                
                #accuracy = (output.round() == targets).float().mean()
                #print(accuracy)
                
            # the squeeze here is needed and fine
            output = torch.squeeze(output,2)
            

            for i in range(output.size(0)):
                list_input_sizes.append(input_sizes[i])
                # collected are all padded data!
                y_pred_rounded.append( output[i].detach().round().numpy().tolist())
                y_pred_not_rounded.append( output[i].detach().numpy().tolist())
                genome_id.append(genome_ids[i])
                
        for idx in range(len(y_pred_rounded)):
            print('***'*20)
            print(f'Annotating genome:{genome_id[idx]} of length: {list_input_sizes[idx]} ')
            slice_index = list_input_sizes[idx]
            preds = y_pred_rounded[idx]
            preds_not_rounded =y_pred_not_rounded[idx]
            # slice the padding!
            preds_not_rounded = preds_not_rounded[:slice_index]
            #print(len(targs))

            
            self.logger.info('Confusion matrix:')
        self.output_annotation(genome_id,y_pred_rounded)



class BinaryClassifierDLNonSequential(BinaryClassifierDLBase):
    def __init__(self, name, algotrithm_type, algorithm, params, train_dataset, valid_dataset, test_dataset, mode='train'):
        print('INTIALIZING DL NON-SEQUENTIAL (UNPADDED!)')
        super().__init__(name, algotrithm_type, algorithm, params, train_dataset, valid_dataset, test_dataset, mode)
        
    def model_train(self, epochs = 15):
        self.logger.info('START OF TRAINING')
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
                
                
                inputs = data['datum']
                #print(inputs)
                output = self.algorithm(inputs)
                
                targets = data['label']
                
                
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
                    inputs = data['datum']
                    output = self.algorithm(inputs)
                    
                    targets = data['label']
                    # calc loss
                    loss = self.loss(output, targets)

                    acc_vl_loss += loss.item()


                
            #acc = metrics.accuracy_score(y_true_list, y_pred_list)
            #print(f'val acc: {acc}')
                
            print(f'acculumative loss: {acc_loss/len(self.train_loader)}')
            print(f'acculumative vl loss: {acc_vl_loss/len(self.valid_loader)}')
            self.logger.info(f'epoch: {self.epoch}')
            self.logger.info(f'\tacculumative loss: {acc_loss/len(self.train_loader)}')
            self.logger.info(f'\tacculumative vl loss: {acc_vl_loss/len(self.valid_loader)}')
            
    
    def model_eval(self):
        #y_true_list = []
        #y_pred_list = []
        targs = []
        preds = []
        preds_not_rounded = []
        
        self.logger.info('START OF EVAL')
        print(len(self.test_loader))
        for (idx, data) in enumerate(self.test_loader):
            self.algorithm.eval()
            with torch.no_grad():
                inputs = data['datum']
                output = self.algorithm(inputs)
                
                targets = data['label']

                
            #print(targets.size())
                #accuracy = (output.round() == targets).float().mean()
                #print(accuracy)
            targs = [*targs, *targets.detach().numpy().tolist()]
            #targs.append(targets.detach().numpy().tolist())
            preds = [*preds, *output.detach().round().numpy().tolist()]
            #preds.append(output.detach().round().numpy().tolist())
            preds_not_rounded = [*preds_not_rounded, *output.detach().numpy().tolist()]
            #preds_not_rounded.append(output.detach().numpy().tolist())

        

                
        print('***'*20)
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
        print(f'preds:\n\t{preds_not_rounded}')
        print('Confusion matrix:')
        print(conf_mat)
        
        self.logger.info(f'test acc: {acc}')
        self.logger.info(f'test prec: {prec}')
        self.logger.info(f'test recall: {reca}')
        self.logger.info(f'test f1: {f1_score}')
        self.logger.info(f'preds:\n\t{preds_not_rounded}')
        self.logger.info('Confusion matrix:')
        self.logger.info(conf_mat)