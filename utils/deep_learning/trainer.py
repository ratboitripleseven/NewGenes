import torch
from torch.utils.data import DataLoader
from sklearn.metrics import classification_report, confusion_matrix, precision_score, accuracy_score, roc_auc_score


class DLClassifier:
    def __init__(self, model, loss, optim, lr ):
        '''
        WARNING!: expected model is for BINARY classifier that has output a sigmoidal output function!!!!!!
        affects eval and training if otherwise!
        
        For now lr is required to initiate 
        need to sit on this first
        '''
        self.model = model
        self.loss = loss
        self.lr =lr
        self.optimizer = optim(self.model.parameters(), self.lr)
        
        self.roc_auc_score = None
        self.precision_score = None
        self.accuracy_score = None
        
    def train(self, ds, bs, epoch, report_loss=True):
        
        dl = DataLoader(ds, batch_size=bs, shuffle=True)
        epoch_number = 0

        epoch_loss = []
        for epoch in range(epoch):
            print('EPOCH {}:'.format(epoch_number + 1))
            acc_tr_loss = 0.
            # set model to train
            self.model.train()
            for (i, data) in enumerate(dl):
                
                inputs, labels = data['gc_signature'], data['hgt']
                
                self.optimizer.zero_grad()
                # Forward pass
                y_pred = self.model(inputs)
                # Compute Loss
                #train_loss = self.loss(y_pred.squeeze(), labels)
                try:
                    train_loss = self.loss(y_pred, torch.unsqueeze(labels, 1))
                except:
                    print(y_pred)
                    print(labels)

                # Backward pass
                train_loss.backward()
                
                # adjust weight
                self.optimizer.step()
                
                acc_tr_loss += train_loss.item()
                #running_loss += loss.item()
                #last_loss = running_loss / 1000 # loss per batch
                # print('  batch {} loss: {}'.format(i + 1, last_loss))
                #running_loss = 0.
            
            acc_tr_loss = acc_tr_loss / (len(dl)+1)
            # epoch_loss.append(acc_loss / len(dl))        
            # set model to not train
            self.model.train(False)
            if report_loss:
                print('LOSS train {}'.format(acc_tr_loss))
            epoch_number += 1
            
                
    def eval(self, ds, report_eval=True):
        self.model.eval()
        dl = DataLoader(ds, batch_size=1)
        
        preds = []
        labels = []
        #softmax_probas = []
        for (i, data) in enumerate(dl):
            input, label = data['gc_signature'], data['hgt']
            
            y_pred = self.model(input)
            #softmax = torch.nn.Softmax(dim=1)
            #softmax_out = softmax(y_pred)

            
            
            #preds.append(y_pred.round().detach().numpy())
            #labels.append(label.detach().numpy())
            # got it from work
            preds = [*preds, *y_pred.round().detach().numpy()]
            labels = [*labels,*label.detach().numpy()]
            #softmax_probas = [*softmax_probas, *softmax_out.detach().cpu().numpy()]
            
        
        self.accuracy_score = accuracy_score(labels, preds)
        self.precision_score = precision_score(labels, preds)
        self.roc_auc_score = roc_auc_score(labels, preds)
        
        if report_eval:
            cm = confusion_matrix(labels, preds)
            print(cm)
            report = classification_report(labels,preds, target_names=["Non-HGT","HGT"])
            print(report)
        
        
def dl_cross_validatior(model, ds_folds,bs_per_fold, epoch_per_fold):
    total_acc = 0
    total_prec = 0
    total_roc_auc= 0 
    for (train_fold, test_fold) in ds_folds:
        fold_classifier = model
        fold_classifier.train(train_fold, bs_per_fold, epoch_per_fold, report_loss=False)
        fold_classifier.eval(test_fold, report_eval=False)
        total_acc += fold_classifier.accuracy_score
        total_prec += fold_classifier.precision_score
        total_roc_auc += fold_classifier.roc_auc_score
        
    print(f'Avg acc: {total_acc/6}')
    print(f'Avg prec: {total_prec/6}')
    print(f'Avg roc_auc: {total_roc_auc/6}')
        
    

    