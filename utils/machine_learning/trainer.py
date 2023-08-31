from utils.machine_learning.evals.report import evaluate_predictions, dummy_probs, plot_PR_curve, plot_ROC_curve
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from utils.machine_learning.data_loader import ml_load_species
from sklearn import metrics
import numpy as np

class MLClassifier:
    def __init__(self, model):
        self.model=model
        self.precision = None
        self.roc_auc = None
        self.accuracy = None
        
    def cv_train(self, cross_val,cross_val_scoring, X_train, Y_train):
        n_scores = cross_val_score(self.model, X_train, Y_train, scoring=cross_val_scoring, cv=cross_val, n_jobs=-1)
        print( cross_val_scoring + ': %.3f (%.3f)' % (n_scores.mean(), n_scores.std()))
        
    def train(self, X_train, Y_train, sample_weight=None):
        #print('Fitting Model')
        self.model.fit(X_train,Y_train, sample_weight=sample_weight)
        #print('Done')
        
    def get_accuracy(self, X_test, Y_test):
        if self.accuracy is None:
            predictions = self.model.predict(X_test)
            self.accuracy = metrics.accuracy_score(Y_test, predictions)
            print(f'acc: {self.accuracy}')
        else:
            print(f'acc: {self.accuracy}')
    
    
    def get_precision(self, X_test, Y_test):
        if self.precision is None:
            predictions = self.model.predict(X_test)
            self.precision = metrics.precision_score(Y_test, predictions)
            print(f'prec: {self.precision}')
        else:
            print(f'prec: {self.precision}')
        
        
    def get_roc_auc(self, X_test, Y_test):
        #predictions = self.model.predict(X_test)
        #self.roc_auc = metrics.roc_auc_score(Y_test, self.model.decision_function(X_test))
        
        if self.roc_auc is None:
            y_proba = self.model.predict_proba(X_test)[:, 1]
            self.roc_auc = metrics.roc_auc_score(Y_test, y_proba)
            print(f'roc_auc: {self.roc_auc}')
        else:
            print(f'roc_auc: {self.roc_auc}')
        
    
        
        
    def eval(self, X_test, Y_test):
        '''
        deprecate tbhis maybe
        
        '''
        # make a predictions for test set
        predictions = self.model.predict(X_test)
        # Evaluate predictions
        evaluate_predictions(Y_test, predictions)

    def plot(self, X_train, Y_train,  X_test,Y_test):
        # prediction probabilities
        predict_prob = self.model.predict_proba(X_test)
        pos_probs_test = predict_prob[:, 1]
        # dummy model
        pos_dummy_probs = dummy_probs(X_train, Y_train, X_test)

        # prediction probabilities for training set
        predict_prob = self.model.predict_proba(X_train)
        pos_probs_train = predict_prob[:, 1]

        # visualization
        #model_label= species+ " data type" + data_type
        plot_ROC_curve(Y_train, pos_probs_train, Y_test, pos_probs_test, pos_dummy_probs, 'ROC')
        plot_PR_curve(Y_train, pos_probs_train, Y_test, pos_probs_test, pos_dummy_probs, 'PR')



def ml_cross_validator(model, ml_load_species_folds_func):
    fold_1_clf = MLClassifier(model)
    fold_2_clf = MLClassifier(model)
    fold_3_clf = MLClassifier(model)
    fold_4_clf = MLClassifier(model)
    fold_5_clf = MLClassifier(model)
    fold_6_clf = MLClassifier(model)
    
    (X_fold1, y_fold1), (X_fold2, y_fold2), (X_fold3, y_fold3), (X_fold4, y_fold4), (X_fold5, y_fold5), (X_fold6, y_fold6) = ml_load_species_folds_func
    
    fold_1_clf.train(np.concatenate((X_fold2, X_fold3, X_fold4, X_fold5, X_fold6)), np.concatenate((y_fold2, y_fold3, y_fold4, y_fold5, y_fold6)))
    fold_2_clf.train(np.concatenate((X_fold1, X_fold3, X_fold4, X_fold5, X_fold6)), np.concatenate((y_fold1, y_fold3, y_fold4, y_fold5, y_fold6)))
    fold_3_clf.train(np.concatenate((X_fold1, X_fold2, X_fold4, X_fold5, X_fold6)), np.concatenate((y_fold1, y_fold2, y_fold4, y_fold5, y_fold6)))
    fold_4_clf.train(np.concatenate((X_fold1, X_fold2, X_fold3, X_fold5, X_fold6)), np.concatenate((y_fold1, y_fold2, y_fold3, y_fold5, y_fold6)))
    fold_5_clf.train(np.concatenate((X_fold1, X_fold2, X_fold3, X_fold4, X_fold6)), np.concatenate((y_fold1, y_fold2, y_fold3, y_fold4, y_fold6)))
    fold_6_clf.train(np.concatenate((X_fold1, X_fold2, X_fold3, X_fold4, X_fold5)), np.concatenate((y_fold1, y_fold2, y_fold3, y_fold4, y_fold5)))
    
    fold_1_clf.get_precision(X_fold1, y_fold1)
    fold_1_clf.get_roc_auc(X_fold1, y_fold1)
    fold_1_clf.get_accuracy(X_fold1, y_fold1)
    fold_2_clf.get_precision(X_fold2, y_fold2)
    fold_2_clf.get_roc_auc(X_fold2, y_fold2)
    fold_2_clf.get_accuracy(X_fold2, y_fold2)
    fold_3_clf.get_precision(X_fold3, y_fold3)
    fold_3_clf.get_roc_auc(X_fold3, y_fold3)
    fold_3_clf.get_accuracy(X_fold3, y_fold3)
    fold_4_clf.get_precision(X_fold4, y_fold4)
    fold_4_clf.get_roc_auc(X_fold4, y_fold4)
    fold_4_clf.get_accuracy(X_fold4, y_fold4)
    fold_5_clf.get_precision(X_fold5, y_fold5)
    fold_5_clf.get_roc_auc(X_fold5, y_fold5)
    fold_5_clf.get_accuracy(X_fold5, y_fold5)
    fold_6_clf.get_precision(X_fold6, y_fold6)
    fold_6_clf.get_roc_auc(X_fold6, y_fold6)
    fold_6_clf.get_accuracy(X_fold6, y_fold6)
    
    total_precision = fold_1_clf.precision + fold_2_clf.precision + fold_3_clf.precision + fold_4_clf.precision + fold_5_clf.precision + fold_6_clf.precision
    total_roc_auc = fold_1_clf.roc_auc + fold_2_clf.roc_auc + fold_3_clf.roc_auc + fold_4_clf.roc_auc + fold_5_clf.roc_auc + fold_6_clf.roc_auc
    total_accuracy = fold_1_clf.accuracy + fold_2_clf.accuracy + fold_3_clf.accuracy + fold_4_clf.accuracy  + fold_5_clf.accuracy + fold_6_clf.accuracy
    
    print(f'Average precision : {total_precision/6}')
    print(f'Average roc_auc : {total_roc_auc/6}')
    print(f'Average accuracy : {total_accuracy/6}')
    

'''
Below should be deprecated but for now leave it

'''
def ml_trainer(model, species, data_type, cross_val, cross_val_scoring='roc_auc', plot=False):
    
    print(f"**** Training on {species} with type {data_type} ****\n")
    # init data
    X, Y = ml_load_species(species, data_type)

    
    # partition (they all are partioned the same way!)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.20, shuffle=False)
    n_scores = cross_val_score(model, X_train, Y_train, scoring=cross_val_scoring, cv=cross_val, n_jobs=-1)
    print( cross_val_scoring + ': %.3f (%.3f)' % (n_scores.mean(), n_scores.std()))
    
    # fit the model on the whole dataset
    model.fit(X_train, Y_train)

    # make a predictions for test set
    predictions = model.predict(X_test)
    # Evaluate predictions
    evaluate_predictions(Y_test, predictions)
    
    if plot:
        # prediction probabilities
        predict_prob = model.predict_proba(X_test)
        pos_probs_test = predict_prob[:, 1]
        # dummy model
        pos_dummy_probs = dummy_probs(X_train, Y_train, X_test)

        # prediction probabilities for training set
        predict_prob = model.predict_proba(X_train)
        pos_probs_train = predict_prob[:, 1]

        # visualization
        model_label= species+ " data type" + data_type
        plot_ROC_curve(Y_train, pos_probs_train, Y_test, pos_probs_test, pos_dummy_probs, model_label)
        plot_PR_curve(Y_train, pos_probs_train, Y_test, pos_probs_test, pos_dummy_probs, model_label)