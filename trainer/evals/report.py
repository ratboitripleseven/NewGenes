from turtle import color
from sklearn.dummy import DummyClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
import seaborn as sns

from sklearn.model_selection import train_test_split


def evaluate_predictions(Y_test, predictions):
    print(confusion_matrix(Y_test, predictions))
    print(classification_report(Y_test, predictions))

def dummy_probs(X_train, Y_train, X_test):
    # create and fit dummy model
    dummy_model = DummyClassifier(strategy='stratified')
    dummy_model.fit(X_train, Y_train)

    # calculate roc auc
    # dummy model
    yhat_dummy = dummy_model.predict_proba(X_test)
    pos_dummy_probs = yhat_dummy[:, 1]
    return(pos_dummy_probs)

def plot_ROC_curve(Y_train, pos_probs_train, Y_test, pos_probs_test, pos_dummy_probs, model_label, add_model_info=""):
    #markersize
    ms=6

    sns.set_style('whitegrid')
    sns.set_palette("Dark2")
    plt.figure(figsize=(15,15))

    # plot no skill roc curve
    plt.plot([0, 1], [0, 1], color="k", label='No skill') # linestyle='--'
    # calculate roc curve for training set for model
    fpr, tpr, _ = roc_curve(Y_train, pos_probs_train)
    # plot model roc curve
    plt.plot(fpr, tpr, marker='.',ms=ms, label="$\it{E. coli\ K12}$ training (%.3f)" % roc_auc_score(Y_train, pos_probs_train))#color="#ce1256",
    # calculate roc curve for test set for model
    fpr, tpr, _ = roc_curve(Y_test, pos_probs_test)
    # plot model roc curve
    plt.plot(fpr, tpr, marker='.',ms=ms, label="$\it{E. coli\ K12}$ test (%.3f)" % roc_auc_score(Y_test, pos_probs_test))#color="#007200"

    plt.axis(xmin=-0.01, ymin=0,ymax=1.01, xmax=1) #ymin=-0.001, xmin=-0.01,
    plt.xticks(size=24)

    # title
    plt.title('ROC curve {} {}'.format(model_label, add_model_info), fontsize=28)
    # axis labels
    plt.xlabel('False Positive Rate', fontsize=24)
    plt.ylabel('True Positive Rate', fontsize=24)
    # show the legend
    plt.legend(fontsize=24, loc="lower right")
    # show the plot
    plt.show()

def plot_PR_curve(Y_train, pos_probs_train, Y_test, pos_probs_test, pos_dummy_probs, model_label, add_model_info="", legend_pos=0):
    #markersize
    ms=2
    sns.set_style('whitegrid')
    sns.set_palette('Dark2')

    plt.figure(figsize=(10,10))

    # plot precision-recall curves
    # calculate the no skill line as the proportion of the positive class
    no_skill = len(Y_test[Y_test==1]) / len(Y_test)
    # plot the no skill precision-recall curve
    plt.plot([0, 1], [no_skill, no_skill], color="k", label='No skill (%.3f)' % average_precision_score(Y_test, pos_dummy_probs))#linewidth=3, color="#fe9929"
    # calculate & plot model precision-recall curve for training set
    precision, recall, _ = precision_recall_curve(Y_train, pos_probs_train)
    plt.plot(recall, precision, marker='.',ms=ms, label="$\it{E. coli\ K12}$ training (%.3f)" % average_precision_score(Y_train, pos_probs_train)) #color="#ce1256",
    # calculate & plot model precision-recall curve for test set
    precision, recall, _ = precision_recall_curve(Y_test, pos_probs_test)
    plt.plot(recall, precision, marker='.', ms=ms, label="$\it{E. coli\ K12}$ test (%.3f)" % average_precision_score(Y_test, pos_probs_test)) #color="#007200"
    # title
    plt.title('PR curve {} {}'.format(model_label, add_model_info), fontsize=16)
    # axis labels
    plt.axis(xmin=0, ymin=0.001)#, ymax= 1.005, xmax=1.005)
    plt.xlabel('Recall', fontsize=12)
    plt.ylabel('Precision', fontsize=12)
    # show the legend
    plt.legend(fontsize=12, loc=legend_pos)
    # show the plot
    plt.show()