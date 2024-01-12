import torch.nn as nn
#from torch.utils import data
from torch.nn.utils.rnn import pad_sequence
import torch
import torch.nn.utils.rnn as rnn
from torch import autograd
# for now only v3 works!



class LSTMHGTTagger(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMHGTTagger, self).__init__()
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        
        self.relu = nn.ReLU()
        self.hidden2hidden=nn.Linear(hidden_dim,hidden_dim)
        self.hidden2tag = nn.Linear(hidden_dim, tagset_size)

    def forward(self, input):
        
        data, seq_length = input
        
        # pack padded sequence.. exp from work
        input = rnn.pack_padded_sequence(data, lengths=seq_length, batch_first=True, enforce_sorted=False)
        
        #input to model
        lstm_out, _ = self.lstm(input)
        
        # unpack
        # apparently this unpacks them? https://gist.github.com/HarshTrivedi/f4e7293e941b17d19058f6fb90ab0fec
        output, input_sizes = rnn.pad_packed_sequence(lstm_out, batch_first=True)
        
        # to hidden, softmax and sigmoid!
        output = self.relu(self.hidden2hidden(output))
        output = torch.sigmoid(self.hidden2tag(output))
        
        # tag_scores = F.log_softmax(tag_space, dim=1)
        # tag_scores = F.softmax(tag_space, dim=1)
        return output, input_sizes
    
    
    
class LSTMHGTTagger_v2(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMHGTTagger_v2, self).__init__()
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.dropout = nn.Dropout(0.50)
        self.relu = nn.ReLU()
        self.hidden2hidden=nn.Linear(hidden_dim,hidden_dim)
        self.hidden2tag = nn.Linear(hidden_dim, tagset_size)

    def forward(self, input):
        
        data, seq_length = input
        
        # pack padded sequence.. exp from work
        input = rnn.pack_padded_sequence(data, lengths=seq_length, batch_first=True, enforce_sorted=False)
        
        #input to model
        lstm_out, _ = self.lstm(input)
        
        # unpack
        # apparently this unpacks them? https://gist.github.com/HarshTrivedi/f4e7293e941b17d19058f6fb90ab0fec
        output, input_sizes = rnn.pad_packed_sequence(lstm_out, batch_first=True)
        
        output = self.dropout(output)
        output = self.relu(self.hidden2hidden(output))
        output = torch.sigmoid(self.hidden2tag(output))
        
        # tag_scores = F.log_softmax(tag_space, dim=1)
        # tag_scores = F.softmax(tag_space, dim=1)
        return output, input_sizes
        #return tag_scores, input_sizes
        

        
class LSTMHGTTagger_v3(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMHGTTagger_v3, self).__init__()
        
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.dropout = nn.Dropout(0.25)
        # self.relu = nn.ReLU()
        self.act_func = nn.Sigmoid()
        self.hidden2hidden=nn.Linear(hidden_dim,hidden_dim)
        self.hidden2tag = nn.Linear(hidden_dim, tagset_size)

    def forward(self, input):
        
        data, seq_length = input
        
        # pack padded sequence.. exp from work
        input = rnn.pack_padded_sequence(data, lengths=seq_length, batch_first=True, enforce_sorted=False)
        
        #input to model
        lstm_out, _ = self.lstm(input)
        
        # unpack
        # apparently this unpacks them? https://gist.github.com/HarshTrivedi/f4e7293e941b17d19058f6fb90ab0fec
        output, input_sizes = rnn.pad_packed_sequence(lstm_out, batch_first=True)
        
        output = self.dropout(output)
        output = self.act_func(self.hidden2hidden(output))
        output = self.act_func(self.hidden2tag(output))
        
        # tag_scores = F.log_softmax(tag_space, dim=1)
        # tag_scores = F.softmax(tag_space, dim=1)
        return output, input_sizes
        #return tag_scores, input_sizes
        
        

class LSTMHGTTagger_v4(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMHGTTagger_v4, self).__init__()
        
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.dropout = nn.Dropout(0.25)
        # self.relu = nn.ReLU()
        self.act_func = nn.Sigmoid()
        self.hidden2hidden_1=nn.Linear(hidden_dim,hidden_dim)
        self.hidden2hidden_2=nn.Linear(hidden_dim,hidden_dim)
        self.hidden2hidden_3=nn.Linear(hidden_dim,hidden_dim)
        self.hidden2tag = nn.Linear(hidden_dim, tagset_size)

    def forward(self, input):
        
        data, seq_length = input
        
        # pack padded sequence.. exp from work
        input = rnn.pack_padded_sequence(data, lengths=seq_length, batch_first=True, enforce_sorted=False)
        
        #input to model
        lstm_out, _ = self.lstm(input)
        
        # unpack
        # apparently this unpacks them? https://gist.github.com/HarshTrivedi/f4e7293e941b17d19058f6fb90ab0fec
        output, input_sizes = rnn.pad_packed_sequence(lstm_out, batch_first=True)
        
        output = self.dropout(output)
        output = self.act_func(self.hidden2hidden_1(output))
        output = self.act_func(self.hidden2hidden_2(output))
        output = self.act_func(self.hidden2hidden_3(output))
        output = self.act_func(self.hidden2tag(output))
        
        # tag_scores = F.log_softmax(tag_space, dim=1)
        # tag_scores = F.softmax(tag_space, dim=1)
        return output, input_sizes
        
