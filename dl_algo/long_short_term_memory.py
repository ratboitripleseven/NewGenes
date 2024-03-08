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
        

class LSTMHGTTagger_v5(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMHGTTagger_v5, self).__init__()
        
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.act_func_1 = nn.ReLU()
        self.act_func_2 = nn.ReLU()
        self.act_func_3 = nn.ReLU()
        self.act_func_last = nn.Sigmoid()
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
        #print(f'1 {lstm_out.shape()}')
        
        # unpack
        # apparently this unpacks them? https://gist.github.com/HarshTrivedi/f4e7293e941b17d19058f6fb90ab0fec
        output, input_sizes = rnn.pad_packed_sequence(lstm_out, batch_first=True)
        #print(f'2 {output.size()}')
        output = self.act_func_1(self.hidden2hidden_1(output))
        output = self.act_func_2(self.hidden2hidden_2(output))
        output = self.act_func_3(self.hidden2hidden_3(output))
        output = self.act_func_last(self.hidden2tag(output))
        #print(f'3 {output.size()}')
        
        # tag_scores = F.log_softmax(tag_space, dim=1)
        # tag_scores = F.softmax(tag_space, dim=1)
        return output, input_sizes
    
class LSTMHGTTagger_v6(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMHGTTagger_v6, self).__init__()
        
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.act_func_last = nn.Sigmoid()
        self.hidden2tag = nn.Linear(hidden_dim, tagset_size)

    def forward(self, input):
        
        data, seq_length = input
        
        # pack padded sequence.. exp from work
        input = rnn.pack_padded_sequence(data, lengths=seq_length, batch_first=True, enforce_sorted=False)
        
        #input to model
        lstm_out, _ = self.lstm(input)
        #print(f'1 {lstm_out.size()}')
        
        # unpack
        # apparently this unpacks them? https://gist.github.com/HarshTrivedi/f4e7293e941b17d19058f6fb90ab0fec
        output, input_sizes = rnn.pad_packed_sequence(lstm_out, batch_first=True)
        #print(f'2 {output.size()}')
        output = self.act_func_last(self.hidden2tag(output))

        return output, input_sizes
    
class LSTMHGTTagger_unpadded_v6(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(LSTMHGTTagger_unpadded_v6, self).__init__()
        
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.act_func_last = nn.Sigmoid()
        self.hidden2tag = nn.Linear(hidden_dim, tagset_size)

    def forward(self, input):
        
        
        #input to model
        lstm_out,_ = self.lstm(input)
        #print(lstm_out)
        #print(lstm_out)
        output = self.act_func_last(self.hidden2tag(lstm_out))

        return output
    
class BiLSTMHGTTagger_v6(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(BiLSTMHGTTagger_v6, self).__init__()
        
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True, bidirectional = True)
        self.act_func_last = nn.Sigmoid()
        self.hidden2tag = nn.Linear(hidden_dim*2, tagset_size)

    def forward(self, input):
        
        data, seq_length = input
        
        # pack padded sequence.. exp from work
        input = rnn.pack_padded_sequence(data, lengths=seq_length, batch_first=True, enforce_sorted=False)
        
        #input to model
        lstm_out, _ = self.lstm(input)
        
        
        # unpack
        # apparently this unpacks them? https://gist.github.com/HarshTrivedi/f4e7293e941b17d19058f6fb90ab0fec
        output, input_sizes = rnn.pad_packed_sequence(lstm_out, batch_first=True)
        #print(f'2 {output.size()}')
        output = self.act_func_last(self.hidden2tag(output))

        return output, input_sizes
    
class LSTMHGTTagger_nofc_v1(nn.Module):

    def __init__(self, embedding_dim, hidden_dim):
        super(LSTMHGTTagger_nofc_v1, self).__init__()
        print('WARNING hidden dim shoul only be 1!')
        if hidden_dim != 1:
            raise Exception("Hidden dime can only be 1!!!!")
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.act_func_last = nn.Sigmoid()
        #self.hidden2tag = nn.Linear(hidden_dim*2, tagset_size)

    def forward(self, input):
        
        data, seq_length = input
        
        # pack padded sequence.. exp from work
        input = rnn.pack_padded_sequence(data, lengths=seq_length, batch_first=True, enforce_sorted=False)
        
        #input to model
        lstm_out, _ = self.lstm(input)
        
        
        # unpack
        # apparently this unpacks them? https://gist.github.com/HarshTrivedi/f4e7293e941b17d19058f6fb90ab0fec
        output, input_sizes = rnn.pad_packed_sequence(lstm_out, batch_first=True)
        #print(f'2 {output.size()}')
        output = self.act_func_last(output)

        return output, input_sizes
    
class BILSTMHGTTagger(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, tagset_size):
        super(BILSTMHGTTagger, self).__init__()
        
        self.last_epoch= 0
        self.hidden_dim = hidden_dim

        # The LSTM takes word embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True, bidirectional= True)
        self.dropout = nn.Dropout(0.25)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        self.hidden2hidden=nn.Linear(hidden_dim*2,hidden_dim)
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
        output = self.sigmoid(self.hidden2tag(output))
        
        # tag_scores = F.log_softmax(tag_space, dim=1)
        # tag_scores = F.softmax(tag_space, dim=1)
        return output, input_sizes