""" This is a set of program ensemble used in predicting
translatable ORFs Program is made using  distributions
from BioPython, Numpy, mlpy, random, matplotlib
etc.. and requires all of these dependencies installed to run
- - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - -
There are  3 classes in this package
class - -
    1) read_datasets
            datasets = read_datasets(CDS training, NORF training, random training)
    2) training_model
            model = training_model(datasets, model_type)
            model_type = string ("Combined", "Positional", "Compositional")
            training the model : model.make_model()
            saving the trained model : model.save()
            get weights: make.give_weights(csv_file)
            csv_file = string ("File name for the weights output")
    3)Predictor
            predictor = predictor(record, positional_model, compositional_model, combined_model)
            record = string (Fasta file of the smORFs for prediction)
            positional_model, compositional_model, combined_model = string (Models output of training_model class)
            
    """

# UTR and DTR lenghts

UTR = int(raw_input("Lenght of the 5' untranslated region which is concidered"))
DTR = int(raw_input("Lenght of the 3' untranslated region which is concidered"))


#importing dependencies

import re, numpy , csv, random, os
import csv , mlpy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import get_features_comp  , get_features_pos, get_features_pos_only_utr

# Assigning utr lenghts and dtr lengths in feature extraction programs
get_features_comp.utr,get_features_comp.dtr =  UTR, DTR
get_features_pos.utr,get_features_pos.dtr =  UTR, DTR
get_features_pos_only_utr.utr,get_features_pos_only_utr.dtr =  UTR, DTR

# Defining functions

def read_fasta(fastafile, label):
    return [[i, label] for i in SeqIO.parse(fastafile, "fasta")]

def give_len (fastafile):
    return len([i for i in SeqIO.parse(fastafile, "fasta")])

#define vetorize
def vectorize_set_value2(j):
    for n , i in enumerate(j):
        try:
            float(i)
        except:
            if i == 'A':
                j[n] = 1.0
            if i == 'T':
                j[n] = 2.0
            if i == 'G':
                j[n] = 3.0
            if i == 'C':
                j[n] = 4.0
            if i == 'TGA':
                j[n] = 1.0
            if i == 'TAG':
                j[n] = 2.0
            if i == 'TAA':
                j[n] = 3.0
            if i not in ['A','T','G','C','TGA','TAG','TAA']:
                j[n] = 0.0
    return j

def select_fet2(seq, clas, type_fet):
    if type_fet in ["Comp" , "composition", "c","C", "compositional", "comp"]:
        fet_lis = get_features_comp.get_fet_predictor(seq)
    #print len(feature_selected)
    elif type_fet in ["pos", "positional", "p", "P", "Pos"]:
        fet_lis = get_features_pos_only_utr.get_fet_predictor(seq) # 
    elif type_fet in ["U","u","positional_utr", "utr", "UTR"]:
        fet_lis = get_features_pos_only_utr.get_fet_predictor(seq)
    else:
        pass
    return vectorize_set_value2(fet_lis)+[clas]

def select_fet1(seq,type_fet):
    if type_fet in ["Comp" , "composition", "c","C", "compositional","comp"]:
        fet_lis = get_features_comp.get_fet_predictor(seq)
    #print len(feature_selected)
    elif type_fet in ["pos", "positional", "p", "P", "Pos"]:
        fet_lis = get_features_pos_only_utr.get_fet_predictor(seq)
    elif type_fet in ["U","u","positional_utr", "utr", "UTR"]:
        fet_lis = get_features_pos_only_utr.get_fet_predictor(seq)
    else:
        pass
    return vectorize_set_value2(fet_lis)

# - - - - - - - - - - -

# - - - - - - - - - - - 
# making classes

class read_datasets:	
    def __init__(self, cds, norf, rorf):
        random.seed(1)
        self.cds = cds
        self.norf = norf
        self.rorf = rorf
        self.combine_data = read_fasta(cds, "1") + read_fasta(norf, "0") + read_fasta(rorf, "0")
        self.random_sample = random.sample(range(give_len(cds)+give_len(norf)+give_len(rorf)), (give_len(cds)+give_len(norf)+give_len(rorf))/2 )
        self.len = give_len(cds)+give_len(norf)+give_len(rorf)
    def training(self):
        return [self.combine_data[-1]] + [self.combine_data[i] for i in self.random_sample]
    def testing(self):
        testing_sample = [i for i in range(self.len) if i not in self.random_sample]
        return [self.combine_data[i] for i in testing_sample]

class training_model:
    def __init__(self, datasets, model_type):
        self.model_type = model_type
        self.datasets = datasets        
    solver_type_s = "l1r_l2loss_svc" # model details
    C_v = 5
    eps_v = 0.001
    weight_d = {}
    logit1 = mlpy.LibLinear(solver_type= solver_type_s, C = C_v , eps = eps_v , weight = weight_d)
    logit2 = mlpy.LibLinear(solver_type= solver_type_s, C = C_v , eps = eps_v , weight = weight_d)
    logit3 = mlpy.LibLinear(solver_type= solver_type_s, C = C_v , eps = eps_v , weight = weight_d)
    fet_names = []
    def classifier_type(self, solver, C_val , eps_val, weight_dict):
        self.solver_type_s = solver
        self.C_v = C_val
        self.eps_v = eps_val
        self.weight_d = weight_dict
    def make_model(self):
        if self.model_type not in ["Combined", "comb","Comb", "overall", "Overall", "o", "O"]:
            training_fet = np.array([select_fet2(i[0], i[1] , self.model_type) for  i in self.datasets.training()])
            testing_fet = np.array([select_fet2(i[0],i[1], self.model_type) for i in self.datasets.testing()])
            self.logit1.learn(training_fet[:,:-1], training_fet[:,-1])
            test_pred = [self.logit1.pred(i[:-1]) for i in testing_fet]
            test_vals = [i[-1] for i in testing_fet]
            cds_test_pred = [self.logit1.pred(i[:-1]) for i in testing_fet if i[-1] == '1']
            cds_test_vals = [i[-1] for i in testing_fet if i[-1] == '1']
            ncds_test_pred = [self.logit1.pred(i[:-1]) for i in testing_fet if i[-1] == '0']
            ncds_test_vals = [i[-1] for i in testing_fet if i[-1] == '0']
            print "=-"*5+"Accuracies" +"=-"*5
            print "Overall accuracy of predictor : %s" %mlpy.accuracy(test_vals , test_pred)
            print "Accuracy of predicting CDS correctly : %s" %mlpy.accuracy(cds_test_vals , cds_test_pred)
            print "Accuracy of predicting NCDS correctly : %s" %mlpy.accuracy(ncds_test_vals , ncds_test_pred)
        else:
            training_fet_1 = np.array([select_fet2(i[0], i[1] , "u") for  i in self.datasets.training()])
            testing_fet_1 = np.array([select_fet2(i[0],i[1], "u") for i in self.datasets.testing()])
            self.logit1.learn(training_fet_1[:,:-1], training_fet_1[:,-1])
            training_fet_2 = np.array([select_fet2(i[0], i[1] , "c") for  i in self.datasets.training()])
            testing_fet_2 = np.array([select_fet2(i[0],i[1], "c") for i in self.datasets.testing()])
            self.logit2.learn(training_fet_2[:,:-1], training_fet_2[:,-1])
            k = len(training_fet_1)
            j = len(testing_fet_1)
            comb_training = [[self.logit1.pred_probability(training_fet_1[i][:-1])[1] , self.logit2.pred_probability(training_fet_2[i][:-1])[1] , training_fet_1[i][-1]] for i in range(k)] # training file
            comb_training = np.array(comb_training, dtype = float)
            self.logit3.learn(comb_training[:,:-1], comb_training[:,-1])
            comb_testing = [[self.logit1.pred_probability(testing_fet_1[i][:-1])[1] , self.logit2.pred_probability(testing_fet_2[i][:-1])[1] , testing_fet_1[i][-1]] for i in range(j)] # testing file
            comb_testing = np.array(comb_testing, dtype = float)
            test_pred = [self.logit3.pred(i[:-1]) for i in comb_testing]
            test_vals = [i[-1] for i in comb_testing]
            cds_test_pred = [self.logit3.pred(i[:-1]) for i in comb_testing if i[-1] == 1]
            cds_test_vals = [i[-1] for i in comb_testing if i[-1] == 1]
            ncds_test_pred = [self.logit3.pred(i[:-1]) for i in comb_testing if i[-1] == 0]
            ncds_test_vals = [i[-1] for i in comb_testing if i[-1] == 0]
            print "=-"*5+"Accuracies" +"=-"*5
            print "Overall accuracy of predictor : %s" %mlpy.accuracy(test_vals , test_pred)
            print "Accuracy of predicting CDS correctly : %s" %mlpy.accuracy(cds_test_vals , cds_test_pred)
            print "Accuracy of predicting NCDS correctly : %s" %mlpy.accuracy(ncds_test_vals , ncds_test_pred)
    def save(self):
        model_name = raw_input("Model name")
        while True:
            mydir = "%s/%s"%(os.getcwd(), model_name)
            try:
                os.makedirs(mydir)
                break
            except OSError, e:
                if e.errno != 17:
                        raise   
                pass
        if self.model_type not in ["Combined", "comb","Comb", "overall", "Overall"]:
            self.logit1.save_model(self.model_type)
        else:
            self.logit1.save_model("%s/1_Pos_utr_model_%s" % (model_name,self.model_type))
            self.logit2.save_model("%s/2_Comp_model_%s" % (model_name,self.model_type))
            self.logit3.save_model("%s/3_Overall_model_%s" % (model_name,self.model_type))
    def give_weights(self, csvfile):
        weighed_feature = []
        if self.model_type in ["Comp" , "composition", "c","C", "compositional","comp"]:
            self.fet_names = get_features_comp.give_names()
        elif self.model_type in ["pos", "positional", "p", "P", "Pos"]:
            self.fet_names = get_features_pos_only_utr.give_names()
        elif self.model_type in ["U","u","positional_utr", "utr", "UTR"]:
            self.fet_names = get_features_pos_only_utr.give_names()
        else :
            print "Weights not available for combined model use compositional or positional models first to get weights"
        if self.fet_names!= []:
            weights = self.logit1.w()
            for i in range(len(self.fet_names)):
                if weights[i] < 0.0:
                    weighed_feature.append([self.fet_names[i] , abs(weights[i]) , "+"])
                else:
                    weighed_feature.append([self.fet_names[i] , abs(weights[i]) , "-"])
            with open(csvfile, "w") as f:
                csvwriter = csv.writer(f, delimiter = ",")
                csvwriter.writerows(weighed_feature)
    
class predictor():
    def __init__(self, record, positional_model, compositional_model, combined_model):
        self.record = record
        self.pos_model = mlpy.LibLinear.load_model(positional_model)
        self.comp_model = mlpy.LibLinear.load_model(compositional_model)
        self.comb_model = mlpy.LibLinear.load_model(combined_model)
    expected_dict = {}
    UTR1 = UTR
    DTR1 = DTR
    expected_file = "%s/expected_values.csv"%os.path.dirname(os.path.abspath(__file__))
    with open(expected_file,"r") as f:
       csvreader = csv.reader(f, delimiter = "\t")
       for i in csvreader:
           expected_dict.update({int(i[0]):float(i[2])})
    record_type = "fasta_file"
    cds_seqs = []
    ncds_seqs = []
    cds_seqs_list  = [["ORF id","Description", "Positional score","Compositional score", "Overall score" ]]
    ncds_seqs_list = [["ORF id","Description", "Positional score","Compositional score", "Overall score" ]]
    def pred(self, Usage):
        if Usage == 1:
            if self.record_type == "fasta_file":
                for i in SeqIO.parse(self.record, "fasta"):
                    try:
                        seq_len = len(i.seq[self.UTR1:-self.DTR1])
                        listk = [self.pos_model.pred_probability(select_fet1(i, "u"))[1], self.comp_model.pred_probability(select_fet1(i, "c"))[1]]
                        comb_pred = self.comb_model.pred_probability(listk)[1]# the last braket value for old models its 0 , for new its 1.
                        if seq_len <= 300:
                            if comb_pred >= self.expected_dict[seq_len] and comb_pred >= 0.5 :
                                self.cds_seqs.append(SeqRecord(seq = i.seq, id = i.id , description = "%s Score - %s" % (i.description, comb_pred)))
                                self.cds_seqs_list.append([i.id , i.description , listk[0] , listk[1] , comb_pred])
                            else:
                                self.ncds_seqs.append(SeqRecord(seq = i.seq, id = i.id , description = "%s Score - %s" % (i.description, comb_pred)))
                                self.ncds_seqs_list.append([i.id , i.description , listk[0] , listk[1] , comb_pred])
                        else:
                            if comb_pred >= 0.5 :
                                self.cds_seqs.append(SeqRecord(seq = i.seq, id = i.id , description = "%s Score - %s" % (i.description, comb_pred)))
                                self.cds_seqs_list.append([i.id , i.description , listk[0] , listk[1] , comb_pred])
                            else:
                                self.ncds_seqs.append(SeqRecord(seq = i.seq, id = i.id , description = "%s Score - %s" % (i.description, comb_pred)))
                                self.ncds_seqs_list.append([i.id , i.description , listk[0] , listk[1] , comb_pred])
                    except:
                        print "error"
        elif Usage == 2:
            if self.record_type == "fasta_file":
                for i in SeqIO.parse(self.record, "fasta"):
                    try:
                        seq_len = len(i.seq[self.UTR1:-self.DTR1])
                        listk = [self.pos_model.pred_probability(select_fet1(i, "u"))[1], self.comp_model.pred_probability(select_fet1(i, "c"))[1]]
                        comb_pred = self.comb_model.pred_probability(listk)[1]# the last braket value for old models its 0 , for new its 1.
                        if comb_pred >= 0.5 :
                            self.cds_seqs.append(SeqRecord(seq = i.seq, id = i.id , description = "%s Score - %s" % (i.description, comb_pred)))
                            self.cds_seqs_list.append([i.id , i.description , listk[0] , listk[1] , comb_pred])
                        else:
                            self.ncds_seqs.append(SeqRecord(seq = i.seq, id = i.id , description = "%s Score - %s" % (i.description, comb_pred)))
                            self.ncds_seqs_list.append([i.id , i.description , listk[0] , listk[1] , comb_pred])
                    except:
                        pass
        else:
            print "Select valid options"
        print " Predictions done"
    def write_pred(self,types, file_name):
        print "Writing output to a file"
        if types == 1:
            SeqIO.write(self.cds_seqs, "%s-CDS"%file_name, "fasta")
            SeqIO.write(self.ncds_seqs, "%s-NCDS"%file_name, "fasta")
        elif types == 2:
            with open("%s-CDS.csv"%file_name, "w") as f:
                csvwriter = csv.writer(f, delimiter = ",")
                csvwriter.writerows(self.cds_seqs_list)
            with open("%s-NCDS.csv"%file_name, "w") as f:
                csvwriter = csv.writer(f, delimiter = ",")
                csvwriter.writerows(self.ncds_seqs_list)            

