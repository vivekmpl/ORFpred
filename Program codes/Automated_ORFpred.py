# ORFpred Automated

"This is program to to run ORFpred"

print "- -"*10
print "ORFpred- Molecular Parasitology Lab BSBE IITB" , u"\u00a9"
print "- -"*10

# import module

import os

# Define classifier 

default = "Comb"
thresholds = 0.5
classifier_type = "l2r_l1r_l2loss_svc"
c_value = 5
eps_value =  5

def train_model(Model_type):
    Datasets = ORFpred_GUI.read_datasets(raw_input("Positive example(CDS)"),raw_input("Negative example(nc RNA ORF)"),raw_input("Negative example(Random ORF)"))
    Model_class = ORFpred_GUI.training_model(Datasets, Model_type)
    Model_class.solver_type_s, Model_class.C_v, Model_class.eps_v = classifier_type, c_value, eps_value
    Model_class.make_model()
    save_q = raw_input("Do you want to save the model(y/n)")
    if save_q == "y":
        Model_class.save()


def predictor():
    model_folder = raw_input("Folder containing the trained models")
    models = ["%s/"%model_folder+"%s" % i for i in os.listdir(model_folder)]
    models.sort()
    ORF_file = raw_input("Fasta file of ORF for prediction")
    predictor = ORFpred_GUI.predictor(ORF_file,models[0],models[1],models[2])
    Usage =  int(raw_input("Use expected thresholds (1) or use default thresholds (2)"))
    predictor.pred(Usage)
    types = int(raw_input("Output as a fasta file (1) or CSV list (2)"))
    file_name = raw_input("Output file name")
    predictor.write_pred(types,file_name)
      
    
def get_weights(Model_type, csv_file):
    Datasets = ORFpred_GUI.read_datasets(raw_input("Positive example(CDS)"),raw_input("Negative example(nc RNA ORF)"),raw_input("Negative example(Random ORF)"))
    Model_class = ORFpred_GUI.training_model(Datasets, Model_type)
    Model_class.solver_type_s, Model_class.C_v, Model_class.eps_v = classifier_type, c_value, eps_value
    Model_class.make_model()
    Model_class.give_weights(csv_file)

Question1 = raw_input("Choose options\nTrain model (1), Predict translatable smORFs (2), Get weights on features (3), Simulate expected frequencies (4)")

if Question1 in ["1","2","3"]:
    from Modules import ORFpred_GUI
    #Define UTR length
    ORFpred_GUI.utr_len = input
    ORFpred_GUI.dtr_len = input
    ORFpred_GUI.Bases_upstream = 0
    if Question1 == "1":
        train_model("Comb")
    elif Question1 == "2":
        predictor()
    elif Question1 =="3":
        Question2 = raw_input("Weight analysis of Positional features (1), Compositional features (2)")
        if Question2 == "2":
            get_weights("comp",raw_input("Output file name(.csv)"))
        if Question2 == "1":
            get_weights("pos",raw_input("Output file name(.csv)"))

elif Question1 in ["4"]:
    from Simulate_expected_value import expected_freq
    model_folder = raw_input("Folder containing the trained models")
    run_number = int(raw_input("Number of runs"))
    CSV_file = raw_input("Output file in csv format")
    expected_freq.simulate_runs(run_number, model_folder, CSV_file)
        


#END
