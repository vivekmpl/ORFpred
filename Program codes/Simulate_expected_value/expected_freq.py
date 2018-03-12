# A monte carlo like simulation to get expected frequency

from Bio import SeqIO
import  tempfile , random , os , re , csv , scipy.stats
import numpy as np
import numpy
import sequence_shuffler, orfpy_new_utr, ORFpred_expected_value

def poi_conf_int(poisson_mean, length):
   return poisson_mean + 1.96 * (numpy.sqrt(poisson_mean/length))

def simulate_expected(model_folder, temp_translatability_scores_file):
    # make random file
    models = ["%s/"%model_folder+"%s" % i for i in os.listdir(model_folder)]
    models.sort() # get the models in 
    seq_temp_file = tempfile.NamedTemporaryFile() # make a temp file to store sequences
    SeqIO.write(sequence_shuffler.shuffle_seq("%s/PlasmoDB-9.2_Pfalciparum3D7_Genome.fasta"%os.path.dirname(os.path.abspath(__file__))),seq_temp_file, "fasta") # write a  random sequence to that temp seq file
    print "Random genome made"
    orf_temp_list = []
    orf_temp_file = tempfile.NamedTemporaryFile() # make a temp file to store ORFs from random genome
    for i in SeqIO.parse(seq_temp_file.name,"fasta"):
        orf_temp_list.append(orfpy_new_utr.extract_orfs(i,9,300,350,100))# extract orfs
        orf_temp_list.append(orfpy_new_utr.extract_orfs_rev(i,9,300,350,100))# extract rev orfs
    seq_temp_file.close()
    for i in orf_temp_list:
        for k in i:
            orf_temp_file.file.write(k) # write ORFs from the list into the temp file
    print "ORFs isolated"
    print "%s  ORFs isolated" % np.sum([len(i) for i in orf_temp_list]) #1
    orf_temp_list = []
    predictor = ORFpred_expected_value.predictor(orf_temp_file.name, models[0],models[1], models[2]) # defined a predictor class
    predictor.give_scores(temp_translatability_scores_file.name)# Write scores and seq lenght to 
    orf_temp_file.close()
    print "prediction of scores done"

def simulate_runs(run_number, model_folder, CSV_file):
    temp_translatability_scores_file = tempfile.NamedTemporaryFile() # Make a temp score file
    # simulate runs
    for j in range(run_number):
       print "Started run number %s"%j
       simulate_expected(model_folder, temp_translatability_scores_file)
       print "Finished run number %s"%j
    #make segregated list
    segregated_list = [] # the final list
    print "started segregating"
    for i in range(9,300,3):
        len_temp_list_pcomp = []
        len_temp_list_ppos = []
        len_temp_list_pcds = []
        with open(temp_translatability_scores_file.name, "r") as f:
            csvreader = csv.reader(f, delimiter = ",")
            for j in csvreader:
                if float(j[-1]) == float(i):
                   len_temp_list_ppos.append(j[0])
                   len_temp_list_pcomp.append(j[1])
                   len_temp_list_pcds.append(j[2])
        len_temp_list_pcds = np.array(len_temp_list_pcds, dtype = float)
        len_temp_list_pcomp = np.array(len_temp_list_pcomp, dtype = float)
        len_temp_list_ppos = np.array(len_temp_list_ppos, dtype = float)
        m_ppos = np.mean(len_temp_list_ppos)
        m_pcomp = np.mean(len_temp_list_pcomp)
        m_pcds = np.mean(len_temp_list_pcds)
        v_ppos = np.var(len_temp_list_ppos)
        v_pcomp = np.var(len_temp_list_pcomp)
        v_pcds = np.var(len_temp_list_pcds)
        ci_ppos = scipy.stats.norm.interval(0.99 , m_ppos, v_ppos)[1]
        ci_pcomp = scipy.stats.norm.interval(0.99 , m_pcomp, v_pcomp)[1]
        ci_pcds = scipy.stats.norm.interval(0.99 , m_pcds, v_pcds)[1]
        #segregated_list.append([i, m_ppos, v_ppos , ci_ppos, m_pcomp, v_pcomp , ci_pcomp, m_pcds, v_pcds , ci_pcds])
        segregated_list.append([i, m_pcds, ci_pcds])
    temp_translatability_scores_file.close()
    print len(segregated_list)
    with open(CSV_file, "w") as f:
       csvwriter = csv.writer(f, delimiter = ",")
       csvwriter.writerows(segregated_list)

