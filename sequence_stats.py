import cPickle
from Bio.Seq import Seq
from collections import Counter
from models import *

def get_largest_empty_region (_features):
    start_empty_region = _features[1].location.end.position
    end_empty_region = _features[2].location.start.position

    for i in xrange(2,len(_features)-1):
        local_start = _features[i].location.end.position
        local_end = _features[i+1].location.start.position
        if (local_end-local_start) > (end_empty_region - start_empty_region):
            end_empty_region = local_end
            start_empty_region = local_start

    local_start = _features[i+1].location.end.position
    local_end = _features[0].location.end.position
    if (local_end-local_start) > (end_empty_region - start_empty_region):
        end_empty_region = local_end
        start_empty_region = local_start

    return (start_empty_region, end_empty_region)



def get_sequence_stats(genome):
    _features = cPickle.loads(genome.features)
    d_loop_start = 0
    d_loop_end = 0

    total_protein_coding_GC_count = 0
    total_protein_coding_seq_length = 0
    total_protein_coding_nucleotide_statistics = Counter()

    #get sequence as string and then as Seq() obj
    sequence_as_string = genome.fasta[genome.fasta.find('\n')+1:].replace('\n','')
    mt_dna = Seq(sequence_as_string)

    total_seq_coding_GC_count = sequence_as_string.count('GC')
    total_seq_coding_seq_length = len(sequence_as_string)

    counter_seq_coding_nucleotide_statistics =  Counter(sequence_as_string)
    total_seq_coding_nucleotide_statistics = {'T':counter_seq_coding_nucleotide_statistics['T'],
                                 'A':counter_seq_coding_nucleotide_statistics['A'],
                                 'G':counter_seq_coding_nucleotide_statistics['G'],
                                 'C':counter_seq_coding_nucleotide_statistics['C']}


    for i in xrange (1,len(_features)-1):
        try:
            if ((_features[i].location.end.position == _features[i+1].location.end.position) and
                ((_features[i].location.start.position == _features[i+1].location.start.position))) :
                #we pass if there feature is duplicated
                pass
            else:
                if _features[i].type == 'D-loop':
                    d_loop_start = _features[i].location.start.position
                    d_loop_end = _features[i].location.end.position
                else:
                    #calcualte statistics for the analysis
                    feature_as_string  = _features[i].extract(mt_dna).tostring()
                    _total_protein_coding_nucleotide_statistics = Counter(feature_as_string)

                    total_protein_coding_seq_length += len(feature_as_string)
                    total_protein_coding_GC_count += feature_as_string.count('GC')

                    nucleotide_protein_coding_statistics = {'T':_total_protein_coding_nucleotide_statistics['T'],
                                     'A':_total_protein_coding_nucleotide_statistics['A'],
                                     'G':_total_protein_coding_nucleotide_statistics['G'],
                                     'C':_total_protein_coding_nucleotide_statistics['C']}

                    total_protein_coding_nucleotide_statistics += Counter (nucleotide_protein_coding_statistics)
        except:
            continue

    #check d_loop range boundaries
    if (d_loop_start == 0) and (d_loop_end == 0):
        try:
            d_loop_start, d_loop_end = get_largest_empty_region(_features)
        except:
            d_loop_start, d_loop_end = (0,1)

    d_loop_as_string = sequence_as_string[d_loop_start:d_loop_end]
    d_loop_GC_count = d_loop_as_string.count('GC')
    d_loop_seq_length = len(d_loop_as_string)

    counter_d_loop_nucleotide_statistics =  Counter(d_loop_as_string)
    total_d_loop__nucleotide_statistics = {'T':counter_d_loop_nucleotide_statistics['T'],
                                 'A':counter_d_loop_nucleotide_statistics['A'],
                                 'G':counter_d_loop_nucleotide_statistics['G'],
                                 'C':counter_d_loop_nucleotide_statistics['C']}


    return ({
                    'Total_seq_coding_GC_%': (total_seq_coding_nucleotide_statistics['G'] +
                                                total_seq_coding_nucleotide_statistics['C']*100.00) /
                                                 total_seq_coding_seq_length,
                    'Total_seq_coding_seq_length':total_seq_coding_seq_length,
                    'Total_seq_coding_G_count' : total_seq_coding_nucleotide_statistics['G'],
                    'Total_seq_coding_C_count' : total_seq_coding_nucleotide_statistics['C'],
                    'Total_seq_coding_A_count' : total_seq_coding_nucleotide_statistics['A'],
                    'Total_seq_coding_T_count' : total_seq_coding_nucleotide_statistics['T'],
                    'Total_protein-coding_seq_length': total_protein_coding_seq_length,
                    'Total_protein-coding_GC_%': (total_protein_coding_nucleotide_statistics['G'] +
                                                total_protein_coding_nucleotide_statistics['C'])*100.00 /
                                                 total_protein_coding_seq_length,
                    'Total_protein-coding_G_count' : total_protein_coding_nucleotide_statistics['G'],
                    'Total_protein-coding_C_count' : total_protein_coding_nucleotide_statistics['C'],
                    'Total_protein-coding_A_count' : total_protein_coding_nucleotide_statistics['A'],
                    'Total_protein-coding_T_count' : total_protein_coding_nucleotide_statistics['T'],
                    'D-loop_GC_%': (total_d_loop__nucleotide_statistics['G']+
                                            total_d_loop__nucleotide_statistics['C'])*100.00 / d_loop_seq_length ,
                    'D-loop_seq_length' : d_loop_seq_length,
                    'D-loop_G_count': total_d_loop__nucleotide_statistics['G'],
                    'D-loop_C_count': total_d_loop__nucleotide_statistics['C'],
                    'D-loop_A_count': total_d_loop__nucleotide_statistics['A'],
                    'D-loop_T_count': total_d_loop__nucleotide_statistics['T'],
                    })
