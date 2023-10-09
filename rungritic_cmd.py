import sys
import argparse
import pandas as pd
from gritic import dataloader,gritictimer,sampletools,gritictimer

#https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def load_copy_number_table(copy_number_path):
    copy_number_table = pd.read_csv(copy_number_path,sep='\t',dtype={'Chromosome':str,'Segment_Start':int,'Segment_End':int,'Major_CN':int,'Minor_CN':int})
    copy_number_table['Total_CN'] = copy_number_table['Major_CN']+copy_number_table['Minor_CN']
    return copy_number_table

def load_mutation_table(mutation_table_path):
    mutation_table = pd.read_csv(mutation_table_path,sep='\t',dtype={'Chromosome':str,'Position':int,'Tumor_Ref_Count':int,'Tumor_Alt_Count':int})
    return mutation_table

def load_subclone_table(subclone_table_path):
    subclone_table = pd.read_csv(subclone_table_path,sep='\t',dtype={'Cluster':str,'Subclone_CCF':float,'Subclone_Fraction':float})
    return subclone_table

if __name__ == '__main__':
    my_parser = argparse.ArgumentParser(allow_abbrev=False)
    
    my_parser.add_argument('--mutation_table', action='store', type=str, required=True, help='A tab separated file containing the following columns: Chromosome,Position,Tumor_Ref_Count,Tumor_Alt_Count')
    
    my_parser.add_argument('--copy_number_table', action='store', type=str, required=True, help='A tab separated file containing the following columns: Chromosome,Segment_Start,Segment_End,Major_CN,Minor_CN. Copy number should be integer values.')
    my_parser.add_argument('--purity', action='store', type=float, required=True, help='The purity of the sample as determined by the copy number caller.')

    my_parser.add_argument('--sample_id', action='store', type=str, required=True)
    my_parser.add_argument('--output', action='store', type=str, required=True, help='The output directory')


    my_parser.add_argument('--wgd_status', action='store', type=str2bool, required=False,default=None,help='If the sample is known to have undergone whole genome duplication, set this to True. If the sample is known to not have undergone whole genome duplication, set this to False. If the sample is of unknown WGD status, leave this blank.')
    my_parser.add_argument('--non_parsimony_penalty', action='store', type=str2bool, required=False,default=False,help='Apply additional penalty to non-parsimonious copy number route histories. Default is False.')
    my_parser.add_argument('--plot_trees', action='store', type=str, required=False,default=True,help='Plot the copy number trees. Default is True.')
    my_parser.add_argument('--subclone_table', action='store', type=str, required=False,default=None,help='A tab separated file containing the following columns: Cluster,Subclone_CCF,Subclone_Fraction (optional).' )
    
    args = my_parser.parse_args()
    sample_id = args.sample_id
    output_dir = args.output
    sample_purity = args.purity
    sample_wgd_status = args.wgd_status
    non_parsimony_penalty = args.non_parsimony_penalty
    plot_trees = args.plot_trees
    
    
    copy_number_table = load_copy_number_table(args.copy_number_table)
    copy_number_table = dataloader.merge_segments(copy_number_table)
    mutation_table = load_mutation_table(args.mutation_table)
    
    if args.subclone_table is None:
        subclone_table = None
    else:
        subclone_table = load_subclone_table(args.subclone_table)
        
 
    sample = sampletools.Sample(mutation_table,copy_number_table,subclone_table,sample_id,sample_purity)
    gritictimer.process_sample(sample,output_dir,plot_trees=plot_trees,wgd_override=sample_wgd_status,non_parsimony_penalty=non_parsimony_penalty)