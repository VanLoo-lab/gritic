
import pandas as pd
from gritic import gritictimer,sampletools



def load_copy_number_table(copy_number_path):
    copy_number_table = pd.read_csv(copy_number_path,sep='\t',dtype={'Chromosome':str,'Segment_Start':int,'Segment_End':int,'Major_CN':int,'Minor_CN':int})
    return copy_number_table

def load_mutation_table(mutation_table_path):
    mutation_table = pd.read_csv(mutation_table_path,sep='\t',dtype={'Chromosome':str,'Position':int,'Tumor_Ref_Count':int,'Tumor_Alt_Count':int})
    return mutation_table

def load_subclone_table(subclone_table_path):
    subclone_table = pd.read_csv(subclone_table_path,sep='\t',dtype={'Cluster':str,'Subclone_CCF':float,'Subclone_Fraction':float})
    return subclone_table

if __name__ == '__main__':
    
    sample_id = 'TEST_ID'
    output_dir = 'examples/output'
    sample_purity = 0.5
    sample_wgd_status = True

    plot_trees = True

    copy_number_table = load_copy_number_table('examples/cn_table_example.tsv')
    mutation_table = load_mutation_table('examples/snv_table_example.tsv')
    #if no subclone table set to None
    subclone_table = load_subclone_table('examples/subclone_table_example.tsv')
    
    sample = sampletools.Sample(mutation_table,copy_number_table,subclone_table,sample_id,sample_purity)
    gritictimer.process_sample(sample,output_dir,plot_trees=plot_trees,wgd_override=sample_wgd_status)