ROUTE_KEY_COLUMNS = ['Sample_ID', 'Segment_ID', 'Route']
SEGMENT_KEY_COLUMNS = ['Sample_ID', 'Segment_ID']
NODE_PHASING_LABELS = ('Major', 'Minor')
TIMING_REPRESENTATION_COLUMN = 'Timing_Representation'
UNIFORM_NO_GAIN_REPRESENTATION = 'Uniform_No_Gain'
ROUTE_PARTICLES_REPRESENTATION = 'Route_Particles'
TIMING_REPRESENTATIONS = (
    UNIFORM_NO_GAIN_REPRESENTATION,
    ROUTE_PARTICLES_REPRESENTATION,
)

SEGMENT_METADATA_COLUMNS = [
    'Sample_ID',
    'Segment_ID',
    'Chromosome',
    'Segment_Start',
    'Segment_End',
    'Major_CN',
    'Minor_CN',
    'Total_CN',
    'N_Mutations',
    'Mutation_Rate',
    'WGD_Status',
]

ROUTE_TABLE_COLUMNS = [
    'Sample_ID',
    'Segment_ID',
    'Route',
    TIMING_REPRESENTATION_COLUMN,
    'Chromosome',
    'Segment_Start',
    'Segment_End',
    'Major_CN',
    'Minor_CN',
    'Total_CN',
    'N_Mutations',
    'Mutation_Rate',
    'Probability',
    'Penalized_Probability',
    'Average_N_Events',
    'Average_Pre_WGD_Losses',
    'Average_Post_WGD_Losses',
    'Time',
    'Density',
    'WGD_Status',
    'WGD_Timing',
    'WGD_Timing_CI_Low',
    'WGD_Timing_CI_High',
]

GAIN_TIMING_TABLE_COLUMNS = [
    'Sample_ID',
    'Segment_ID',
    'Route',
    'Node',
    'Node_Phasing',
    'Timing',
    'Timing_CI_Low',
    'Timing_CI_High',
]

GAIN_DRAW_COLUMNS = SEGMENT_METADATA_COLUMNS + [
    'Posterior_Sample_Index',
    'Route',
    'Node',
    'Node_Phasing',
    'Gain_Timing',
    'WGD_Timing',
    'Gain_Index',
]

ROUTE_DRAW_COLUMNS = SEGMENT_METADATA_COLUMNS + [
    'Posterior_Sample_Index',
    'Route',
    'WGD_Timing',
]

SUMMARY_VALUE_COLUMNS = [
    'Gain_Index',
    'Proportion',
    'Timing_Median',
    'Timing_Low_CI',
    'Timing_High_CI',
    'Pre_WGD_Probability',
    'Post_WGD_Probability',
    'WGD_Timing_Median',
    'WGD_Timing_Low_CI',
    'WGD_Timing_High_CI',
]

SUMMARY_COLUMNS = SEGMENT_METADATA_COLUMNS + SUMMARY_VALUE_COLUMNS
