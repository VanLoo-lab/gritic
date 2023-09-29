#thanks to 
#https://www.geeksforgeeks.org/find-the-point-where-maximum-intervals-overlap/
#for the basic idea
def get_ids_with_maximum_overlap(segment_ci_store,segment_width_store):

    segment_ids_sorted_by_ci_low = sorted(segment_ci_store.keys(),key=lambda segment_id:segment_ci_store[segment_id][0])
    segment_ids_sorted_by_ci_high= sorted(segment_ci_store.keys(),key=lambda segment_id:segment_ci_store[segment_id][1])
    
    current_segments_in = [segment_ids_sorted_by_ci_low[0]]
    segments_with_max_overlap = list(current_segments_in)
    max_segment_overlap_width = sum([segment_width_store[segment_id] for segment_id in current_segments_in])

    n_segments = len(segment_ci_store.keys())
    i = 1
    j = 0
    best_overlap_timing =None
    while i<n_segments and j <n_segments:
        next_arrival_time = segment_ci_store[segment_ids_sorted_by_ci_low[i]][0]
        next_leaving_time = segment_ci_store[segment_ids_sorted_by_ci_high[j]][1]


        if next_arrival_time <= next_leaving_time:
            joining_segment_id = segment_ids_sorted_by_ci_low[i]

            current_segments_in.append(joining_segment_id)
            
            current_segment_overlap_width = sum([segment_width_store[segment_id] for segment_id in current_segments_in])
            if current_segment_overlap_width > max_segment_overlap_width:
                segments_with_max_overlap = list(current_segments_in)
                max_segment_overlap_width = current_segment_overlap_width
                best_overlap_timing = next_arrival_time
            i += 1
        else:
            leaving_segment_id = segment_ids_sorted_by_ci_high[j]
            current_segments_in.remove(leaving_segment_id)
            j+=1
    return segments_with_max_overlap,max_segment_overlap_width,best_overlap_timing