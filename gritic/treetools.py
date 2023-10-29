import networkx as nx
import itertools
import hashlib
import numpy as np
#problem posed  https://leetcode.com/problems/all-possible-full-binary-trees/
#code from https://www.youtube.com/watch?v=nZtrZPTTCAo
# Definition for a binary tree node.

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (11,6)
#I modded to make isomorphic
class TreeNode:
     def __init__(self, left=None, right=None):
        self.left = left
        self.right = right
def allPossibleFBT(n):

    #return the list of fbt with n nodes
    def backtrack(n):
        if n  ==0:
            return []
        if n  ==1:
            return [TreeNode()]
        
        res = []
        for left_leaves in range(1,n):
            right_leaves = n-left_leaves
            if right_leaves < left_leaves:
                continue
            leftTrees,rightTrees = backtrack(left_leaves),backtrack(right_leaves)
        
            for i,t1 in enumerate(leftTrees):
                for j, t2 in enumerate(rightTrees):

                    if left_leaves == right_leaves and j>i:
                        break
                    res.append(TreeNode(t1,t2))

        return res
    return backtrack(n)

def convert_to_plotting_tree(route_tree,wgd_timing,route_timing_summary,node_order,tol=1e-7):
    plotting_tree = route_tree.copy()
    nx.set_node_attributes(plotting_tree,False,'Loss_Symbol')
    for node in route_tree.nodes:
        wgd_in_path = False
        if route_tree.nodes[node]['WGD_Symbol']:
            wgd_in_path = True
        current_timing = route_timing_summary[node_order.index(node)]
        for ancestor in nx.ancestors(route_tree,node):
            if route_tree.nodes[ancestor]['WGD_Symbol']:
                wgd_in_path = True
        for descendant in nx.descendants(route_tree,node):
            if route_tree.nodes[descendant]['WGD_Symbol']:
                wgd_in_path = True

        if not wgd_in_path:
            predecessors = list(route_tree.predecessors(node))
            if len(predecessors)==0:
                predecessor = None
                prev_timing = 0
            else:
                predecessor = predecessors[0]
                prev_timing = route_timing_summary[node_order.index(predecessor)]
                
            if (current_timing > wgd_timing-tol) and prev_timing <= (wgd_timing+tol):
                new_predecessor = max(list(plotting_tree.nodes))+1
                new_loss =  max(list(plotting_tree.nodes))+2
                plotting_tree.add_edge(new_predecessor,node)
                plotting_tree.add_edge(new_predecessor,new_loss)

                plotting_tree.nodes[new_predecessor]['Loss_Symbol']=False
                plotting_tree.nodes[new_predecessor]['WGD_Symbol']=True
                
                plotting_tree.nodes[new_loss]['WGD_Symbol']=False
                plotting_tree.nodes[new_loss]['Loss_Symbol']=True
                if predecessor is not None:
                    
                    plotting_tree.remove_edge(predecessor,node)
                    plotting_tree.add_edge(predecessor,new_predecessor)

    return plotting_tree

def get_tree_hash(tree):
    connected_nodes = list(nx.connected_components(tree.to_undirected()))
    hashes = []
    for nodes in connected_nodes:
        sub_tree = tree.subgraph(nodes)
        tree_hash = nx.weisfeiler_lehman_graph_hash(sub_tree,node_attr='WGD_Symbol')
        hashes.append(tree_hash)
    joint_hash = '-'.join(sorted(hashes))
    final_hash = hashlib.md5(joint_hash.encode('utf-8')).hexdigest()
    return final_hash
    
def convert_to_nx_tree(tree,D=None,current_node_id=0):
    if D is None:
        D = nx.DiGraph()
        D.add_node(current_node_id,WGD_Symbol=False,Terminal_Node=False)
    if tree.left is not None:            
        new_node_id = len(D.nodes)
        D.add_node(new_node_id,WGD_Symbol=False,Terminal_Node=False)
        D.add_edge(current_node_id,new_node_id)
        convert_to_nx_tree(tree.left,D,new_node_id)
        
        new_node_id = len(D.nodes)
        D.add_node(new_node_id,WGD_Symbol=False,Terminal_Node=False)
        D.add_edge(current_node_id,new_node_id)
        convert_to_nx_tree(tree.right,D,new_node_id)
    else:
        D.nodes[current_node_id]['Terminal_Node'] = True
    
    return D

def convert_to_wgd_tree(tree,wgd_nodes):
    wgd_tree = tree.copy()
    for node in wgd_nodes:
        wgd_tree.nodes[node]['WGD_Symbol'] = True
    return wgd_tree


def valid_wgd_node_combination(tree,wgd_nodes,wgd_trees_status):
    wgd_nodes = set(wgd_nodes)
    
    if len(wgd_nodes)==0:
        if wgd_trees_status == 'Only_WGD' or wgd_trees_status == 'Only_WGD_2+2':
            return False
        return True
    
    if len(wgd_nodes)==1:
        if wgd_trees_status == 'No_WGD' or wgd_trees_status == 'Only_WGD_2+2':
            return False
        return True
    
    for node in wgd_nodes:
        descendants = nx.descendants(tree,node)
        if len(descendants.intersection(wgd_nodes))>0:
            return False
    return True
#this method is gross and highly inefficient and should be changed
def get_wgd_trees(tree,wgd_trees_status):
    internal_nodes  = [node for node in tree.nodes if not tree.nodes[node]['Terminal_Node'] ]

    wgd_trees = []
    wgd_tree_hashes = set()
    for node_length in range(len(internal_nodes)+1):
        for wgd_nodes in itertools.combinations(internal_nodes,node_length):
            # one node combination is always acceptable
            if not valid_wgd_node_combination(tree,wgd_nodes,wgd_trees_status):
                continue
            
            wgd_tree = convert_to_wgd_tree(tree,wgd_nodes)
            #crass way of removing the isomoprhic wgd trees
            wgd_hash = get_tree_hash(wgd_tree)

            if not wgd_hash in wgd_tree_hashes:

                wgd_trees.append(wgd_tree)
                wgd_tree_hashes.add(wgd_hash)
    return wgd_trees


def check_tree(tree):
    for node in tree:
        if len(tree.out_edges(node)) != 0 and len(tree.out_edges(node)) != 2:
            return False
    return True


def convert_to_nx_trees(trees):
    nx_trees = []
    for tree in trees:
        nx_tree = convert_to_nx_tree(tree)
        assert check_tree(nx_tree)
        nx_trees.append(nx_tree)
    return nx_trees
def get_nx_trees(major_cn,minor_cn,wgd_status,wgd_trees_status):
    major_cn_trees = convert_to_nx_trees(allPossibleFBT(major_cn))
    minor_cn_trees = convert_to_nx_trees(allPossibleFBT(minor_cn))

    if minor_cn ==0:
        combined_trees = major_cn_trees
    else:
        combined_trees = []
        for i,major_tree in enumerate(major_cn_trees):
            for j,minor_tree in enumerate(minor_cn_trees):
                if major_cn == minor_cn and j>i:
                    break
                combined_trees.append(nx.disjoint_union(major_tree,minor_tree))

    
    if wgd_status and not wgd_trees_status =='No_WGD':
        all_wgd_trees = []
        for nx_tree in combined_trees:
            wgd_trees = get_wgd_trees(nx_tree,wgd_trees_status)
            all_wgd_trees.extend(wgd_trees)
        full_trees = all_wgd_trees
    else:
        full_trees = combined_trees

    tree_store = {}
    for tree in full_trees:
        tree_id = get_tree_hash(tree)
        #tree_hash = nx.weisfeiler_lehman_graph_hash(tree,node_attr='WGD_Symbol')
        #tree_id = hashlib.md5(f'{major_cn}_{minor_cn}_{tree_hash}'.encode('utf-8')).hexdigest()
        tree_store[tree_id] = tree
    
    return tree_store


def get_possible_paths_recursive(tree,current_node,path):
    if  len(tree.out_edges(current_node))==0:
        return [path]

    possible_paths = []
    for edge in tree.out_edges(current_node):

        new_path = list(path)
        new_path.append(edge[1])
        possible_paths.extend(get_possible_paths_recursive(tree,edge[1],new_path))
    
    return possible_paths
#https://stackoverflow.com/questions/4122390/getting-the-root-head-of-a-digraph-in-networkx-python
def get_possible_paths(tree):
    root_nodes = [node for node,in_degree in tree.in_degree() if in_degree==0] 
    assert len(root_nodes)<=2
    all_possible_paths = []
    for root_node in root_nodes:
        
        all_possible_paths.extend(get_possible_paths_recursive(tree,root_node,[root_node]))
    return all_possible_paths


def get_wgd_paths_recursive(tree,current_node,path):
    if tree.nodes[current_node]['WGD_Symbol']:
        return [path]
    possible_paths = []
    for edge in tree.out_edges(current_node):

        new_path = list(path)
        new_path.append(edge[1])
        possible_paths.extend(get_wgd_paths_recursive(tree,edge[1],new_path))
    return possible_paths
#https://stackoverflow.com/questions/4122390/getting-the-root-head-of-a-digraph-in-networkx-python
def get_wgd_paths(tree):
    root_nodes = [node for node,in_degree in tree.in_degree() if in_degree==0] 
    assert len(root_nodes)<=2
    all_possible_paths = []
    for root_node in root_nodes:
        all_possible_paths.extend(get_wgd_paths_recursive(tree,root_node,[root_node]))
    return all_possible_paths

def get_node_phasing_tree(tree):
    node_phasing = {}
    connected_components = sorted(list(nx.connected_components(tree.to_undirected())),key=lambda x:len(x))[::-1]
    assert len(connected_components) <=2
    if len(connected_components) == 1:
        phasing = (np.nan,)
    elif len(connected_components[0]) == len(connected_components[1]):
        phasing = ('A','B')
    else:
        phasing = ('Major','Minor')
    for i,connected_component in enumerate(connected_components):
        for node in connected_component:
            node_phasing[node] = phasing[i]
    return node_phasing
def get_node_attributes(tree,wgd_status):
    node_attributes = {}
    node_phasing = get_node_phasing_tree(tree)
    for node in tree.nodes:
        node_attribute = {}
        
        predecessors = list(tree.predecessors(node))
        if len(predecessors) ==0:
            predecessor = None
        else:
            predecessor = predecessors[0]
        node_attribute['Predecessor'] = predecessor
        node_attribute['Successors'] = list(tree.successors(node))
        node_attribute['Ancestors'] = list(nx.ancestors(tree,node))
        node_attribute['WGD'] = tree.nodes[node]['WGD_Symbol']
        node_attribute['Phasing'] = node_phasing[node]
        descendants = set(nx.descendants(tree,node))
        n_final_ancestors = len([node for node in descendants if len(tree.adj[node]) ==0])
        final_mult = n_final_ancestors if n_final_ancestors > 0 else 1
        node_attribute['Multiplicity'] = final_mult

        if not wgd_status:
            node_attribute['WGD_Ordering'] = 'NA'
        elif tree.nodes[node]['WGD_Symbol']:
            node_attribute['WGD_Ordering'] = 'WGD'
        else:
            for descendant in nx.descendants(tree,node):
                if tree.nodes[descendant]['WGD_Symbol']:
                    node_attribute['WGD_Ordering'] = 'Pre'
                    break
            #gross
            if not 'WGD_Ordering' in node_attribute.keys():
                for ancestor in nx.ancestors(tree,node):
                    if tree.nodes[ancestor]['WGD_Symbol']:
                        node_attribute['WGD_Ordering'] = 'Post'
                        break
            if not 'WGD_Ordering' in node_attribute.keys():
                node_attribute['WGD_Ordering'] = 'Calculate'
        node_attributes[node] = node_attribute

    
    return node_attributes
def split_tree(tree):
        connected_components = sorted(list(nx.connected_components(tree.to_undirected())),key=lambda x:len(x))[::-1]
        major_tree = tree.subgraph(connected_components[0])
        minor_tree = nx.empty_graph(0,create_using=nx.DiGraph())
        if len(connected_components) >1:
            minor_tree = tree.subgraph(connected_components[1])
        return major_tree,minor_tree
def write_tree(route_tree,output_path=None):
    nodes_to_delete = []
    for node in route_tree.nodes():
        del route_tree.nodes[node]['Full_Timing']
        del route_tree.nodes[node]['Terminal_Node']
        if len(list(route_tree.successors(node)))==0:
            nodes_to_delete.append(node)
    
    route_tree.remove_nodes_from(nodes_to_delete)
    
    nx.write_graphml(route_tree,output_path)
def plot_tree(route_tree,title,output_path=None):
    node_colors = []
    for node in route_tree.nodes:
        if 'Loss_Symbol' in route_tree.nodes[node] and route_tree.nodes[node]['Loss_Symbol']:
            node_colors.append('#4d4d4d')
        elif len(list(route_tree.successors(node)))==0:
            node_colors.append('#F7483B')
        elif route_tree.nodes[node]['WGD_Symbol']:
            node_colors.append('#E8BF5E')
        else:
            node_colors.append('#509BCE')
    
    pos = get_combined_hierarchy_pos(route_tree)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.margins(0.1)
    plt.title(title)
    labels = nx.get_node_attributes(route_tree,'Label')
    
 
    
    
    nx.draw(route_tree,pos=pos,with_labels=True,labels=labels,node_color=node_colors,font_weight='bold',ax=ax1)

    if output_path is not None:
        plt.savefig(output_path)
    else:
        plt.show()
    plt.close(fig)
def get_combined_hierarchy_pos(route_tree):
    connected_components = sorted(list(nx.connected_components(route_tree.to_undirected())),key=lambda x:len(x))[::-1]
    major_component = route_tree.subgraph(connected_components[0])
    major_pos = hierarchy_pos(major_component)
    if len(connected_components) >1:
        minor_component = route_tree.subgraph(connected_components[1])
        minor_pos = hierarchy_pos(minor_component,x_offset=1)
    else:
        minor_pos = {}
    
    return {**major_pos,**minor_pos}
#https://stackoverflow.com/questions/29586520/can-one-get-hierarchical-graphs-from-networkx-with-python-3
def hierarchy_pos(G, root=None, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5,x_offset=0):

    '''
    From Joel's answer at https://stackoverflow.com/a/29597209/2966723.  
    Licensed under Creative Commons Attribution-Share Alike 
    
    If the graph is a tree this will return the positions to plot this in a 
    hierarchical layout.
    
    G: the graph (must be a tree)
    
    root: the root node of current branch 
    - if the tree is directed and this is not given, 
      the root will be found and used
    - if the tree is directed and this is given, then 
      the positions will be just for the descendants of this node.
    - if the tree is undirected and not given, 
      then a random choice will be used.
    
    width: horizontal space allocated for this branch - avoids overlap with other branches
    
    vert_gap: gap between levels of hierarchy
    
    vert_loc: vertical location of root
    
    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5, pos = None, parent = None,x_offset=0):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''

        if pos is None:
            pos = {root:(xcenter+x_offset,vert_loc)}
        else:
            pos[root] = (xcenter+x_offset, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)  
        if len(children)!=0:
            dx = width/len(children) 
            nextx = xcenter - width/2 - dx/2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G,child, width = dx, vert_gap = vert_gap, 
                                    vert_loc = vert_loc-vert_gap, xcenter=nextx,
                                    pos=pos, parent = root,x_offset=x_offset)
        return pos

            
    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter,x_offset=x_offset)