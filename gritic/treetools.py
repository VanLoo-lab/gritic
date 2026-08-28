import networkx as nx
import itertools
import hashlib
import random
import warnings

from gritic.tableschemas import NODE_PHASING_LABELS


ALLELE_ATTRIBUTE = 'Allele'
EMPTY_MINOR_COMPONENT_HASH = '<empty-minor>'
_ORDERED_TREE_HASH_VERSION = 'ordered-alleles-v1'
_ROUTE_IDENTITY_NODE_ATTRIBUTES = (
    ALLELE_ATTRIBUTE,
    'WGD_Symbol',
    'Terminal_Node',
)

#problem posed  https://leetcode.com/problems/all-possible-full-binary-trees/
#code from https://www.youtube.com/watch?v=nZtrZPTTCAo
# Definition for a binary tree node.

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
                plotting_tree.nodes[new_predecessor][ALLELE_ATTRIBUTE] = (
                    route_tree.nodes[node][ALLELE_ATTRIBUTE]
                )
                
                plotting_tree.nodes[new_loss]['WGD_Symbol']=False
                plotting_tree.nodes[new_loss]['Loss_Symbol']=True
                plotting_tree.nodes[new_loss][ALLELE_ATTRIBUTE] = (
                    route_tree.nodes[node][ALLELE_ATTRIBUTE]
                )
                if predecessor is not None:
                    
                    plotting_tree.remove_edge(predecessor,node)
                    plotting_tree.add_edge(predecessor,new_predecessor)

    return plotting_tree

def _get_allele_node_sets(tree):
    allele_nodes = {allele: set() for allele in NODE_PHASING_LABELS}
    for node, attributes in tree.nodes(data=True):
        allele = attributes.get(ALLELE_ATTRIBUTE)
        if allele not in allele_nodes:
            raise ValueError(
                f'Every route node must have {ALLELE_ATTRIBUTE}=Major or '
                f'Minor; node {node!r} has {allele!r}'
            )
        allele_nodes[allele].add(node)

    major_nodes = allele_nodes['Major']
    minor_nodes = allele_nodes['Minor']
    if not major_nodes:
        raise ValueError('A route tree must contain a Major allele component')

    for allele, nodes in allele_nodes.items():
        if not nodes:
            continue
        component = tree.subgraph(nodes).to_undirected()
        if not nx.is_connected(component):
            raise ValueError(
                f'The {allele} allele must form exactly one connected component'
            )

    for source, target in tree.edges:
        if tree.nodes[source][ALLELE_ATTRIBUTE] != tree.nodes[target][ALLELE_ATTRIBUTE]:
            raise ValueError('Route-tree edges cannot connect different alleles')

    return major_nodes, minor_nodes


def _get_component_hash(component):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message=(
                r'^The hashes produced for directed graphs changed in '
                r'version v3\.5 due to a bugfix to track in and out edges '
                r'separately \(see documentation\)\.$'
            ),
            category=UserWarning,
        )
        return nx.weisfeiler_lehman_graph_hash(
            component,
            node_attr='WGD_Symbol',
        )


def get_ordered_component_hashes(tree):
    """Return the Major and Minor topology hashes in explicit allele order."""
    major_nodes, minor_nodes = _get_allele_node_sets(tree)
    major_hash = _get_component_hash(tree.subgraph(major_nodes))
    minor_hash = EMPTY_MINOR_COMPONENT_HASH
    if minor_nodes:
        minor_hash = _get_component_hash(tree.subgraph(minor_nodes))
    return major_hash, minor_hash


def _combine_ordered_component_hashes(major_hash, minor_hash):
    identity = (
        f'{_ORDERED_TREE_HASH_VERSION}|Major:{major_hash}|Minor:{minor_hash}'
    )
    return hashlib.md5(identity.encode('utf-8')).hexdigest()


def get_tree_hash(tree):
    """Hash a route forest while preserving Major/Minor component identity."""
    return _combine_ordered_component_hashes(
        *get_ordered_component_hashes(tree)
    )


def _route_trees_are_isomorphic(first_tree, second_tree):
    """Compare every node attribute that contributes to route identity."""
    missing_attribute = object()

    def node_match(first_attributes, second_attributes):
        return all(
            first_attributes.get(attribute, missing_attribute)
            == second_attributes.get(attribute, missing_attribute)
            for attribute in _ROUTE_IDENTITY_NODE_ATTRIBUTES
        )

    return nx.is_isomorphic(
        first_tree,
        second_tree,
        node_match=node_match,
    )


def _register_ordered_tree(tree_store, tree_id, tree, context):
    """Store one isomorphism representative and reject true hash collisions."""
    retained_tree = tree_store.get(tree_id)
    if retained_tree is None:
        tree_store[tree_id] = tree
        return True

    if not _route_trees_are_isomorphic(retained_tree, tree):
        raise ValueError(
            f'Ordered route hash collision for {tree_id!r} while {context}: '
            'non-isomorphic route trees share one identifier'
        )
    return False


def _get_component_leaf_count(tree, component_nodes):
    return sum(tree.out_degree(node) == 0 for node in component_nodes)


def get_mirror_tree(tree):
    """Return the balanced route obtained by exchanging its allele identities."""
    major_nodes, minor_nodes = _get_allele_node_sets(tree)
    major_leaf_count = _get_component_leaf_count(tree, major_nodes)
    minor_leaf_count = _get_component_leaf_count(tree, minor_nodes)
    if (
        major_leaf_count <= 0
        or minor_leaf_count <= 0
        or major_leaf_count != minor_leaf_count
    ):
        raise ValueError(
            'Only balanced two-allele routes with equal positive terminal '
            'copy counts have a mirror route'
        )

    mirror_tree = tree.copy()
    for node in major_nodes:
        mirror_tree.nodes[node][ALLELE_ATTRIBUTE] = 'Minor'
    for node in minor_nodes:
        mirror_tree.nodes[node][ALLELE_ATTRIBUTE] = 'Major'
    return mirror_tree


def get_mirror_tree_hash(tree):
    """Return the ordered route ID of a balanced route's mirror."""
    return get_tree_hash(get_mirror_tree(tree))


def _get_new_node_id(tree):
    node_id = len(tree)
    while node_id in tree:
        node_id += 1
    return node_id


def convert_to_nx_tree(tree,allele,D=None,current_node_id=0):
    if allele not in NODE_PHASING_LABELS:
        raise ValueError(f'allele must be Major or Minor, not {allele!r}')
    if D is None:
        D = nx.DiGraph()
    if current_node_id not in D:
        D.add_node(current_node_id)

    current_attributes = D.nodes[current_node_id]
    current_allele = current_attributes.get(ALLELE_ATTRIBUTE)
    if current_allele is not None and current_allele != allele:
        raise ValueError(
            f'Cannot extend {current_allele!r} allele node {current_node_id!r} '
            f'with a {allele!r} allele tree'
        )
    current_attributes[ALLELE_ATTRIBUTE] = allele
    current_attributes.setdefault('WGD_Symbol', False)
    current_attributes['Terminal_Node'] = tree.left is None

    if tree.left is not None:
        new_node_id = _get_new_node_id(D)
        D.add_node(
            new_node_id,
            WGD_Symbol=False,
            Terminal_Node=False,
            Allele=allele,
        )
        D.add_edge(current_node_id,new_node_id)
        convert_to_nx_tree(tree.left,allele,D,new_node_id)
        
        new_node_id = _get_new_node_id(D)
        D.add_node(
            new_node_id,
            WGD_Symbol=False,
            Terminal_Node=False,
            Allele=allele,
        )
        D.add_edge(current_node_id,new_node_id)
        convert_to_nx_tree(tree.right,allele,D,new_node_id)
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
#this method inefficient and should be changed
def get_wgd_trees(tree,wgd_trees_status):
    internal_nodes  = [node for node in tree.nodes if not tree.nodes[node]['Terminal_Node'] ]

    wgd_trees = []
    wgd_tree_store = {}
    for node_length in range(len(internal_nodes)+1):
        for wgd_nodes in itertools.combinations(internal_nodes,node_length):
            # one node combination is always acceptable
            if not valid_wgd_node_combination(tree,wgd_nodes,wgd_trees_status):
                continue
            
            wgd_tree = convert_to_wgd_tree(tree,wgd_nodes)
            #inelegant way of removing the isomoprhic wgd trees
            wgd_hash = get_tree_hash(wgd_tree)

            if _register_ordered_tree(
                wgd_tree_store,
                wgd_hash,
                wgd_tree,
                f'deduplicating {wgd_trees_status} WGD routes',
            ):
                wgd_trees.append(wgd_tree)
    return wgd_trees


def check_tree(tree):
    for node in tree:
        if len(tree.out_edges(node)) != 0 and len(tree.out_edges(node)) != 2:
            return False
    return True


def convert_to_nx_trees(trees,allele):
    nx_trees = []
    for tree in trees:
        nx_tree = convert_to_nx_tree(tree,allele)
        assert check_tree(nx_tree)
        nx_trees.append(nx_tree)
    return nx_trees
def get_nx_trees(major_cn,minor_cn,wgd_status,wgd_trees_status):
    major_cn_trees = convert_to_nx_trees(allPossibleFBT(major_cn),'Major')
    minor_cn_trees = convert_to_nx_trees(allPossibleFBT(minor_cn),'Minor')

    if minor_cn ==0:
        combined_trees = major_cn_trees
    else:
        combined_trees = []
        for major_tree in major_cn_trees:
            for minor_tree in minor_cn_trees:
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
        _register_ordered_tree(
            tree_store,
            tree_id,
            tree,
            (
                f'building the {major_cn}+{minor_cn} route store '
                f'(wgd_status={wgd_status}, '
                f'wgd_trees_status={wgd_trees_status})'
            ),
        )
    
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
    _get_allele_node_sets(tree)
    return {
        node: attributes[ALLELE_ATTRIBUTE]
        for node, attributes in tree.nodes(data=True)
    }


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
    major_nodes, minor_nodes = _get_allele_node_sets(tree)
    major_tree = tree.subgraph(major_nodes)
    minor_tree = nx.empty_graph(0, create_using=nx.DiGraph())
    if minor_nodes:
        minor_tree = tree.subgraph(minor_nodes)
    return major_tree, minor_tree


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
    import matplotlib.pyplot as plt

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
    fig = plt.figure(figsize=(11, 6))
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
    major_component, minor_component = split_tree(route_tree)
    major_pos = hierarchy_pos(major_component)
    if minor_component.number_of_nodes():
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
