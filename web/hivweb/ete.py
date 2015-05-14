def orderByDistance(t):
    """
    Reorder the nodes of a tree in-place.
    Subtrees with leaves further from the root are placed after 
    siblings with smaller distance to the leaves.
    NOTE: this routine modifies the tree IN-PLACE!!
    """
    def distance(n):
        _, d = n.get_farthest_leaf()
        print(d)
        return d
    
    # preorder traversal of the tree, ordering children by number of descendants
    for node in t.traverse("preorder"):
        # sort the children in the desired order
        node.children = sorted(node.children, key=distance)


@app.route('/ete/<transmit>/<tsi_donor>/<tsi_acceptor>/<clockmodel>/mcc.svg')
def mcc_ete_tree_svg(transmit, tsi_donor, tsi_acceptor, clockmodel):
    # tree = process.tree_svg(os.path.join("../sims/runs", transmit, tsi_donor, tsi_acceptor, clockmodel, 'mcc.tree'), compress=False)
    # resp = make_response(tree)
    # resp.headers['Content-Type'] = 'image/svg+xml'
    # resp.headers['Content-Encoding'] = 'gzip'
    resp = None
    from ete2 import Tree
    import tempfile

    treefile = os.path.join("../sims/runs", transmit, tsi_donor, tsi_acceptor, clockmodel, 'mcc.tree')
    with open(treefile) as tf:
        tree = Tree(tf.read())
    orderByDistance(tree)
    
    svg, _ = tree.render("%%return", w=6.5, h=4, units="in", dpi=80)
    resp = make_response(svg)
    resp.headers['Content-Type'] = 'image/svg+xml'
    resp.headers['Content-Encoding'] = 'gzip'
        
    return resp
 
