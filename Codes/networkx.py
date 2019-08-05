import networkx as nx
import sys
import argparse
parser = argparse.ArgumentParser(description='Calculate centrality metrics for PCST .sif network')
args = parser.parse_args()
parser.add_argument('-sif','--siffile', help='SIF, Network file. Can be output from OmicsIntegrator',required=True)
parser.add_argument('-o','--out', help='Output txt file',required=True)
g = nx.read_edgelist(args.siffile, create_using=None, delimiter='\tpp\t', nodetype=str, data=True, edgetype=None, encoding='utf-8')
output_file = str(args.out)
sys.stdout=open(args.out,"w")
degree_c = nx.degree_centrality(g)
pageR = nx.pagerank(g)
eigen_c = nx.eigenvector_centrality_numpy(g)
closeness_c = nx.closeness_centrality(g)
betweenness_c = nx.betweenness_centrality(g)
print "Node Name" + '\t' + "Degree Centrality" + '\t' + "Closeness Centrality" + '\t' + "Betweenness Centrality" + '\t' + "EigenVector Centrality" + '\t' + "Pagerank"
for node in degree_c:
    print str(node) + '\t' + str(degree_c[node]) + '\t' + str(closeness_c[node]) + '\t' + str(betweenness_c[node]) + '\t' + str(eigen_c[node]) + '\t' + str(pageR[node])
