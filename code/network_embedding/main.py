
import argparse
import numpy as np
import networkx as nx
import gene2vec
import warnings
warnings.filterwarnings(action='ignore', category=UserWarning, module='gensim')
from gensim.models import Word2Vec

def parse_args():
	'''
	Parses the gene2vec arguments.
	'''
	parser = argparse.ArgumentParser(description="Run gene2vec.")

	parser.add_argument('--input', nargs='?', default='data/breast_input_graph.edgelist',
	                    help='Input graph path')

	parser.add_argument('--input_na', nargs='?', default='data/breast_input_na.edgelist',
						help='Input graph_na path')

	parser.add_argument('--output', nargs='?', default='data/breast_out.txt',
	                    help='Embeddings path')

	parser.add_argument('--dimensions', type=int, default=128,
	                    help='Number of dimensions. Default is 128.')

	parser.add_argument('--walk-length', type=int, default=80,
	                    help='Length of walk per source. Default is 80.')

	parser.add_argument('--num-walks', type=int, default=10,
	                    help='Number of walks per source. Default is 10.')

	parser.add_argument('--window-size', type=int, default=5,
                    	help='Context size for optimization. Default is 10.')

	parser.add_argument('--iter', default=5, type=int,
                      help='Number of epochs in SGD')

	parser.add_argument('--workers', type=int, default=8,
	                    help='Number of parallel workers. Default is 8.')

	parser.add_argument('--p', type=float, default=1,
	                    help='Return hyperparameter. Default is 1.')

	parser.add_argument('--q', type=float, default=3,
	                    help='Inout hyperparameter. Default is 1.')

	parser.add_argument('--weighted', dest='weighted', action='store_true',
	                    help='Boolean specifying (un)weighted. Default is unweighted.')
	parser.add_argument('--unweighted', dest='unweighted', action='store_false')
	parser.set_defaults(weighted=False)

	parser.add_argument('--directed', dest='directed', action='store_true',
	                    help='Graph is (un)directed. Default is undirected.')
	parser.add_argument('--undirected', dest='undirected', action='store_false')
	parser.set_defaults(directed=False)

	return parser.parse_args()

def read_graph(inputdata):
	'''
	Reads the input network in networkx.
	'''
	if args.weighted:
		G = nx.read_edgelist(inputdata, nodetype=str, data=(('weight',float),), create_using=nx.DiGraph())
	else:
		G = nx.read_edgelist(inputdata, nodetype=str, create_using=nx.DiGraph())
		for edge in G.edges():
			G[edge[0]][edge[1]]['weight'] = 1

	if not args.directed:
		G = G.to_undirected()

	return G



def _filter_walks(walks, node_num):
	filter_walks = []
	for walk in walks:
		if int(walk[0]) <= node_num:
			fwalks = [nid for nid in walk if int(nid) <= node_num]
			filter_walks.append(fwalks)
	lf1 = len(filter_walks)
	return filter_walks

def learn_embeddings(walks):
	'''
	Learn embeddings by optimizing the Skipgram objective using SGD.
	'''
	walks = [map(str, walk) for walk in walks]

	model = Word2Vec(walks, size=args.dimensions, window=args.window_size, min_count=0, sg=1, workers=args.workers, iter=args.iter)
	model.wv.save_word2vec_format(args.output)
	
	return

def main(args):


    nx_G = read_graph(args.input)
    G = gene2vec.Graph(nx_G, args.directed, args.p, args.q)
    G.preprocess_transition_probs()
    walks_graph = G.simulate_walks(args.num_walks, args.walk_length)
    lnodes=len(G.alias_nodes)
    nx_G = read_graph(args.input_na)
    G = gene2vec.Graph(nx_G, args.directed, args.p ,args.q)
    G.preprocess_transition_probs()
    walks_na = G.simulate_walks(args.num_walks, args.walk_length*2)
    filter_na = G.filter_walk(walks_na, lnodes)
    walks = walks_graph+filter_na
    learn_embeddings(walks)



if __name__ == "__main__":
    args = parse_args()
    main(args)
