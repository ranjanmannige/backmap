# MCA
import os,copy,random


feature_dict = {i:label for i,label in zip(
				range(4),
				  ('sepal length in cm',
				  'sepal width in cm',
				  'petal length in cm',
				  'petal width in cm', ))}

import pandas as pd

df = pd.io.parsers.read_csv(
	filepath_or_buffer='iris.csv',
	#header=None,
	sep=',',
	)

if 0:
	df = pd.io.parsers.read_csv(
		filepath_or_buffer='https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data',
		header=None,
		sep=',',
		)
	
	df.columns = [l for i,l in sorted(feature_dict.items())] + ['class label']
	
	df.dropna(how="all", inplace=True) # to drop the empty line at file-end
	df.tail()
if 1:
	f = open('rcode_to_scop_class_rstep0.1.csv','r')
	lines = f.read().rstrip().split("\n")
	f.close()
	
	
	print "READING"
	newblock = lines[0]+"\n"
	for l in lines[1:]:
		
		linesplit = l.split(",")
		class_name = ".".join(linesplit[-1].split(".")[:1])
		if class_name in "abc":
			newblock += ",".join(linesplit[:-1])+","+class_name+"\n"
	print "DONE"
	
	
	feature_dict = {}
	header_line = lines[0].lstrip().rstrip().lstrip(",").split(",")
	for i in range(len(header_line[:-1])):
		feature_dict[i] = float(header_line[i])
	
	f = open("temp.csv","w")
	f.write(newblock.rstrip("\n"))
	f.close()
	
	df = pd.io.parsers.read_csv(
		filepath_or_buffer="temp.csv",
		#header=None,
		sep=',',
		)

from sklearn.preprocessing import LabelEncoder

start_column = 0
if list(df.columns.values)[0] == 'Unnamed: 0':
	start_column = 1

#X = df[[0,1,2,3]].values
X = df[range(start_column,df.shape[1]-1)].values # df.shape == (150, 5); range(df.shape[1]-1) = [0,1,2,3]
y = df['class label'].values

oldy = copy.deepcopy(y)
enc = LabelEncoder()
label_encoder = enc.fit(y)
y = label_encoder.transform(y) + 1

#label_dict = {1: 'Setosa', 2: 'Versicolor', 3:'Virginica'}
label_dict = {}
for a,b, in zip(y,oldy):
	if not a in label_dict:
		label_dict[a] = b

from matplotlib import pyplot as plt
import numpy as np
import math


column_names = list(df.columns.values)[start_column:-1]

nrows = int(np.sqrt(len(column_names)))
ncols = int(np.ceil(float(len(column_names))/nrows))

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12,6))

species_keys   = sorted(label_dict.keys())
species_colors = []
species_scatter_markers = []

marker_count = 0
possible_markers = ['o','+','x','|','_']   
#('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
for i in species_keys:
	species_colors.append((random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)))
	species_scatter_markers.append(possible_markers[marker_count])
	marker_count+=1

for ax,cnt in zip(axes.ravel(), range(len(column_names))):
	
	# set bin sizes
	min_b = math.floor(np.min(X[:,cnt]))
	max_b = math.ceil(np.max(X[:,cnt]))
	bins = np.linspace(min_b, max_b, 25)

	# plottling the histograms
	
	
	#for lab,col in zip(range(1,4), ('blue', 'red', 'green')):
	for lab,col in zip(species_keys,species_colors):
		ax.hist(X[y==lab, cnt],
				   color=col,
				   label='class %s' %label_dict[lab],
				   bins=bins,
				   alpha=0.5,)
	ylims = ax.get_ylim()
	
	# plot annotation
	leg = ax.legend(loc='upper right', fancybox=True, fontsize=8)
	leg.get_frame().set_alpha(0.5)
	ax.set_ylim([0, max(ylims)+2])
	ax.set_xlabel(feature_dict[cnt])
	ax.set_title('#%s' %str(cnt+1))

	# hide axis ticks
	ax.tick_params(axis="both", which="both", bottom="off", top="off",  
			labelbottom="on", left="off", right="off", labelleft="on")

	# remove axis spines
	ax.spines["top"].set_visible(False)  
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.spines["left"].set_visible(False)

axes[0][0].set_ylabel('count')
axes[1][0].set_ylabel('count')

fig.tight_layout()	   
plt.show()

np.set_printoptions(precision=4)

mean_vectors = []
#for cl in range(1,4):
for cl in species_keys:
	mean_vectors.append(np.mean(X[y==cl], axis=0))
	print('Mean Vector class %s: %s\n' %(cl, mean_vectors[cl-1]))


#S_W = np.zeros((4,4))
S_W = np.zeros((len(column_names),len(column_names)))
#for cl,mv in zip(range(1,4), mean_vectors):
for cl,mv in zip(range(1,len(column_names)), mean_vectors):
	class_sc_mat = np.zeros((len(column_names),len(column_names)))				  # scatter matrix for every class
	for row in X[y == cl]:
		row, mv = row.reshape(len(column_names),1), mv.reshape(len(column_names),1) # make column vectors
		class_sc_mat += (row-mv).dot((row-mv).T)
	S_W += class_sc_mat							 # sum class scatter matrices
print('within-class Scatter Matrix:\n', S_W)

overall_mean = np.mean(X, axis=0)

S_B = np.zeros((len(column_names),len(column_names)))
for i,mean_vec in enumerate(mean_vectors):  
	n = X[y==i+1,:].shape[0]
	mean_vec = mean_vec.reshape(len(column_names),1) # make column vector
	overall_mean = overall_mean.reshape(len(column_names),1) # make column vector
	S_B += n * (mean_vec - overall_mean).dot((mean_vec - overall_mean).T)

print('between-class Scatter Matrix:\n', S_B)

eig_vals, eig_vecs = np.linalg.eig(np.linalg.inv(S_W).dot(S_B))

for i in range(len(eig_vals)):
	eigvec_sc = eig_vecs[:,i].reshape(len(column_names),1)   
	print('\nEigenvector {}: \n{}'.format(i+1, eigvec_sc.real))
	print('Eigenvalue {:}: {:.2e}'.format(i+1, eig_vals[i].real))
	
for i in range(len(eig_vals)):
	eigv = eig_vecs[:,i].reshape(len(column_names),1)
	np.testing.assert_array_almost_equal(np.linalg.inv(S_W).dot(S_B).dot(eigv),eig_vals[i] * eigv,decimal=6, err_msg='', verbose=True)
print('ok')

# Make a list of (eigenvalue, eigenvector) tuples
eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]

# Sort the (eigenvalue, eigenvector) tuples from high to low
eig_pairs = sorted(eig_pairs, key=lambda k: k[0], reverse=True)

# Visually confirm that the list is correctly sorted by decreasing eigenvalues

print('Eigenvalues in decreasing order:\n')
for i in eig_pairs:
	print(i[0])
	
	
print('Variance explained:\n')
eigv_sum = sum(eig_vals)
for i,j in enumerate(eig_pairs):
	print('eigenvalue {0:}: {1:.2%}'.format(i+1, (j[0]/eigv_sum).real))

W = np.hstack((eig_pairs[0][1].reshape(len(column_names),1), eig_pairs[1][1].reshape(len(column_names),1)))
print('Matrix W:\n', W.real)


X_lda = X.dot(W)
assert X_lda.shape == (df.shape[0],2), "The matrix is not 2x150 dimensional."


from matplotlib import pyplot as plt

def plot_step_lda():

	ax = plt.subplot(111)
	
	for label,marker,color in zip(species_keys,species_scatter_markers,species_colors):
		plt.scatter(x=X_lda[:,0].real[y == label],
				y=X_lda[:,1].real[y == label],
				marker=marker,
				color=color,
				alpha=1,
				label=label_dict[label],
				)

	plt.xlabel('LD1')
	plt.ylabel('LD2')

	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.5)
	plt.title('LDA: Iris projection onto the first 2 linear discriminants')

	# hide axis ticks
	plt.tick_params(axis="both", which="both", bottom="off", top="off",  
			labelbottom="on", left="off", right="off", labelleft="on")

	# remove axis spines
	ax.spines["top"].set_visible(False)  
	ax.spines["right"].set_visible(False)
	ax.spines["bottom"].set_visible(False)
	ax.spines["left"].set_visible(False)	

	plt.grid()
	plt.tight_layout
	plt.show()

plot_step_lda()