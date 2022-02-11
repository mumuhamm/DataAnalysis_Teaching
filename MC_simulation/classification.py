import matplotlib.pyplot as plt
import itertools
import numpy as np
from sklearn import svm , datasets
from sklearn.metrics import classification_report
plt.style.use('seaborn-poster')

iris = datasets.load_iris()
print(iris.feature_names)
print(iris.data[:20])
print('the %d data sample has  %d features'%(iris.data.shape[0], iris.data.shape[1]))

print(iris.target_names)
print(set(iris.target))
X = iris.data[:, [0, 2]] # only two features
y = iris.target
target_names = iris.target_names
feature_names = iris.feature_names
# get the classes
n_class = len(set(y))
print('The %d classes in the data'%(n_class))

colors = ['b', 'g', 'r']
symbols = ['o', 'o', 'o']
plt.figure(figsize = (6,6))
for i, c, s in (zip(range(n_class), colors, symbols)):
    ix = y == i
    plt.scatter(X[:, 0][ix], X[:, 1][ix], color = c, marker = s, s = 60, label = target_names[i])
plt.legend(loc = 2, scatterpoints = 1)
plt.xlabel('Feature 1 - ' + feature_names[0])
plt.ylabel('Feature 2 - ' + feature_names[2])
plt.show()

clf = svm.SVC(kernel ='linear')
clf.fit(X,y)
print(clf.predict(X))

def des_boundary(X, y, clf, title = None):
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.01), np.arange(y_min, y_max, 0.01))
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    plt.figure(figsize = (10, 8))
    plt.contourf(xx, yy, Z, alpha=0.4)
    for i, c, s in (zip(range(n_class), colors, symbols)):
        ix = y == i
        plt.scatter(X[:, 0][ix], X[:, 1][ix], color = c, marker = s, s = 60, label = target_names[i])
    if title is not None:
        plt.title(title)
    plt.xlabel('Feature 1')
    plt.ylabel('Feature 2')
    plt.show()
des_boundary(X, y, clf)
