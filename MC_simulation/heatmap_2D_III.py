import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns; sns.set()
import SeabornFig2Grid as sfg


iris = sns.load_dataset("iris")
g1 = sns.PairGrid(iris, hue="species")
g1.map(plt.scatter, s=5)
# A FacetGrid
#g2 = sns.FacetGrid(tips, col="time",  hue="smoker")
#g2.map(plt.scatter, "total_bill", "tip", edgecolor="w")
# A JointGrid
g3 = sns.jointplot("sepal_width", "petal_length", data=iris,
                   kind="kde", space=0, color="g")


fig = plt.figure(figsize=(13,8))
gs = gridspec.GridSpec(1, 2)
mg1 = sfg.SeabornFig2Grid(g1, fig, gs[0])
#mg2 = sfg.SeabornFig2Grid(g2, fig, gs[2])
mg3 = sfg.SeabornFig2Grid(g3, fig, gs[1])

gs.tight_layout(fig)
#gs.update(top=0.7)

plt.show()
