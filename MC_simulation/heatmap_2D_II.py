# library & dataset
import seaborn as sns
import matplotlib.pyplot as plt
df = sns.load_dataset('iris')
 
# Custom the inside plot: options are: “scatter” | “reg” | “resid” | “kde” | “hex”
scatter = sns.jointplot(x=df["sepal_length"], y=df["sepal_width"], kind='scatter')
hex     = sns.jointplot(x=df["sepal_length"], y=df["sepal_width"], kind='hex')
kde     = sns.jointplot(x=df["sepal_length"], y=df["sepal_width"], kind='kde')
plot_type = [scatter, hex, kde]
for i in plot_type:
    plt.savefig('heatmap'+i+'.png')
    plt.show()

sns.jointplot(x=df["sepal_length"], y=df["sepal_width"], kind='scatter', s=200, color='m', edgecolor="skyblue", linewidth=2)
 
# Custom the color
sns.set(style="white", color_codes=True)
sns.jointplot(x=df["sepal_length"], y=df["sepal_width"], kind='kde', color="skyblue")

plt.savefig('heatmap4.png')
plt.show()

