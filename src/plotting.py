import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


class Radar(object):

    def __init__(self, fig, titles):
        rect = [0.0, 0.0, 1.0, 1.0]
        self.n = len(titles)
        self.angles = map(lambda x: x%360, np.arange(90, 90+360, 360.0/self.n))
        self.axes = [fig.add_axes(rect, projection="polar", label="axes%d" % i) for i in range(self.n)]

        self.ax = self.axes[0]
        res, labels = self.ax.set_thetagrids(self.angles, labels=titles, fontsize=12, frac=1.2)

        labels[1].set_rotation(75)
        labels[4].set_rotation(285)
        labels[2].set_rotation(330)
        labels[3].set_rotation(30)

        for ax in self.axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)

        for ax, angle in zip(self.axes, self.angles):
            ax.set_rgrids([0.25, 0.5, 0.75, 1.0], angle=angle)
            ax.spines["polar"].set_visible(False)
            ax.set_ylim(0, 1.0)

    def plot(self, values, marker, fill, *args, **kw):
        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        self.ax.plot(angle, values, marker, alpha=0.9, *args, **kw)

        if fill:
            self.ax.fill(angle, values, alpha=0.3, *args, **kw)

    def set_title(self, title):
        self.ax.set_title(title, y=1.15, fontsize=14)


def spider_plot_compare(data, labels, title=None, styles=None):
    fig = plt.figure(figsize=(4, 4))
    fields = ['F1-score', 'Precision', 'Sensitivity', 'Specificity', 'Accuracy']
    radar = Radar(fig, fields)

    if not styles:
        colors = sns.color_palette('deep', n_colors=len(data))
        markers = '-'*len(data)
        fill = [True]*len(data)
        styles = zip(colors, markers, fill)

    for (color, marker, fill), entry in zip(styles, data):
        points = [entry[field] for field in fields]
        radar.plot(points, marker, fill, lw=2, color=color)

    if title:
        radar.set_title(title)

    radar.ax.legend(labels, loc=(0.9, 0.85))

