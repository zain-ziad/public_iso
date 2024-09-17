import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import io
from PIL import Image

def make_bar_plots(self, data, x, y, upper_lim, lower_lim, xlabel, ylabel):
    plots = []
    fullsize_images = []
    try:
        clipped_data = data[(data >= lower_lim) & (data <= upper_lim)]
        fig = plt.Figure(figsize=tuple(i/100 for i in (x, y)), dpi=200)
        ax = fig.add_subplot(111)
        with sns.axes_style("darkgrid"):
            sns.histplot(clipped_data, kde=True, color='#1C72AE', bins = 'auto', ax=ax)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.axis('tight')

        buf = io.BytesIO()                                                          # Initialize a memory buffer
        fig.savefig(buf, format='png')                                              # save the fig as a png in buffer
        buf.seek(0)                                                                 # rewind to beginning of file
        original_img = Image.open(buf)                                              # set that as original_img
        fullsize_images.append(original_img)                                        # save the original image in a list to retrive it for closeup

        new_dimensions = (int(original_img.width / 2), int(original_img.height / 2))
        resized_img = original_img.resize(new_dimensions, Image.Resampling.LANCZOS)
        plots.append(resized_img)

    except Exception as e:
        print(e)

    finally:
        args = (plots, fullsize_images)
        self.root.command_queue.put((f'filterdata_finalize_plotting', args))
        self.root.command_queue.put(('destroy_progress_window', {}))

def make_violin_plots(self, data, x, y, upper_lim, lower_lim, xlabel):
    plots = []
    fullsize_images = []
    try:
        clipped_data = data[(data >= lower_lim) & (data <= upper_lim)]
        fig = plt.Figure(figsize=tuple(i/100 for i in (x, y)), dpi=200)
        ax = fig.add_subplot(111)
        with sns.axes_style("darkgrid"):
            sns.violinplot(clipped_data, color='#1C72AE', alpha = 1, linewidth = 0.5, ax=ax)
            sns.stripplot(clipped_data, jitter=0.4, size=1, color='black', alpha = 0.5, ax=ax)
        ax.set(xlabel=xlabel, ylabel='')
        ax.axis('tight')

        buf = io.BytesIO()                                                          # Initialize a memory buffer
        fig.savefig(buf, format='png')                                              # save the fig as a png in buffer
        buf.seek(0)                                                                 # rewind to beginning of file
        original_img = Image.open(buf)                                              # set that as original_img
        fullsize_images.append(original_img)                                        # save the original image in a list to retrive it for closeup

        new_dimensions = (int(original_img.width / 2), int(original_img.height / 2))
        resized_img = original_img.resize(new_dimensions, Image.Resampling.LANCZOS)
        plots.append(resized_img)

    except Exception as e:
        print(e)

    finally:
        args = (plots, fullsize_images)
        self.root.command_queue.put((f'filterdata_finalize_plotting', args))
        self.root.command_queue.put(('destroy_progress_window', {}))
