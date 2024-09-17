from squidpy.pl import spatial_scatter
import io
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as patches
from PIL import Image
from holoviews.element.chart import Chart
from holoviews.plotting.bokeh import PointPlot
from holoviews import opts
import holoviews as hv
import param
import numpy as np
import datashader as ds
import datashader.transfer_functions as tf
from matplotlib import colormaps
import customtkinter as ctk
from core.assets.themecolors import COLORS

import scipy.sparse as sp
matplotlib.use('Agg')

class Circle(Chart):
    group = param.String(default='Circle', constant=True)
    size = param.Integer()

class CirclePlot(PointPlot):
    _plot_methods = dict(single='circle', batched='circle')
    style_opts = ['radius' if so == 'size' else so for so in PointPlot.style_opts if so != 'marker']

# A function to make plotting containers
def make_container(self, title, row, col, index, spatial_container = False, frame_grid_ypad = (0,0)):
    if not hasattr(self, 'container_parent'):
        raise AttributeError(f"container_parent not defined for class: {self}.")
    
    if not hasattr(self, 'listofcanvases'):
        raise AttributeError(f"listofcanvases not defined for class: {self}.")
    
    if not hasattr(self, 'listofCB_canvases') and spatial_container:
        raise AttributeError(f"listofCB_canvases not defined for class: {self}.")
    
    canvas_pack_ypad = (0,0) if spatial_container else (0,12)
        
    frame = ctk.CTkFrame(self.container_parent, fg_color=COLORS['white'], border_color=COLORS['border'], corner_radius=20, border_width=1)
    frame.grid(column = col, row = row, sticky = 'nsew', padx = (8,0), pady = frame_grid_ypad)
    label = ctk.CTkLabel(frame, text_color=COLORS['neutral 1'], font=('Poppins', 18), text=title)
    label.pack(side = 'top', pady = 10, fill = 'x', padx = 8)
    button = ctk.CTkButton(frame, text='', image=(ctk.CTkImage(light_image = Image.open('app\\core\\assets\\viewdata\\open_in_new.png'))),
                            height=24, width=24, corner_radius=8, fg_color=COLORS['white'], hover_color=COLORS['background'], command= lambda t = title, i = index: self.plot_closeup(title = t, index=i))
    button.place(relx=1.0, x=-16-50, y = 10)
    canvas = ctk.CTkCanvas(frame, background = COLORS['white'], highlightthickness = 0)
    canvas.pack(side = 'top', fill = 'both', expand = True, padx=2, pady=canvas_pack_ypad)
    self.listofcanvases.append(canvas)
    
    if spatial_container:
        self.cb_canvas_height = int(48 * self.root.scaling_factor)
        cb_canvas = ctk.CTkCanvas(frame, background = COLORS['white'], highlightthickness = 0, height = self.cb_canvas_height)
        cb_canvas.pack(fill = 'x', padx=2, pady=(0, 20))
        self.listofCB_canvases.append(cb_canvas)

def make_squidpy_plots(self, x, y, plotting_args_list, imgsize, finalize, colorbar_loc = None, cb_y = 0):
    plots = []
    fullsize_images = []
    colorbars = []
    try:
        for args in plotting_args_list:
            fig = plt.Figure(figsize=tuple(i/100 for i in imgsize), dpi=100)
            ax = fig.add_subplot(111)
            args['adata'] = self.parent.adata
            args['ax'] = ax
            scatter = spatial_scatter(**args)
            fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
            ax.axis('tight')
            ax.axis('off')

            if colorbar_loc is not None:
                # Extract the colorbar as an image
                # Create a new figure for the colorbar

                colorbar_width = imgsize[0] / 100 *0.5
                colorbar_height = 0.4

                fig_colorbar, ax_colorbar = plt.subplots(figsize=(colorbar_width, colorbar_height), dpi=100)

                for collection in ax.collections:
                    if hasattr(collection, 'get_array'):
                        aax = collection

                cb = fig.colorbar(
                    mappable=aax,   
                    cax=ax_colorbar,
                    orientation='horizontal'
                )
                
                # Extend colorbar
                bot = -0.03
                top = 1.03

                # Upper bound
                xy = np.array([[0, 1], [0, top], [1, top], [1, 1]])
                if colorbar_loc in ["top", "bottom"]:
                    xy = xy[:, ::-1]

                # Make Bezier curve
                curve = [
                    mpath.Path.MOVETO,
                    mpath.Path.CURVE4,
                    mpath.Path.CURVE4,
                    mpath.Path.CURVE4,
                ]

                # Upper patch
                color = cb.cmap(cb.norm(cb._values[-1]))
                patch = patches.PathPatch(
                    mpath.Path(xy, curve),
                    facecolor=color,
                    linewidth=0,
                    antialiased=True,
                    transform=cb.ax.transAxes,
                    clip_on=False,
                )
                cb.ax.add_patch(patch)

                # Lower bound
                xy = np.array([[0, 0], [0, bot], [1, bot], [1, 0]])
                if colorbar_loc in ["top", "bottom"]:
                    xy = xy[:, ::-1]

                # Lower patch
                color = cb.cmap(cb.norm(cb._values[0]))
                patch = patches.PathPatch(
                    mpath.Path(xy, curve),
                    facecolor=color,
                    linewidth=0,
                    antialiased=True,
                    transform=cb.ax.transAxes,
                    clip_on=False,
                )
                cb.ax.add_patch(patch)

                # Hide the outline
                cb.outline.set_visible(False)

                # Save the colorbar to an image file or buffer
                cb_buf = io.BytesIO()
                fig_colorbar.savefig(cb_buf, format='png', bbox_inches='tight', pad_inches=0.1)
                cb_buf.seek(0)
                colorbar_img_fullsize = Image.open(cb_buf)
                cbimg_x, cbimg_y = colorbar_img_fullsize.width, colorbar_img_fullsize.height
                cbratio = min(x / cbimg_x, cb_y / cbimg_y)
                new_cbdimensions = (int(cbimg_x * cbratio), int(cbimg_y * cbratio))
                resized_cbimg = colorbar_img_fullsize.resize(new_cbdimensions, Image.Resampling.LANCZOS)
                colorbars.append(resized_cbimg)


            buf = io.BytesIO()                                                          # Initialize a memory buffer
            fig.savefig(buf, format='png')                                              # save the fig as a png in buffer
            buf.seek(0)                                                                 # rewind to beginning of file
            original_img = Image.open(buf)                                              # set that as original_img
            fullsize_images.append(original_img)                                        # save the original image in a list to retrive it for closeup

            img_x, img_y = original_img.width, original_img.height
            ratio = min(x / img_x, y / img_y)
            new_dimensions = (int(img_x * ratio), int(img_y * ratio))
            resized_img = original_img.resize(new_dimensions, Image.Resampling.LANCZOS)
            plots.append(resized_img)

    except Exception as e:
        print(e)

    finally:
        if finalize == 'viewdata':
            self.current_plot_method = 'squidpy'

        args = (plots, fullsize_images)

        if finalize in ['normalization', 'filterdata'] :
            args = args + (colorbars,)

        self.root.command_queue.put((f'{finalize}_finalize_plotting', args))
        self.root.command_queue.put(('destroy_progress_window', {}))
            
def make_interactive_plots(self, x, y, 
                           plotting_args_list):
    display_imgs = []
    interactive_plots = []

    try:
        for args in plotting_args_list:
            args['gene'] = args.pop('color')
            args['title'] = args['gene']
            args['spot_alpha'] = args.pop('alpha')
            args.pop('colorbar')

            img_interactive , original_plot_img = interactive_gene_plot(x, y, self.parent.adata, **args)
            display_imgs.append(original_plot_img)
            interactive_plots.append(img_interactive)

    except Exception as e:
        print(e)

    finally:
        self.current_plot_method = 'interactive'
        self.root.command_queue.put(('viewdata_finalize_plotting', (display_imgs, interactive_plots)))
        self.root.command_queue.put(('destroy_progress_window', {}))

def interactive_gene_plot(x, y,
                          adata, 
                          gene, 
                          title, 
                          size, 
                          cmap, 
                          spot_alpha, 
                          img_alpha,
                          vmin,
                          vmax):
    hv.Store.register({Circle: CirclePlot}, 'bokeh')
    hv.extension('bokeh')
    data = adata.obs.loc[:, ['imagecol', 'imagerow']]
    data.rename(columns={'imagecol': 'x', 'imagerow': 'y'}, inplace=True)
    tissue_name = next(iter(adata.uns['spatial']))
    spot_diam = (adata.uns['spatial'][tissue_name]['scalefactors']['spot_diameter_fullres'] * adata.uns['spatial'][tissue_name]['scalefactors']['tissue_hires_scalef']) * size

    img = adata.uns['spatial'][tissue_name]['images']['hires']
    width = img.shape[1]
    height = img.shape[0]
    data['y'] = height - data['y']
    data = data.reindex(adata.obs.index)
    if gene == 'in_tissue':
        data[gene] = 1
    else:
        idx = adata.var.index.get_loc(gene)
        filter_data = adata.X[:, idx]
        data[gene] =  np.clip((filter_data.toarray() if sp.isspmatrix_csr(filter_data) else filter_data), vmin, vmax)

    scatter = Circle(data, kdims='x', vdims=['y', gene])
    scatter.opts(height = height, width=width, title=title, color = gene, cmap = cmap, radius = spot_diam / 2, alpha=spot_alpha)
    display_img = hv.RGB(img, bounds=(0,0,width,height))
    display_img.opts(alpha=img_alpha)
    img_interactive = (display_img * scatter).opts(opts.Circle(tools=['hover']))

    ratio = min(x / width, y / height)
    new_dimensions = (int(width * ratio), int(height * ratio))
    fig = plt.figure(figsize=(tuple(i/100 for i in new_dimensions)), dpi=100)
    ax = fig.add_subplot(111)
    ax.imshow(img, extent=[0, width, 0, height], alpha = img_alpha)
    ax.scatter(data['x'], data['y'], s=(spot_diam / 2)*ratio, c=data[gene], alpha=spot_alpha, cmap=cmap)
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    ax.axis('tight')
    ax.axis('off')
    buf = io.BytesIO()                                              
    fig.savefig(buf, format='png')                                        
    buf.seek(0)                                                              
    original_plot_img = Image.open(buf)

    return img_interactive, original_plot_img

def make_datashader_plots(self, x, y, 
                          plotting_args_list):
    plots = []
    fullsize_images = []
    try:
        for args in plotting_args_list:
            args['gene'] = args.pop('color')
            args['title'] = args['gene']
            args['spot_alpha'] = args.pop('alpha')
            args.pop('colorbar')

            original_img = datashader_geneplot(self.parent.adata, **args)
            fullsize_images.append(original_img)

            print(original_img)
            img_x, img_y = original_img.size
            ratio = min(x / img_x, y / img_y)
            new_dimensions = (int(img_x * ratio), int(img_y * ratio))
            resized_img = original_img.resize(new_dimensions, Image.Resampling.LANCZOS)
            plots.append(resized_img)

    except Exception as e:
        print(e)

    finally:
        self.current_plot_method = 'datashader'
        self.root.command_queue.put(('viewdata_finalize_plotting', (plots, fullsize_images)))
        self.root.command_queue.put(('destroy_progress_window', {}))

def datashader_geneplot(adata, 
                        gene,
                        size, 
                        cmap, 
                        spot_alpha, 
                        img_alpha,
                        vmin,
                        vmax,
                        title):
    
    data = adata.obs.loc[:, ['imagecol', 'imagerow']]
    data.rename(columns={'imagecol': 'x', 'imagerow': 'y'}, inplace=True)
    tissue_name = next(iter(adata.uns['spatial']))
    spot_diam = (adata.uns['spatial'][tissue_name]['scalefactors']['spot_diameter_fullres'] * adata.uns['spatial'][tissue_name]['scalefactors']['tissue_hires_scalef']) * size

    timg_arr = adata.uns['spatial'][tissue_name]['images']['hires'] * 255
    white_arr = np.full_like(timg_arr, 255)
    adjusted_img_arr = timg_arr * img_alpha + white_arr * (1 - img_alpha)
    adjusted_img_arr = np.clip(adjusted_img_arr, 0, 255).astype(np.uint8)
    timg = Image.fromarray(adjusted_img_arr)    

    width = timg_arr.shape[1]
    height = timg_arr.shape[0]
    data['y'] = height - data['y']
    data = data.reindex(adata.obs.index)
    if gene == 'in_tissue':
        data[gene] = 1
    else:
        idx = adata.var.index.get_loc(gene)
        filter_data = adata.X[:, idx]
        data[gene] =  np.clip((filter_data.toarray() if sp.isspmatrix_csr(filter_data) else filter_data), vmin, vmax)

    canvas = ds.Canvas(plot_width=width, plot_height=height)
    agg = canvas.points(data, 'x', 'y', agg=ds.mean(gene))
    cmap = colormaps[cmap]
    img = tf.shade(agg, cmap=cmap, alpha = spot_alpha*255, how='linear')
    spread_img = tf.spread(img, px=int((spot_diam/2)*size)).to_pil()
    timg.paste(spread_img, (0, 0), spread_img)
    return timg