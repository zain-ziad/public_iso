#       Import Packages     #
import customtkinter as ctk
from core.assets.themecolors import COLORS
from core.customwidgets import *
from core.backend.viewdata_backend import *
import threading
import numpy as np
from tkinter import PhotoImage
from core.validatecommands import *
import math
from bokeh.plotting import show


#       ViewData Class     #
class ViewData(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent, fg_color=COLORS['background'])
        self.parent = parent
        self.root = self.parent.parent
        self.columnconfigure(1, weight=1, uniform=1)
        self.columnconfigure(0, weight=3, uniform=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=1)

        self.selected_plot = ctk.StringVar(value='None')
        self.plot_method = ctk.StringVar(value='Squidpy')
        self.current_plot_method = 'squidpy'
        self.current_plot_toplevels = []
        self.plots_list = []
        self.plots_frame = {}                   #Store the frames for different plots
        self.celltypelist = []
        self.currentgenes = []

        self.initialize()

    def initialize(self):
        #   Functions   #
        def plotentry_on_enter(*args):
            if len(self.plots_list) == 6:
                self.root.lockAppWindow()
                response = self.root.build_msgmodal(msgbox_args = {'height' : 372,
                                                                    'no_button':False,
                                                                    'cancel_button':False,
                                                                    'yes_text': 'Continue',
                                                                    'title':'Warning',
                                                                    'icon':'warning',
                                                                    'title_color':COLORS['warning'],
                                                                    'message' : 'Displaying more than 6 plots simultaneously \n may increase processing time.'})
                self.root.releaseAppWindow()
            
            name = self.plot_name_entry.get()
            if name == '':
                num = len(self.plots_list) + 1
                name = "Plot " + str(num)
            
            original_name = name
            counter = 1
            while name in self.plots_list:
                name = f"{original_name}_{counter}"
                counter += 1

            self.plot_name_entry.delete(0, tk.END)
            self.plot_name_entry.configure(placeholder_text='Enter Name...')

            self.plots_list.append(name)
            if self.plot_optionmenu.cget('state') == 'disabled':
                self.plot_optionmenu.configure(state = 'normal')

            self.plot_optionmenu.configure(values = self.plots_list)
            self.plots_frame[name] = self.Plot(parent=self.options_subframe)
            self.selected_plot.set(name) 
            self.focus()
        
        def plot_trash(*args):
            name = self.selected_plot.get()     #Get name of selected plot in the option menu
            if name == 'None':                  #If its None, do nothing
                return
            self.root.lockAppWindow()
            response = self.root.build_msgmodal(msgbox_args = {'height' : 342,
                                                                    'no_button':True,
                                                                    'cancel_button':False,
                                                                    'yes_text': 'Yes',
                                                                    'title':'Confirmation',
                                                                    'icon':'info',
                                                                    'title_color':COLORS['primary 1'],
                                                                    'message' : f"Do you want to permanently delete the plot '{self.selected_plot.get()}' ?"})
            self.root.releaseAppWindow()
            if response == 'No':
                return
            selection_index = self.plots_list.index(name) if self.plots_list.index(name) != len(self.plots_list)-1 else  self.plots_list.index(name) -1   #First, we need to know which plot we'll select once the current one is removed. We selected the plot previous to it. 
            self.plots_list.remove(name)                        #Remove the current plot.
            self.plot_optionmenu.configure(values = self.plots_list)    #Update the option menu list.
            self.selected_plot.set(self.plots_list[selection_index] if selection_index >= 0  else 'None')    #Update the selected plot to the index we calculated. If its already 0, use 'None'
        
        def switch_plot(*args):
            to = self.selected_plot.get()

            if hasattr(self, 'subclass'):
                self.subclass.grid_forget()

            if to == 'None':
                return
            
            self.subclass = self.plots_frame[to]
            self.subclass.grid(row = 8, column = 0, columnspan = 2, sticky = 'we')
            self.subclass.colormap_widget.resize()

        #   Top Bar  #
        topbar_frame = ctk.CTkFrame(self, fg_color='transparent')
        customTkLabel(topbar_frame, text='View Data', font=('Poppins Medium', 40), height=54, textMoveY= -20, text_color=COLORS['primary 1'], bg_color=COLORS['background']).pack(anchor = 'w', side = 'top')
        ctk.CTkLabel(topbar_frame, textvariable = self.parent.dataset_displayname, font = ('Poppins Medium', 18), text_color=COLORS['neutral 2']).pack(anchor = 'w', side = 'bottom', pady = 6)
        topbar_frame.grid(row = 0, column = 0, columnspan = 2, sticky = 'NSWE')

        #   No Plot    #
        self.no_plot_frame = ctk.CTkFrame(self, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])

        #   Options     #
        options_frame = ctk.CTkFrame(self, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])

        #   Create the line #
        line_canvas = tk.Canvas(options_frame, height = 3, bg='white', highlightthickness=0)
        line_canvas.create_line((0,0,999,1), dash=(4,8), fill = COLORS['neutral 2'])

        #   Create the buttons for plotting #
        buttons_frame = ctk.CTkFrame(options_frame, fg_color='white', corner_radius=0, height=60)
        self.export_button = customCTkOptionMenu(buttons_frame, 
                            values= ['All Plots'],
                            text = 'Export', 
                            font=('Poppins Semibold', 18),
                            fg_color=COLORS['neutral 4'],
                            button_color=COLORS['neutral 4'],
                            button_hover_color=COLORS['neutral 4'],
                            add_border= False,
                            corner_radius=16,
                            button_fg_color = COLORS['neutral 4'])
        self.export_button.pack(pady = 12, side = 'left', fill = 'y', expand = True)
        self.refresh_button = ctk.CTkButton(buttons_frame,
                                            text='Refresh',
                                            font=('Poppins Semibold', 18),
                                            corner_radius=16,
                                            fg_color=COLORS['primary 1'],
                                            hover_color=COLORS['primary 1 hover'],
                                            height=64,
                                            text_color= COLORS['white'],
                                            command=self.initiate_plotting)
        self.refresh_button.pack(pady = 12, side = 'right', fill = 'y', expand = True)


        #   Options Subframe    #
        self.options_subframe = ctk.CTkScrollableFrame(options_frame, 
                                                  fg_color=COLORS['white'],
                                                  corner_radius=0, 
                                                  scrollbar_button_color=COLORS['neutral 3'],
                                                  scrollbar_button_hover_color=COLORS['neutral 3'])
        self.options_subframe.columnconfigure(1, weight=1, uniform=1)
        self.options_subframe.columnconfigure(0, weight=4, uniform=1)

        ctk.CTkLabel(self.options_subframe,
                     text='Manage Plots',
                     text_color= COLORS['primary 1'],
                     font = ('Poppins Medium', 20)).grid(row = 0, column = 0, columnspan = 2, sticky = 'w', pady = (0, 40))
        
        ctk.CTkLabel(self.options_subframe,
                     text='Plot Method',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 1, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
        
        self.plotmethod_optionmenu = customCTkOptionMenu(self.options_subframe,
                            variable = self.plot_method,
                            height=56, image_padx = 16,
                            values = ['Squidpy', 'Interactive', 'Datashader'])
    
        self.plotmethod_optionmenu.grid(row = 2, column = 0, columnspan = 2, sticky = 'we', padx = (0, 8), pady = (0, 40))
        
        ctk.CTkLabel(self.options_subframe,
                     text='Add New Plot',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 3, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))

        self.plot_name_entry = ctk.CTkEntry(self.options_subframe,
                     placeholder_text='Enter Name...',
                     border_width=1,
                     border_color=COLORS['border'],
                     corner_radius=16,
                     height=56,
                     font = ('Poppins Medium', 18),
                     text_color= COLORS['neutral 1'],
                     fg_color=COLORS['white'],
                     placeholder_text_color=COLORS['neutral 2'])
        self.plot_name_entry.grid(row = 4, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))
        
        ctk.CTkButton(self.options_subframe,
                      text= '+',
                      text_color=COLORS['primary 1'],
                      height=56,
                      width=56,
                      fg_color=COLORS['primary 2'],
                      corner_radius=16,
                      hover_color=COLORS['primary 2 hover'],
                      font=('Poppins', 20),
                      command=plotentry_on_enter).grid(row = 4, column = 1, sticky = 'e', pady = (0, 40), padx = (0, 8))
        
        ctk.CTkLabel(self.options_subframe,
                     text='Select Plot',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 5, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
        
        self.plot_optionmenu = customCTkOptionMenu(self.options_subframe,
                            variable = self.selected_plot,
                            height=56, image_padx = 16, state = 'disabled')
        
        self.plot_optionmenu.grid(row = 6, column = 0, sticky = 'we', padx = (0, 8))

        ctk.CTkButton(self.options_subframe,
                      text= '',
                      image= ctk.CTkImage(light_image = Image.open('app\\core\\assets\\trash.png'), size=(24,24)),
                      height=56,
                      width=56,
                      fg_color=COLORS['error 2'],
                      corner_radius=16,
                      hover_color=COLORS['error 2 hover'],
                      command=plot_trash).grid(row = 6, column = 1, sticky = 'e', padx = (0, 8))    
        
        ctk.CTkLabel(self.options_subframe,
                     text='Plotting Options',
                     text_color= COLORS['primary 1'],
                     font = ('Poppins Medium', 20)).grid(row = 7, column = 0, columnspan = 2, sticky = 'w', pady = (64, 40))

        self.plot_name_entry.bind("<Return>", plotentry_on_enter)
        self.selected_plot.trace_add('write', switch_plot)

        self.no_plot_frame.grid(row = 1, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')
        options_frame.grid(row = 1, column = 1, pady = (24, 0), sticky='NSWE')
        buttons_frame.pack(side = 'bottom', fill = 'x', padx = 18, pady = (0, 10))
        line_canvas.pack(side = 'bottom', fill = 'x', padx = 2)
        self.options_subframe.pack(side = 'top', fill = 'both', expand = True, padx = (24,12), pady = (40, 2))

    def plot_closeup(self, title, index):
        img = self.fullsize_plots[index]

        if self.current_plot_method in ['squidpy', 'datashader']:
            for toplevel in self.current_plot_toplevels:
                if toplevel.name == title:
                    toplevel.deiconify()
                    toplevel.lift()
                    toplevel.focus_force()
                    return
            
            toplevel = PlotTopLevel(name = title, parent= self, geometry=self.root.geometry(), img = img)
            self.current_plot_toplevels.append(toplevel)

        elif self.current_plot_method == 'interactive':
            show(hv.render(img))

    def initiate_plotting(self):
        #Plotting is done in three steps: Initiate, Process, and Finalize. 
        #This is the Initiate function.
        #This function cleans the area of any existing plots, removes the "No plots displayed" message, builds containers for the plots and sends the process request to the Plot subclass.
    
        self.root.lockAppWindow()                                                                                              #Lock the app window
        n = len(self.plots_list)                                                                                               #How many plots to make
        self.listofcanvases = []                                                                                               #A list to store the plot canvases of all the containers. 

        if hasattr(self, 'container_parent'):                                                                                  #Check and destroy if the container frame already exists
            self.container_parent.destroy()

        if n == 0:                                                                                                             #If the user hasen't added any plots. 
            response = self.root.build_msgmodal(msgbox_args = {'height' : 350,
                                                                'no_button':False,
                                                                'cancel_button':False,
                                                                'yes_text': 'Continue',
                                                                'title':'Error',
                                                                'icon':'fail',
                                                                'title_color':COLORS['error'],
                                                                'message' : 'No plots added.'})
            
            self.no_plot_frame.grid(row = 1, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')
            self.root.releaseAppWindow()
            return
        
        if self.no_plot_frame.winfo_exists():                                                                                   #Check and destroy the no plots frame
            self.no_plot_frame.grid_forget()

        self.container_parent = ctk.CTkFrame(self, fg_color=COLORS['background'])                                               #Make the parent container frame
        self.container_parent.grid(row = 1, column = 0, padx = (0, 24), pady = (16, 0), sticky='NSWE')
        grid_size = math.ceil(math.sqrt(n))                                                                                     #Calculate the grid size based on the number of plots

        #Here, we set the weight of the rows and the columns of the parent container frame so each container is of equal size. 
        for col in range(grid_size):
            self.container_parent.columnconfigure(col, weight=1)
        num_rows = (len(self.plots_list) + grid_size - 1) // grid_size
        for row in range(num_rows):
            self.container_parent.rowconfigure(row, weight=1)


        #Now make all the individual containers. 
        for i, plot in enumerate(self.plots_list):
            row = i // grid_size                                                                                                # Integer division to determine the row index
            col = i % grid_size                                                                                                 # Modulus to determine the column index
            make_container(self, title=plot, row= row, col = col, index = i, frame_grid_ypad = (6, 0))
      

        self.root.update_idletasks()                                                                                            #Update the containers on screen. 
        y = self.listofcanvases[0].winfo_height()                                                                               #Get a container height. 
        x = self.listofcanvases[0].winfo_width()                                                                                #Get a container width. 
        plotting_args_list = []                                                                                                 #Make a list to store all of the plotting arguments of different Plot subclasses. 


        #Get the plotting arguments from each Plot subclass. 
        for key, frame in self.plots_frame.items():
            plotting_args_list.append(frame.build_plotting_args())

        #Build an inderminate progress modal. 
        self.root.build_IndProgressModal()

        #Start a secondary thread based on the method selected. (Step 2: Process)
        if self.plot_method.get() == 'Squidpy':
            secondary_thread = threading.Thread(target=make_squidpy_plots, args=[self, x, y, plotting_args_list, self.parent.imgsize, 'viewdata'])
            secondary_thread.start()
        elif self.plot_method.get() == 'Interactive':
            secondary_thread = threading.Thread(target=make_interactive_plots, args=[self, x, y, plotting_args_list])
            secondary_thread.start()            
        elif self.plot_method.get() == 'Datashader':
            secondary_thread = threading.Thread(target=make_datashader_plots, args=[self, x, y, plotting_args_list])
            secondary_thread.start()

    def finalize_plotting(self, plots, fullsize_plots):
        self.plot_images = []
        self.fullsize_plots = fullsize_plots

        for i in range(len(self.listofcanvases)):
            plot_image = ImageTk.PhotoImage(plots[i])  # Convert the i-th image to a PhotoImage
            self.plot_images.append(plot_image)  # Keep a reference to the PhotoImage
            canvas = self.listofcanvases[i]

            # Get canvas width and height
            canvas_width = canvas.winfo_width()
            canvas_height = canvas.winfo_height()

            # Calculate coordinates to place the image at the center of the canvas
            x = (canvas_width - plot_image.width()) // 2
            y = (canvas_height - plot_image.height()) // 2

            canvas.create_image(x, y, image=self.plot_images[i], anchor='nw')

    def add_gene_buttons(self, gene):
        for key, frame in self.plots_frame.items():
            frame.make_gene_button(gene = gene)
        self.currentgenes.append(gene)

    class Plot(ctk.CTkFrame):
        def __init__(self, parent):
            super().__init__(parent, fg_color='transparent')
            self.parent = parent                    #options_subframe
            self.viewdata = self.parent.master.master.master.master           #viewdata object
            self.mainwindow = self.viewdata.master
            self.root = self.viewdata.root          #root
            self.vcmdFloat = (self.register(callbackFloat))

            #Initilizing some variables we will use in this class
            self.original_img = None
            self.vmin = ctk.StringVar(value='Auto')
            self.vmax = ctk.StringVar(value='Auto')
            self.cmap = ctk.StringVar(value='Viridis')
            self.tissue_opacity = ctk.DoubleVar(value=1.0)
            self.spot_opacity = ctk.DoubleVar(value=1.0)
            self.spot_size = ctk.DoubleVar(value=1.0)
            self.selected_gene = ctk.StringVar(value='in_tissue')
            self.selected_celltype = ctk.StringVar(value='')
            self.plot = None

            self.options()
            self.make_celltype_buttons(celltypelist=self.viewdata.celltypelist)

            for gene in self.viewdata.currentgenes:
                self.make_gene_button(gene = gene)

        def reset_vmin_vmax(self, *args):
            self.vmin.set('Auto')
            self.vmax.set('Auto')

        def make_celltype_buttons(self, celltypelist):
            pass
        
        def make_gene_button(self, gene):
            button = ctk.CTkRadioButton(self.gene_frame, text=gene, 
                                        font=('Poppins Medium', 18), 
                                        variable = self.selected_gene, 
                                        value=gene,
                                        fg_color=COLORS['primary 1'],
                                        text_color=COLORS['neutral 1'],
                                        radiobutton_height=24,
                                        radiobutton_width=24,
                                        border_width_unchecked = 2,
                                        border_width_checked = 6,
                                        border_color=COLORS['neutral 2'],
                                        hover_color=COLORS['primary 1 hover'],
                                        command=self.reset_vmin_vmax)
            button.pack(side='top', anchor= 'nw', pady = (4,0), padx = (12,0))

        def gene_selection(self, *args):
            current_gene = self.gene_name_entry.get()
            if "," in current_gene:
                error = False
                currentlist = [element.strip() for element in current_gene.split(',') if element != '']
                for gene in currentlist:
                    if gene in self.mainwindow.gene_names and gene not in self.viewdata.currentgenes:
                        self.viewdata.add_gene_buttons(gene = gene)
                    else:
                        error = True
                self.gene_name_entry.delete(0, ctk.END)
                if error:
                    self.root.lockAppWindow()
                    response = self.root.build_msgmodal(msgbox_args = {'height' : 384,
                                                                        'no_button':False,
                                                                        'cancel_button':False,
                                                                        'yes_text': 'Continue',
                                                                        'title':'Addition Partially Successful',
                                                                        'icon':'fail',
                                                                        'title_color':COLORS['error'],
                                                                        'message' : 'Some genes were not added. \n Check if they exist in the dataset or \n if they have already been added.'})
                    self.root.releaseAppWindow()
                
            elif current_gene in self.mainwindow.gene_names:
                if current_gene not in self.viewdata.currentgenes:
                    self.viewdata.add_gene_buttons(gene = current_gene)
                    self.gene_name_entry.delete(0, ctk.END)
                else:
                    self.root.lockAppWindow()
                    response = self.root.build_msgmodal(msgbox_args = {'height' : 360,
                                                                        'no_button':False,
                                                                        'cancel_button':False,
                                                                        'yes_text': 'Continue',
                                                                        'title':'Gene Already in List',
                                                                        'icon':'fail',
                                                                        'title_color':COLORS['error'],
                                                                        'message' : f'Gene "{current_gene}" already in the list'})
                    self.root.releaseAppWindow()
                    self.gene_name_entry.delete(0, ctk.END)
            else:
                self.root.lockAppWindow()
                response = self.root.build_msgmodal(msgbox_args = {'height' : 360,
                                                                    'no_button':False,
                                                                    'cancel_button':False,
                                                                    'yes_text': 'Continue',
                                                                    'title':'Gene Not Found',
                                                                    'icon':'fail',
                                                                    'title_color':COLORS['error'],
                                                                    'message' : f'Gene not found in the dataset.'})
                self.root.releaseAppWindow()
                self.gene_name_entry.delete(0, ctk.END)

        def options(self):
            mainframe = ctk.CTkFrame(self, fg_color=COLORS['white'], corner_radius=0)
            mainframe.columnconfigure(1, weight=1, uniform=1)
            mainframe.columnconfigure(0, weight=3, uniform=1)

            #Tissue Opacity
            def update_to_entry(*args):
                to_entry.delete(0, tk.END)
                to_entry.insert(0, str(round(self.tissue_opacity.get(),2)))

            def focusout_to_entry(*args):
                if to_entry.get() == '' or float(to_entry.get()) > 1:
                    to_entry.delete(0, tk.END)
                    to_entry.insert(0, '1.0')
                    self.tissue_opacity.set(1.0)
                else:
                    self.tissue_opacity.set(float(to_entry.get()))
                self.focus()

            ctk.CTkLabel(mainframe,
                        text='Tissue Opacity',
                        text_color= COLORS['neutral 2'],
                        font = ('Poppins Medium', 16)).grid(row = 1, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
            
            to_entry = ctk.CTkEntry(mainframe,
                     border_width=1,
                     border_color=COLORS['border'],
                     corner_radius=16,
                     height=56,
                     font = ('Poppins Medium', 18),
                     text_color= COLORS['neutral 1'],
                     fg_color=COLORS['white'],
                     validatecommand=(self.vcmdFloat, '%P'),
                     validate='all',
                     placeholder_text= '1.0',
                     placeholder_text_color= COLORS['neutral 1']
                    )
            
            to_slider = ctk.CTkSlider(mainframe, from_ = 0, 
                                      to = 1, number_of_steps=100, 
                                      variable=self.tissue_opacity,
                                      progress_color = COLORS['primary 1'],
                                      button_color = COLORS['primary 1'],
                                      fg_color=COLORS['background'],
                                      button_hover_color=COLORS['primary 1 hover'],
                                      command=update_to_entry)
            
            #Spot Opacity
            def update_so_entry(*args):
                so_entry.delete(0, tk.END)
                so_entry.insert(0, str(round(self.spot_opacity.get(),2)))

            def focusout_so_entry(*args):
                if so_entry.get() == '' or float(so_entry.get()) > 1:
                    so_entry.delete(0, tk.END)
                    so_entry.insert(0, '1.0')
                    self.spot_opacity.set(1.0)
                else:
                    self.spot_opacity.set(float(so_entry.get()))
                self.focus()

            ctk.CTkLabel(mainframe,
                        text='Spot Opacity',
                        text_color= COLORS['neutral 2'],
                        font = ('Poppins Medium', 16)).grid(row = 4, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
            
            so_entry = ctk.CTkEntry(mainframe,
                     border_width=1,
                     border_color=COLORS['border'],
                     corner_radius=16,
                     height=56,
                     font = ('Poppins Medium', 18),
                     text_color= COLORS['neutral 1'],
                     fg_color=COLORS['white'],
                     validatecommand=(self.vcmdFloat, '%P'),
                     validate='all',
                     placeholder_text= '1.0',
                     placeholder_text_color= COLORS['neutral 1']
                    )
            
            so_slider = ctk.CTkSlider(mainframe, from_ = 0, 
                                      to = 1, number_of_steps=100, 
                                      variable=self.spot_opacity,
                                      progress_color = COLORS['primary 1'],
                                      button_color = COLORS['primary 1'],
                                      fg_color=COLORS['background'],
                                      button_hover_color=COLORS['primary 1 hover'],
                                      command=update_so_entry)
            
            #Spot Size
            def update_ss_entry(*args):
                ss_entry.delete(0, tk.END)
                ss_entry.insert(0, str(round(self.spot_size.get(),2)))

            def focusout_ss_entry(*args):
                if ss_entry.get() == '' or float(ss_entry.get()) > 1:
                    ss_entry.delete(0, tk.END)
                    ss_entry.insert(0, '1.0')
                    self.spot_size.set(1.0)
                else:
                    self.spot_size.set(float(ss_entry.get()))
                self.focus()

            ctk.CTkLabel(mainframe,
                        text='Spot Size',
                        text_color= COLORS['neutral 2'],
                        font = ('Poppins Medium', 16)).grid(row = 7, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
            
            ss_entry = ctk.CTkEntry(mainframe,
                     border_width=1,
                     border_color=COLORS['border'],
                     corner_radius=16,
                     height=56,
                     font = ('Poppins Medium', 18),
                     text_color= COLORS['neutral 1'],
                     fg_color=COLORS['white'],
                     validatecommand=(self.vcmdFloat, '%P'),
                     validate='all',
                     placeholder_text= '1.0',
                     placeholder_text_color= COLORS['neutral 1']
                    )
            
            ss_slider = ctk.CTkSlider(mainframe, from_ = 0, 
                                      to = 2, number_of_steps=200, 
                                      variable=self.spot_size,
                                      progress_color = COLORS['primary 1'],
                                      button_color = COLORS['primary 1'],
                                      fg_color=COLORS['background'],
                                      button_hover_color=COLORS['primary 1 hover'],
                                      command=update_ss_entry)
            

            mainframe.pack(fill = 'both', expand = True)
            to_entry.grid(row = 2, column = 0, columnspan = 2, sticky = 'ew', pady = (0, 16), padx = (0, 8))
            to_slider.grid(row = 3, column = 0, columnspan = 2, sticky = 'ew', pady = (0, 16), padx = (0, 8))
            to_entry.bind('<Return>', focusout_to_entry)
            to_entry.bind("<FocusOut>", focusout_to_entry)

            so_entry.grid(row = 5, column = 0, columnspan = 2, sticky = 'ew', pady = (0, 16), padx = (0, 8))
            so_slider.grid(row = 6, column = 0, columnspan = 2, sticky = 'ew', pady = (0, 16), padx = (0, 8))
            so_entry.bind('<Return>', focusout_so_entry)
            so_entry.bind("<FocusOut>", focusout_so_entry)

            ss_entry.grid(row = 8, column = 0, columnspan = 2, sticky = 'ew', pady = (0, 16), padx = (0, 8))
            ss_slider.grid(row = 9, column = 0, columnspan = 2, sticky = 'ew', pady = (0, 16), padx = (0, 8))
            ss_entry.bind('<Return>', focusout_ss_entry)
            ss_entry.bind("<FocusOut>", focusout_ss_entry)


            #Colormap
            ctk.CTkLabel(mainframe,
                        text='Color Map',
                        text_color= COLORS['neutral 2'],
                        font = ('Poppins Medium', 16)).grid(row = 10, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))            

            customCTkOptionMenu(mainframe,
                                variable = self.cmap,
                                values= ['Viridis', 'Plasma', 'Inferno', 'Magma', 'Cividis'],
                                height=56, image_padx = 16).grid(row = 11, column = 0, columnspan = 2, sticky = 'we', pady = (0, 16), padx = (0, 5))
        
            self.colormap_widget = ColorMapWidget(mainframe, colormap=self.cmap, root = self.root)
            self.colormap_widget.grid(row = 12, column = 0, columnspan = 2, sticky = 'we', pady = (0, 16), padx = (0, 4))

            self.cl_entry_min = ctk.CTkEntry(mainframe,
            border_width=1,
            border_color=COLORS['border'],
            corner_radius=16,
            height=56,
            width=82,
            font = ('Poppins Medium', 18),
            text_color= COLORS['neutral 1'],
            fg_color=COLORS['white'],
            textvariable=self.vmin
            ).grid(row = 13, column = 0, columnspan = 1, sticky = 'w', pady = (0, 16))

            self.cl_entry_max = ctk.CTkEntry(mainframe,
            border_width=1,
            border_color=COLORS['border'],
            corner_radius=16,
            height=56,
            width=82,
            font = ('Poppins Medium', 18),
            text_color= COLORS['neutral 1'],
            fg_color=COLORS['white'],
            textvariable=self.vmax
            ).grid(row = 13, column = 1, columnspan = 1, sticky = 'e', pady = (0, 16))


            #Gene Frame
            ctk.CTkLabel(mainframe,
                        text='Gene',
                        text_color= COLORS['neutral 2'],
                        font = ('Poppins Medium', 16)).grid(row = 14, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))

            self.gene_name_entry = ctk.CTkEntry(mainframe,
                        placeholder_text='Enter Gene Name',
                        border_width=1,
                        border_color=COLORS['border'],
                        corner_radius=16,
                        height=56,
                        font = ('Poppins Medium', 18),
                        text_color= COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        placeholder_text_color=COLORS['neutral 2'])
            self.gene_name_entry.grid(row = 15, column = 0, sticky = 'we', pady = (0, 16))

            ctk.CTkButton(mainframe,
                        text= '+',
                        text_color=COLORS['primary 1'],
                        height=56,
                        width=56,
                        fg_color=COLORS['primary 2'],
                        corner_radius=16,
                        hover_color=COLORS['primary 2 hover'],
                        font=('Poppins', 20),
                        command=self.gene_selection
                        ).grid(row = 15, column = 1, sticky = 'e', pady = (0, 16), padx = (0, 8))
            
            self.gene_frame = ctk.CTkScrollableFrame(mainframe,
                                           fg_color=COLORS['white'],
                                           border_color=COLORS['border'],
                                           border_width=1,
                                           corner_radius=16,
                                           scrollbar_button_color=COLORS['neutral 3'],
                                           height=150)
            self.gene_frame.grid(row = 16, column = 0, columnspan = 2, sticky = 'we', pady = (0, 18))
            self.gene_frame._scrollbar.configure(height=0)
            self.gene_frame._scrollbar.grid_configure(padx = 2)

            ctk.CTkLabel(mainframe,
                        text='Deconvoluted Cell Types',
                        text_color= COLORS['neutral 2'],
                        font = ('Poppins Medium', 16)).grid(row = 17, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
            
            self.decon_disabled_label = ctk.CTkLabel(mainframe,
                                                     image= ctk.CTkImage(light_image = Image.open("app\\core\\assets\\viewdata\\warning.png"), size=(22,19)),
                                                     text="   Not currently available.",
                                                     font = ('Poppins Medium', 14),
                                                     text_color= COLORS['neutral 2'],
                                                     compound='left')
            self.decon_disabled_label.grid(row = 18, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))


            #   Bindings    #
            self.gene_name_entry.bind("<Return>", lambda _: self.gene_selection())

        def build_plotting_args(self):
            #This function turns the selections into plotting parameters. 

            #Plot only gene
            if self.selected_gene.get() in self.mainwindow.gene_names + ['in_tissue'] and self.selected_celltype.get() == '':
                if self.selected_gene.get() in self.mainwindow.gene_names:                                             
                    gene_idx = self.viewdata.parent.adata.var_names.get_loc(self.selected_gene.get())
                    gene_counts = self.viewdata.parent.adata.X[:, gene_idx]
                    try:
                        vmin = float(self.vmin.get())
                    except:
                        vmin = (round(gene_counts.min()))
                        self.vmin.set(vmin)
                    try:
                        vmax = float(self.vmax.get())
                    except:
                        vmax = (round(gene_counts.max()))
                        self.vmax.set(vmax)
                else:
                    vmin = None
                    vmax = None
                    self.vmin.set('Auto')
                    self.vmax.set('Auto')

                plotting_args = {
                    'color': self.selected_gene.get(),
                    'img_alpha': self.tissue_opacity.get(),
                    'alpha': self.spot_opacity.get(),
                    'size': self.spot_size.get(),
                    'vmin' : vmin, 
                    'vmax' : vmax, 
                    'cmap' : self.cmap.get().lower(),
                    'colorbar' : False
                }
                if self.mainwindow.xenium_imgpath == None and self.mainwindow.isXenium:
                    plotting_args['img'] = None

            #Add other things to plot here

            else:
                self.root.lockAppWindow()
                response = self.root.build_msgmodal(msgbox_args = {'height' : 384,
                                                                    'no_button':False,
                                                                    'cancel_button':False,
                                                                    'yes_text': 'Continue',
                                                                    'title':'Gene Not Found in Dataset',
                                                                    'icon':'fail',
                                                                    'title_color':COLORS['error'],
                                                                    'message' : 'The gene you have selected could not be located within the current dataset. \n This issue typically arises when the specified gene has been excluded \n during the data filtration process.'})
                self.root.releaseAppWindow()
                vmin = None
                vmax = None
                self.vmin.set('Auto')
                self.vmax.set('Auto')
                return
            
            return plotting_args
