#       Import Packages     #
import customtkinter as ctk
from core.assets.themecolors import COLORS
from core.customwidgets import *
from core.backend.viewdata_backend import *
from core.backend.normalization_backend import *
import numpy as np
from tkinter import PhotoImage
from core.validatecommands import *
import math
from bokeh.plotting import show
import threading

#       Normalization Class     #
class Normalization(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent, fg_color=COLORS['background'])
        self.parent = parent
        self.root = self.parent.parent
        self.columnconfigure(1, weight=1, uniform=1)
        self.columnconfigure(0, weight=3, uniform=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=1)

        #   Init Variables   #
        self.normalization_method = ctk.StringVar(value='None')
        self.do_log1p = ctk.BooleanVar(value=False)
        self.do_scale = ctk.BooleanVar(value=False)
        self.plot_type = ctk.StringVar(value='Spatial Plot')
        self.state = ctk.StringVar(value='Side by side')
        self.plot_enabled = False
        self.current_plot_toplevels = []

        #   Build Main Components   #
        self.initialize()

    def initialize(self):
        #   Functions   #


        #   Top Bar  #
        topbar_frame = ctk.CTkFrame(self, fg_color='transparent')
        customTkLabel(topbar_frame, text='Normalization', font=('Poppins Medium', 40), height=54, textMoveY= -20, text_color=COLORS['primary 1'], bg_color=COLORS['background']).pack(anchor = 'w', side = 'top')
        ctk.CTkLabel(topbar_frame, textvariable = self.parent.dataset_displayname, font = ('Poppins Medium', 18), text_color=COLORS['neutral 2']).pack(anchor = 'w', side = 'bottom', pady = 6)
        topbar_frame.grid(row = 0, column = 0, columnspan = 2, sticky = 'NSWE')

        #   No Plot    #
        self.no_plot_frame = ctk.CTkFrame(self, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])
        norm_icon = ctk.CTkImage(light_image = Image.open('app\\core\\assets\\normalization\\normalization.png'),size=(48,48))
        self.no_plot_frame_centered = ctk.CTkFrame(self.no_plot_frame, fg_color='transparent')
        ctk.CTkLabel(self.no_plot_frame_centered, text="", image=norm_icon).pack(pady = (0, 16))
        self.disabled_text = ctk.CTkLabel(self.no_plot_frame_centered, text="Gene plotting available after dataset normalization.", font = ('Poppins Medium', 16), text_color=COLORS['neutral 2'])
        self.disabled_text.pack()
        self.no_plot_frame_centered.pack(expand = True)

        #   Options     #
        options_frame = ctk.CTkFrame(self, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])

        #   Create the line #
        line_canvas = tk.Canvas(options_frame, height = 3, bg='white', highlightthickness=0)
        line_canvas.create_line((0,0,999,1), dash=(4,8), fill = COLORS['neutral 2'])

        #   Create the buttons for plotting #
        buttons_frame = ctk.CTkFrame(options_frame, fg_color='white', corner_radius=0, height=60)
        self.export_button = ctk.CTkButton(buttons_frame, 
                            text = 'Export', 
                            font=('Poppins Semibold', 18),
                            fg_color=COLORS['neutral 4'],
                            hover_color=COLORS['neutral 4'],
                            corner_radius=16,
                            state='disabled')
        self.export_button.pack(pady = 12, side = 'left', fill = 'y', expand = True)
        self.refresh_button = ctk.CTkButton(buttons_frame,
                                            text='Refresh',
                                            font=('Poppins Semibold', 18),
                                            corner_radius=16,
                                            fg_color=COLORS['primary 1'],
                                            hover_color=COLORS['primary 1 hover'],
                                            height=64,
                                            text_color= COLORS['white'],
                                            command=self.init_plot,
                                            state='disabled')
        self.refresh_button.pack(pady = 12, side = 'right', fill = 'y', expand = True)

        #   Options Subframe    #
        self.options_subframe = ctk.CTkScrollableFrame(options_frame, 
                                                  fg_color=COLORS['white'],
                                                  corner_radius=0, 
                                                  scrollbar_button_color=COLORS['neutral 3'],
                                                  scrollbar_button_hover_color=COLORS['neutral 3'])
        self.options_subframe.columnconfigure(0, weight=1, uniform=1)

        ctk.CTkLabel(self.options_subframe,
                     text='Normalization Options',
                     text_color= COLORS['primary 1'],
                     font = ('Poppins Medium', 20)).grid(row = 0, column = 0, sticky = 'w', pady = (0, 40))

        ctk.CTkLabel(self.options_subframe,
                     text='Method',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 1, column = 0, sticky = 'w', pady = (0, 16))
        
        self.norm_method_optionmenu = customCTkOptionMenu(self.options_subframe,
                            variable = self.normalization_method,
                            height=56, image_padx = 16,
                            values = ['None', 'Normalize Total', 'stSME Normalization'])
        
        self.norm_method_optionmenu.grid(row = 2, column = 0, sticky = 'we', padx = (0, 8), pady = (0, 40))

        ctk.CTkLabel(self.options_subframe,
                     text='Other',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 3, column = 0, sticky = 'w', pady = (0, 16))
        
        self.log1p_button = OtherButton(self.options_subframe,
                                        text='Log1p',
                                        var=self.do_log1p)
        
        self.log1p_button.grid(row = 4, column = 0, sticky = 'we', pady = (0, 16), padx = (0, 8))

        
        self.scale_button = OtherButton(self.options_subframe,
                                        text='Scale data',
                                        var = self.do_scale)
        
        self.scale_button.grid(row = 5, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))

        self.apply_button = ctk.CTkButton(self.options_subframe,
                                            text='Apply',
                                            font=('Poppins Semibold', 18),
                                            corner_radius=16,
                                            fg_color=COLORS['primary 1'],
                                            hover_color=COLORS['primary 1 hover'],
                                            height=56,
                                            text_color= COLORS['white'],
                                            command=self.init_normalization_process)
        
        self.apply_button.grid(row = 6, column = 0, sticky = 'we', pady = (0, 64), padx = (0, 8))

        ctk.CTkLabel(self.options_subframe,
                     text='Plotting Options',
                     text_color= COLORS['primary 1'],
                     font = ('Poppins Medium', 20)).grid(row = 7, column = 0, sticky = 'w', pady = (0, 40))

        ctk.CTkLabel(self.options_subframe,
                     text='Gene',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 8, column = 0, sticky = 'w', pady = (0, 16))
        
        self.gene_name_entry = ctk.CTkEntry(self.options_subframe,
                     placeholder_text='Enter Gene...',
                     border_width=1,
                     border_color=COLORS['border'],
                     corner_radius=16,
                     height=56,
                     font = ('Poppins Medium', 18),
                     text_color= COLORS['neutral 1'],
                     fg_color=COLORS['white'],
                     placeholder_text_color=COLORS['neutral 2'])
        self.gene_name_entry.grid(row = 9, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))

        ctk.CTkLabel(self.options_subframe,
                     text='Plot Type',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 10, column = 0, sticky = 'w', pady = (0, 16))
        
        self.plot_type_optionmenu = customCTkOptionMenu(self.options_subframe,
                            variable = self.plot_type,
                            height=56, image_padx = 16,
                            values = ['Spatial Plot', 'Histogram'])
        
        self.plot_type_optionmenu.grid(row = 11, column = 0, sticky = 'we', padx = (0, 8), pady = (0, 40))

        ctk.CTkLabel(self.options_subframe,
                     text='State',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 12, column = 0, sticky = 'w', pady = (0, 16))
        
        self.state_optionmenu = customCTkOptionMenu(self.options_subframe,
                            variable = self.state,
                            height=56, image_padx = 16,
                            values = ['Side by side', 'Raw counts', 'Normalized counts'])
        
        self.state_optionmenu.grid(row = 13, column = 0, sticky = 'we', padx = (0, 8), pady = (0, 40))               
        
        #   Place Frames    #
        self.no_plot_frame.grid(row = 1, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')
        options_frame.grid(row = 1, column = 1, pady = (24, 0), sticky='NSWE')
        buttons_frame.pack(side = 'bottom', fill = 'x', padx = 18, pady = (0, 10))
        line_canvas.pack(side = 'bottom', fill = 'x', padx = 2)
        self.options_subframe.pack(side = 'top', fill = 'both', expand = True, padx = (24,12), pady = (40, 2)) 

        #   Bindings    #
        self.gene_name_entry.bind("<Return>", lambda _: self.root.focus_set())

    def check_thread(self, thread, next_step):
        if thread.is_alive():
            self.after(100, lambda: self.check_thread(thread, next_step))
        else:
            next_step()

    def init_normalization_process(self):
        #a variable for keeping track of what's already done
        self.already_done = [name for name, value in [('Normalize', self.parent.isNormalized.get()), ('Log1p', self.parent.isLog1p.get()), ('Scale', self.parent.isScaled.get())] if value]

        #a variable for keeping track of what's selected
        self.selected = [name for name, value in [('Normalize', (self.normalization_method.get() != 'None')), ('Log1p', self.do_log1p.get()), ('Scale', self.do_scale.get())] if value]

        #a variable for keeping track of what's already done out of selection
        self.done_from_selection = [name for name in self.selected if name in self.already_done]

        if len(self.done_from_selection) == len(self.selected):
            #everything is done
            self.root.lockAppWindow()
            response = self.root.build_msgmodal(msgbox_args = {'no_button':False,
                                                                'cancel_button':False,
                                                                'yes_text': 'Okay',
                                                                'title':'Error',
                                                                'icon':'fail',
                                                                'title_color':COLORS['fail'],
                                                                'message' : 'All of the selected actions have already been applied.'})
            if response:
                self.root.releaseAppWindow()
            
        else:
            self.root.lockAppWindow()

            if len(self.done_from_selection) == 0:
                #nothing is done
                self.list_to_check = self.selected
            
            else:
                #partially done
                mapping = {'Normalize': 'Normalized', 'Log1p': 'Lograthamized', 'Scale': 'Scaled'}
                whats_done = ', '.join(mapping[item] for item in self.done_from_selection if item in mapping)

                response = self.root.build_msgmodal(msgbox_args={
                    'no_button': True,
                    'cancel_button': False,
                    'yes_text': 'Yes',
                    'title': 'Warning',
                    'icon': 'warning',
                    'title_color': COLORS['primary 1'],
                    'message': f'This dataset has already been {whats_done}.\nContinue with the other steps?',
                    'auto_size': True
                })
                if response == 'No':
                    self.root.releaseAppWindow()
                    return
                
                self.list_to_check = [name for name in self.selected if name not in self.already_done]
            
            self.final_message = []

            if self.list_to_check[0] == 'Normalize':
                #starts with normalization

                if self.normalization_method.get() == 'Normalize Total':
                    self.init_normalize_total()

                elif self.normalization_method.get() == 'stSME Normalization':
                    self.init_stSME_normalization()

            elif self.list_to_check[0] == 'Log1p':
                #starts with Log1p

                self.root.build_IndProgressModal()
                self.init_log1p()

            elif self.list_to_check[0] == 'Scale':
                #starts with scale

                self.root.build_IndProgressModal()
                self.init_scale()

            else:
                self.root.releaseAppWindow()

    def init_normalize_total(self):
        self.root.build_IndProgressModal()
        secondary_thread = threading.Thread(target=NormalizeTotal, args=(self, self.parent.adata))
        secondary_thread.start()
        self.check_thread(secondary_thread, self.init_log1p)

    def init_stSME_normalization(self):
        
        if 'X_morphology' not in self.parent.adata.obsm:
            msg = (
                "stSME Normalization uses a pre-trained CNN to extract and normalize features.\nThis process is computationally intensive and may take 10-30 minutes,\npotentially reducing your computer's responsiveness.\nIt's a one-time procedure for all stlearn tasks in the application.\nDo you want to proceed with stSME Normalization?"
            )
            response = self.root.build_msgmodal(msgbox_args={
                'no_button': True,
                'cancel_button': False,
                'yes_text': 'Yes',
                'title': 'stSME Normalization Process Confirmation',
                'icon': 'info',
                'title_color': COLORS['primary 1'],
                'message': msg,
                'auto_size': True
            })
            if response == 'Yes':
                self.root.build_DetProgressModal()
                secondary_thread = threading.Thread(target=run_stSMEnormalization, args=(self, self.parent.adata))
                secondary_thread.start()
            else:
                self.root.releaseAppWindow()
                return

        else:
            self.root.build_IndProgressModal()
            secondary_thread = threading.Thread(target=run_stSMEnormalization_noFeatureExtract, args=(self, self.parent.adata))
            secondary_thread.start()                    
                
        self.check_thread(secondary_thread, self.init_log1p)        

    def init_log1p(self):
        if 'Log1p' in self.list_to_check:
            secondary_thread = threading.Thread(target=Log1p_Transform, args=(self, self.parent.adata))
            secondary_thread.start()
            self.check_thread(secondary_thread, self.init_scale)  
        else:
            self.init_scale()

    def init_scale(self):
        if 'Scale' in self.list_to_check:
            secondary_thread = threading.Thread(target=Scale_Transform, args=(self, self.parent.adata))
            secondary_thread.start()
            self.check_thread(secondary_thread, self.on_completion)  
        else:
            self.on_completion()

    def on_completion(self):
        #if the UI text isnt changed
        if self.disabled_text.cget("text") == "Gene plotting available after dataset normalization.":
            self.disabled_text.configure(text = "Enter gene and click refresh to view plot.")
        
        if self.refresh_button.cget('state') == 'disabled':
            self.refresh_button.configure(state = 'normal')

        if self.export_button.cget('state') == 'disabled':
            self.export_button.configure(state = 'normal')

        result = ', '.join(self.final_message[:-1]) + (' and ' if len(self.final_message) > 1 else '') + self.final_message[-1] if self.final_message else ''
        msg = f"Dataset successfully {result} !"
        self.root.destroy_ProgressModal(title = 'Successfully Done', msg = msg, show = 'success')

    def plot_closeup(self, title, index):
        img = self.fullsize_plots[index]

        for toplevel in self.current_plot_toplevels:
            if toplevel.name == title:
                toplevel.deiconify()
                toplevel.lift()
                toplevel.focus_force()
                return
            
        toplevel = PlotTopLevel(name = title, parent= self, geometry=self.root.geometry(), img = img)
        self.current_plot_toplevels.append(toplevel)

    def init_plot(self):
        
        def generate_args(layer):
            if self.plot_type.get() == "Spatial Plot":
                args = {
                    'color': self.gene_name_entry.get(),
                    'colorbar' : False,
                }

                if layer == 'Raw':
                    args["layer"] = 'raw_counts'

            elif self.plot_type.get() == "Histogram":
                args = {'kde' : True,
                        'bins' : 30}              
                if layer == 'Raw':
                    args["data"] = self.parent.adata[:, gene].layers['raw_counts'].toarray().flatten()
                else:
                    args["data"] = self.parent.adata[:, gene].X.toarray().flatten()

            return args
        
        # lock app window
        
        self.root.lockAppWindow()
        
        # first check if the typed gene is in the dataset

        gene = self.gene_name_entry.get()

        if gene not in self.root.mainwindow.gene_names:
            response = self.root.build_msgmodal(msgbox_args={
                'no_button': False,
                'cancel_button': False,
                'yes_text': 'Okay',
                'title': 'Error',
                'icon': 'fail',
                'title_color': COLORS['error'],
                'message': "Gene could not be found in the dataset.\nPlease make sure the name is spelled correctly and then try again.",
                'auto_size': True
            })
            if response == 'Yes':
                if hasattr(self, 'container_parent'):
                    self.container_parent.destroy()
                    self.no_plot_frame.grid(row = 1, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')

                self.root.releaseAppWindow()
                return
    
        self.listofcanvases = []
        self.listofCB_canvases = []
        plotting_args_list = [] 

        if hasattr(self, 'container_parent'):                                                            #Check and destroy if the container frame already exists
            self.container_parent.destroy()

        if self.no_plot_frame.winfo_exists():                                                            #Check and remove the no plots frame
            self.no_plot_frame.grid_forget()

        self.container_parent = ctk.CTkFrame(self, fg_color=COLORS['background'])                         #Make the parent container frame
        self.container_parent.grid(row = 1, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')

        sbs = False
        if self.state.get() == 'Side by side':
            sbs = True                                                                                     # Use this value to properly set a side by side layout
        
        self.container_parent.columnconfigure(0, weight=1)
        if sbs:
            self.container_parent.columnconfigure(1, weight=1)
        self.container_parent.rowconfigure(0, weight=1)
                                                
        #make first left container
        name = "Raw" if self.state_optionmenu.get() == "Raw counts" or sbs else "Normalized"
        spatial_container = True if self.plot_type.get() == "Spatial Plot" else False
        
        make_container(self, title=f"{gene} {name} Counts", row= 0, col = 0, index = 0, spatial_container=spatial_container)

        if sbs:
            #make second right container
            make_container(self, title=f"{gene} Normalized Counts", row= 0, col = 1, index = 1, spatial_container=spatial_container)

        self.root.update_idletasks()                                                                                            #Update the containers on screen. 
        y = self.listofcanvases[0].winfo_height()                                                                               #Get a container height. 
        x = self.listofcanvases[0].winfo_width()                                                                                #Get a container width. 

        plotting_args_list.append(generate_args(layer=name))
        if sbs:
            plotting_args_list.append(generate_args(layer="Normalized"))

        self.root.build_IndProgressModal()
        
        if self.plot_type.get() == "Spatial Plot":
            secondary_thread = threading.Thread(target=make_squidpy_plots, args=[self, x, y, plotting_args_list, self.parent.imgsize, 'normalization', 'bottom', self.cb_canvas_height])
            secondary_thread.start()
        
    def finalize_plotting(self, plots, fullsize_plots, colorbars):
        self.plot_images = []
        self.cb_images = []
        self.fullsize_plots = fullsize_plots

        for i in range(len(self.listofcanvases)):

            # For colorbar

            cb_image = ImageTk.PhotoImage(colorbars[i])
            self.cb_images.append(cb_image)
            cb_canvas = self.listofCB_canvases[i]

            cb_canvas_width = cb_canvas.winfo_width()
            cb_canvas_height = cb_canvas.winfo_height()

            cb_x = (cb_canvas_width - cb_image.width()) // 2
            cb_y = (cb_canvas_height - cb_image.height()) // 2   

            cb_canvas.create_image(cb_x, cb_y, image=self.cb_images[i], anchor='nw')

            # For the main plot image 

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


        

