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

        #   Build Main Components   #
        self.initialize()


        #   Bindings    #
        self.parent.isNormalized.trace_add('write', self.finalize_normalization)
        self.parent.isLog1p.trace_add('write', self.finalize_normalization)
        self.parent.isScaled.trace_add('write', self.finalize_normalization)

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
        self.disabled_text = ctk.CTkLabel(self.no_plot_frame_centered, text="Gene plotting will become available after you've normalized the dataset.", font = ('Poppins Medium', 16), text_color=COLORS['neutral 2'])
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
                            corner_radius=16)
        self.export_button.pack(pady = 12, side = 'left', fill = 'y', expand = True)
        self.refresh_button = ctk.CTkButton(buttons_frame,
                                            text='Refresh',
                                            font=('Poppins Semibold', 18),
                                            corner_radius=16,
                                            fg_color=COLORS['primary 1'],
                                            hover_color=COLORS['primary 1 hover'],
                                            height=64,
                                            text_color= COLORS['white'])
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
                     text='Normalization Options',
                     text_color= COLORS['primary 1'],
                     font = ('Poppins Medium', 20)).grid(row = 0, column = 0, columnspan = 2, sticky = 'w', pady = (0, 40))

        ctk.CTkLabel(self.options_subframe,
                     text='Method',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 1, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
        
        self.norm_method_optionmenu = customCTkOptionMenu(self.options_subframe,
                            variable = self.normalization_method,
                            height=56, image_padx = 16,
                            values = ['None', 'Normalize Total', 'stSME Normalization'])
        
        self.norm_method_optionmenu.grid(row = 2, column = 0, columnspan = 2, sticky = 'we', padx = (0, 8), pady = (0, 40))

        ctk.CTkLabel(self.options_subframe,
                     text='Other',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 3, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
        
        self.log1p_button = OtherButton(self.options_subframe,
                                        text='Log1p',
                                        var=self.do_log1p)
        
        self.log1p_button.grid(row = 4, column = 0, columnspan = 2, sticky = 'we', pady = (0, 16), padx = (0, 8))

        
        self.scale_button = OtherButton(self.options_subframe,
                                        text='Scale data',
                                        var = self.do_scale)
        
        self.scale_button.grid(row = 5, column = 0, columnspan = 2, sticky = 'we', pady = (0, 40), padx = (0, 8))

        self.apply_button = ctk.CTkButton(self.options_subframe,
                                            text='Apply',
                                            font=('Poppins Semibold', 18),
                                            corner_radius=16,
                                            fg_color=COLORS['primary 1'],
                                            hover_color=COLORS['primary 1 hover'],
                                            height=56,
                                            text_color= COLORS['white'],
                                            command=self.initiate_normalization)
        
        self.apply_button.grid(row = 6, column = 0, columnspan = 2, sticky = 'we', pady = (0, 64), padx = (0, 8))

        ctk.CTkLabel(self.options_subframe,
                     text='Plotting Options',
                     text_color= COLORS['primary 1'],
                     font = ('Poppins Medium', 20)).grid(row = 7, column = 0, columnspan = 2, sticky = 'w', pady = (0, 40))

        ctk.CTkLabel(self.options_subframe,
                     text='Gene',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 8, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
        
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
        self.gene_name_entry.grid(row = 9, column = 0, columnspan = 2, sticky = 'we', pady = (0, 40), padx = (0, 8))

        ctk.CTkLabel(self.options_subframe,
                     text='Plot Type',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 10, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
        
        self.plot_type_optionmenu = customCTkOptionMenu(self.options_subframe,
                            variable = self.plot_type,
                            height=56, image_padx = 16,
                            values = ['Spatial Plot', 'Histogram'])
        
        self.plot_type_optionmenu.grid(row = 11, column = 0, columnspan = 2, sticky = 'we', padx = (0, 8), pady = (0, 40))

        ctk.CTkLabel(self.options_subframe,
                     text='State',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 12, column = 0, columnspan = 2, sticky = 'w', pady = (0, 16))
        
        self.state_optionmenu = customCTkOptionMenu(self.options_subframe,
                            variable = self.state,
                            height=56, image_padx = 16,
                            values = ['Side by side', 'Raw Counts', 'Normalized Counts'])
        
        self.state_optionmenu.grid(row = 13, column = 0, columnspan = 2, sticky = 'we', padx = (0, 8), pady = (0, 40))               
        
        #   Place Frames    #
        self.no_plot_frame.grid(row = 1, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')
        options_frame.grid(row = 1, column = 1, pady = (24, 0), sticky='NSWE')
        buttons_frame.pack(side = 'bottom', fill = 'x', padx = 18, pady = (0, 10))
        line_canvas.pack(side = 'bottom', fill = 'x', padx = 2)
        self.options_subframe.pack(side = 'top', fill = 'both', expand = True, padx = (24,12), pady = (40, 2)) 

        #   Bindings    #
        self.gene_name_entry.bind("<Return>", lambda _: self.root.focus_set())

    def initiate_normalization(self):
        self.root.lockAppWindow()
        adata = self.parent.adata

        if self.normalization_method.get() == 'None' and not self.do_log1p.get() and not self.do_scale.get():
            response = self.root.build_msgmodal(msgbox_args = {'no_button':False,
                                                                'cancel_button':False,
                                                                'yes_text': 'Okay',
                                                                'title':'Error',
                                                                'icon':'fail',
                                                                'title_color':COLORS['error'],
                                                                'message' : 'No options selected. Please select an option and try again.',
                                                                'auto_size': True})
            if response == 'Yes':
                self.root.releaseAppWindow()
                return
            
        actions = []
        self.selected_actions = 0

        # Check each condition
        if self.normalization_method.get() != 'None':
            self.selected_actions += 1
            if self.parent.isNormalized.get():
                actions.append('Normalized')

        if self.do_log1p.get():
            self.selected_actions += 1
            if self.parent.isLog1p.get():
                actions.append('Log1p transformed')

        if self.do_scale.get():
            self.selected_actions += 1
            if self.parent.isScaled.get():
                actions.append('Scaled')

        # Check if all selected actions are already done
        if len(actions) == self.selected_actions and self.selected_actions > 0:
            response = self.root.build_msgmodal(msgbox_args={
                'no_button': False,
                'cancel_button': False,
                'yes_text': 'Okay',
                'title': 'Error',
                'icon': 'fail',
                'title_color': COLORS['error'],
                'message': 'All selected options have already been applied to the dataset.',
                'auto_size': True
            })

            if response:
                self.root.releaseAppWindow()
                return

        if actions:
            whats_done = ", ".join(actions)
            response = self.root.build_msgmodal(msgbox_args={
                'no_button': False,
                'cancel_button': True,
                'yes_text': 'Continue',
                'title': 'Warning',
                'icon': 'warning',
                'title_color': COLORS['warning'],
                'message': f'This dataset has already been {whats_done}.\nContinue with the other steps?',
                'auto_size': True
            })

            if response == 'Cancel':
                self.root.releaseAppWindow()
                return

        def check_thread(thread, next_step):
            print(thread)
            if thread is None:
                next_step()
            else:
                if thread.is_alive():
                    self.after(100, lambda: check_thread(thread, next_step))
                else:
                    next_step()

        def start_normalization():
            secondary_thread = None

            if self.normalization_method.get() != 'None':
                self.selected_actions-=1
                if not self.parent.isNormalized.get():
                    if self.selected_actions > 0:
                        action_left = True
                    else:
                        action_left = False
                    if self.normalization_method.get() == 'Normalize Total':
                        self.root.build_IndProgressModal()
                        secondary_thread = threading.Thread(target=NormalizeTotal, args=(self, adata, action_left,))                                                                #Start stSME + extract_feature
                        secondary_thread.start()

                    if self.normalization_method.get() == 'stSME Normalization':                                                                                                #stLearn based stSME Normalization
                        if 'X_morphology' not in adata.obsm:
                            msg = "stSME Normalization uses a pre-trained CNN to extract and normalize features.\nThis process is computationally intensive and may take 10-30 minutes,\npotentially reducing your computer's responsiveness.\nIt's a one-time procedure for all stlearn tasks in the application.\nDo you want to proceed with stSME Normalization?"
                            response = self.root.build_msgmodal(msgbox_args = {'no_button':True,
                                                                        'cancel_button':False,
                                                                        'yes_text': 'Yes',
                                                                        'title':'stSME Normalization Process Confirmation',
                                                                        'icon':'info',
                                                                        'title_color':COLORS['primary 1'],
                                                                        'message' : msg,
                                                                        'auto_size': True})
                            if response == 'Yes':
                                self.root.build_DetProgressModal()
                                secondary_thread = threading.Thread(target=run_stSMEnormalization, args=(self, adata, action_left,))                                                                #Start stSME + extract_feature
                                secondary_thread.start()
                            else:
                                self.root.releaseAppWindow()
                                return
                        else:  
                            self.root.build_IndProgressModal()
                            secondary_thread = threading.Thread(target=run_stSMEnormalization_noFeatureExtract, args=(self, adata, action_left,))                                                                #Start stSME + extract_feature
                            secondary_thread.start()

            check_thread(secondary_thread, start_log1p)

        def start_log1p():
            secondary_thread = None
            print('started log1p')

            if self.do_log1p.get():
                self.selected_actions-=1

                if not self.parent.isLog1p.get():
                    if self.selected_actions > 0:
                        action_left = True
                    else:
                        action_left = False

                    if not hasattr(self.root, 'indprogressmodal') or not self.root.indprogressmodal.winfo_exists():
                        self.root.build_IndProgressModal()

                    print(self.selected_actions, action_left)

                    secondary_thread = threading.Thread(target=Log1p_Transform, args=(self, adata, action_left,))
                    secondary_thread.start()

            check_thread(secondary_thread, start_scaling)

        def start_scaling():
            secondary_thread = None
            print('starting scaling')

            if self.do_scale.get():
                self.selected_actions-=1

                if not self.parent.isScaled.get(): 
                    if self.selected_actions > 0:
                        action_left = True
                    else:
                        action_left = False    

                    if not hasattr(self.root, 'indprogressmodal') or not self.root.indprogressmodal.winfo_exists():
                        self.root.build_IndProgressModal()

                    print(self.selected_actions, action_left)

                    secondary_thread = threading.Thread(target=Scale_Transform, args=(self, adata, action_left,))
                    secondary_thread.start()

        start_normalization()

    def finalize_normalization(self, *args):
        if not self.plot_enabled:
            self.disabled_text.configure(text = 'Refresh a gene plot to visualize here.')
            self.plot_enabled = True
        
            