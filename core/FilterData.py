#       Import Packages     #
import customtkinter as ctk
from core.assets.themecolors import COLORS
from core.customwidgets import *
from core.validatecommands import *
from core.ctk_rangeslider import *
from core.backend.filterdata_backend import *
from core.backend.viewdata_backend import *
import threading

#       FilterData Class     #
class FilterData(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent, fg_color=COLORS['background'])
        self.parent = parent
        self.root = self.parent.parent
        self.vcmdInt = (self.register(callback))



        self.columnconfigure(1, weight=1, uniform=1)
        self.columnconfigure(0, weight=2, uniform=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=0)
        self.rowconfigure(2, weight=1)


        #   Init Variables   #
        self.n_genes_by_counts_min = ctk.IntVar(value=0)
        self.n_genes_by_counts_max = ctk.IntVar(value=20)
        self.total_counts_min = ctk.IntVar(value=0)
        self.total_counts_max = ctk.IntVar(value=20)
        self.mit_counts_min = ctk.IntVar(value=0)
        self.mit_counts_max = ctk.IntVar(value=20)
        self.metric = ctk.StringVar(value='Expressed Genes')
        self.plot_type = ctk.StringVar(value='Spatial Plot')
        self.current_plot_toplevels = []

        #   Build Main Components   #
        self.initialize()

    def initialize(self):
        #   Functions   #

        def make_stat(parent, var, name, padx=(0,0)):
            frame = ctk.CTkFrame(parent, fg_color='transparent')
            count_label = ctk.CTkLabel(frame, textvariable = var, font=('Poppins SemiBold', 24), text_color=COLORS['primary 1'])
            name_label = ctk.CTkLabel(frame, text=name, font=('Poppins Medium', 14), text_color=COLORS['neutral 2'])
            count_label.pack(side='top', pady = (0, 8))
            name_label.pack(side='bottom')
            frame.pack(side='left', padx = padx)

        def auto_filter_command(*args):
            pass

        def update_min_spots_entry(var, *args):
            try:
                current_val = int(self.min_spots_entry.get())
            except:                                                 #Entry is empty
                current_val = 0
            if var == 'inc':
                new_val = current_val + 1
            elif var == 'dec':
                new_val = current_val - 1
                if new_val < 0:
                    new_val = str(0)
            self.min_spots_entry.delete(0, tk.END)
            self.min_spots_entry.insert(0, str(new_val))

        #   FOR EXPRESSED GENES     #
        def update_eg_entry(*args):
            self.max_exp_genes_entry.delete(0, tk.END)
            self.max_exp_genes_entry.insert(0, self.n_genes_by_counts_max.get())
            self.min_exp_genes_entry.delete(0, tk.END)
            self.min_exp_genes_entry.insert(0, self.n_genes_by_counts_min.get())                  
            
        def eg_min_focusout(*args):                                                                                                             #A function that is called when user enters the min value for filter by gene expression.
            if self.min_exp_genes_entry.get() == '' or int(self.min_exp_genes_entry.get()) > self.n_genes_by_counts_max.get():                  #Check if the entry is more than the max value
                self.n_genes_by_counts_min.set(0)                                                                                               #Set the value n_genes_by_counts_min to 0
                self.min_exp_genes_entry.delete(0, tk.END)                                                                                      #Delete the entry
                self.min_exp_genes_entry.insert(0, 0)                                                                                           #Make it 0
            else:
                self.n_genes_by_counts_min.set(int(self.min_exp_genes_entry.get()))                                                             #Otherwise, just set the entry to n_genes_by_counts_min
            self.focus()                                                                                                                        #Focus out of the text box

        def eg_max_focusout(*args):                                                                                                             #A function that is called when user enters the max value for filter by gene expression.
            if self.max_exp_genes_entry.get() == '' or int(self.max_exp_genes_entry.get()) <= self.n_genes_by_counts_min.get() or int(self.max_exp_genes_entry.get()) > self.parent.ngenesbycounts.get():                               #Check if it is less than min or more than n_genes_by_counts
                self.n_genes_by_counts_max.set(self.parent.ngenesbycounts.get())
                self.max_exp_genes_entry.delete(0, tk.END)
                self.max_exp_genes_entry.insert(0, self.parent.ngenesbycounts.get())
            else:
                self.n_genes_by_counts_max.set(int(self.max_exp_genes_entry.get()))
            self.focus()

        #   FOR TOTAL COUNTS    #
        def update_tc_entry(*args):
            self.max_total_counts_entry.delete(0, tk.END)
            self.max_total_counts_entry.insert(0, self.total_counts_max.get())
            self.min_total_counts_entry.delete(0, tk.END)
            self.min_total_counts_entry.insert(0, self.total_counts_min.get()) 

        def tc_min_focusout(*args):
            if self.min_total_counts_entry.get() == '' or int(self.min_total_counts_entry.get()) > self.total_counts_max.get():
                self.total_counts_min.set(0)
                self.min_total_counts_entry.delete(0, tk.END) 
                self.min_total_counts_entry.insert(0, 0)
            else:
                self.total_counts_min.set(int(self.min_total_counts_entry.get()))
            self.focus()

        def tc_max_focusout(*args):
            if self.max_total_counts_entry.get() == '' or int(self.max_total_counts_entry.get()) <= self.total_counts_min.get() or int(self.max_total_counts_entry.get()) > self.parent.total_counts.get():
                self.total_counts_max.set(self.parent.total_counts.get())
                self.max_total_counts_entry.delete(0, tk.END)
                self.max_total_counts_entry.insert(0, self.parent.total_counts.get())
            else:
                self.total_counts_max.set(int(self.max_total_counts_entry.get()))
            self.focus()  

        #   FOR MIT COUNTS  #   
        def trace_mc(*args):
            if self.parent.mit_counts.get() == 0:
                self.slider_mit_counts.configure(state = 'disabled')
                self.min_mit_counts_entry.configure(state = 'disabled')
                self.max_mit_counts_entry.configure(state = 'disabled')
                self.mit_counts_max.set(value=0)
                self.mit_counts_min.set(value=0)
            else:
                self.slider_mit_counts.configure(state = 'normal')
                self.min_mit_counts_entry.configure(state = 'normal')
                self.max_mit_counts_entry.configure(state = 'normal')
                self.slider_mit_counts.configure(to=self.parent.mit_counts.get())
                self.mit_counts_max.set(self.parent.mit_counts.get())
            update_mc_entry()

        def update_mc_entry(*args):
            self.max_mit_counts_entry.delete(0, tk.END)
            self.max_mit_counts_entry.insert(0, self.mit_counts_max.get())
            self.min_mit_counts_entry.delete(0, tk.END)
            self.min_mit_counts_entry.insert(0, self.mit_counts_min.get()) 

        def mc_min_focusout(*args):
            if self.min_mit_counts_entry.get() == '' or int(self.min_mit_counts_entry.get()) > self.mit_counts_max.get():
                self.mit_counts_min.set(0)
                self.min_mit_counts_entry.delete(0, tk.END) 
                self.min_mit_counts_entry.insert(0, 0)
            else:
                self.mit_counts_min.set(int(self.min_mit_counts_entry.get()))
            self.focus()

        def mc_max_focusout(*args):
            if self.max_mit_counts_entry.get() == '' or int(self.max_mit_counts_entry.get()) < self.mit_counts_min.get() or int(self.max_mit_counts_entry.get()) > self.parent.mit_counts.get():
                self.mit_counts_max.set(self.parent.mit_counts.get())
                self.max_mit_counts_entry.delete(0, tk.END)
                self.max_mit_counts_entry.insert(0, self.parent.mit_counts.get())
            else:
                self.mit_counts_max.set(int(self.max_mit_counts_entry.get()))
            self.focus() 



        #   Top Bar  #
        topbar_frame = ctk.CTkFrame(self, fg_color='transparent')
        customTkLabel(topbar_frame, text='Preprocessing & Filtering', font=('Poppins Medium', 40), height=54, textMoveY= -20, text_color=COLORS['primary 1'], bg_color=COLORS['background']).pack(anchor = 'w', side = 'top', fill='x')
        ctk.CTkLabel(topbar_frame, textvariable = self.parent.dataset_displayname, font = ('Poppins Medium', 18), text_color=COLORS['neutral 2']).pack(anchor = 'w', side = 'bottom', pady = 6)
        topbar_frame.grid(row = 0, column = 0, columnspan = 2, sticky = 'NSWE')

        #   Options     #
        options_frame = ctk.CTkFrame(self, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])

        #   Create the line #
        line_canvas = tk.Canvas(options_frame, height = 3, bg='white', highlightthickness=0)
        line_canvas.create_line((0,0,999,1), dash=(4,8), fill = COLORS['neutral 2'])

        #   Create the buttons for plotting #
        buttons_frame = ctk.CTkFrame(options_frame, fg_color='white', corner_radius=0, height=60)
        self.export_button = ctk.CTkButton(buttons_frame, 
                            text = 'Export Plot', 
                            font=('Poppins Semibold', 18),
                            fg_color=COLORS['neutral 4'],
                            hover_color=COLORS['neutral 4'],
                            corner_radius=16,
                            state='disabled')
        self.export_button.pack(pady = 12, side = 'left', fill = 'both', expand = True, padx = (0, 24))
        self.refresh_button = ctk.CTkButton(buttons_frame,
                                            text='Refresh Plot',
                                            font=('Poppins Semibold', 18),
                                            corner_radius=16,
                                            fg_color=COLORS['primary 1'],
                                            hover_color=COLORS['primary 1 hover'],
                                            height=64,
                                            text_color= COLORS['white'],
                                            command=self.init_plot)
        self.refresh_button.pack(pady = 12, side = 'right', fill = 'both', expand = True)

        #   Filter Stats    #
        self.filter_stats_frame = ctk.CTkFrame(self, fg_color='transparent')
        filter_stats_col1 = ctk.CTkFrame(self.filter_stats_frame, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])
        filter_stats1_mainframe = ctk.CTkFrame(filter_stats_col1, fg_color='transparent')
        filter_stats_col2 = ctk.CTkFrame(self.filter_stats_frame, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])
        filter_stats2_mainframe = ctk.CTkFrame(filter_stats_col2, fg_color='transparent')

        make_stat(filter_stats1_mainframe, var=self.parent.total_genes, name = 'TOTAL GENES', padx=(36,0))
        make_stat(filter_stats1_mainframe, var=self.parent.default_total_genes, name = 'STARTING GENES', padx=(36,36))
        make_stat(filter_stats1_mainframe, var=self.parent.removed_total_genes, name = 'REMOVED GENES', padx=(0,36))

        make_stat(filter_stats2_mainframe, var=self.parent.total_spots, name = 'TOTAL SPOTS', padx=(36,0))
        make_stat(filter_stats2_mainframe, var=self.parent.default_total_spots, name = 'STARTING SPOTS', padx=(36,36))
        make_stat(filter_stats2_mainframe, var=self.parent.removed_total_spots, name = 'REMOVED SPOTS', padx=(0,36))

        self.filter_stats_frame.columnconfigure(0, weight=1, uniform=1)
        self.filter_stats_frame.columnconfigure(1, weight=1, uniform=1)
        filter_stats1_mainframe.pack(anchor = 'center', pady = 40)
        filter_stats2_mainframe.pack(anchor = 'center', pady = 40)
        filter_stats_col1.grid(row = 0, column = 0, padx = (0, 24), sticky='WE')
        filter_stats_col2.grid(row = 0, column = 1, sticky='WE')

        #   No Plot    #
        self.no_plot_frame = ctk.CTkFrame(self, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])
        norm_icon = ctk.CTkImage(light_image = Image.open('app\\core\\assets\\filterdata\\filter_data.png'),size=(48,48))
        self.no_plot_frame_centered = ctk.CTkFrame(self.no_plot_frame, fg_color='transparent')
        ctk.CTkLabel(self.no_plot_frame_centered, text="", image=norm_icon).pack(pady = (0, 16))
        self.disabled_text = ctk.CTkLabel(self.no_plot_frame_centered, text="Refresh to view selected plot.", font = ('Poppins Medium', 16), text_color=COLORS['neutral 2'])
        self.disabled_text.pack()
        self.no_plot_frame_centered.pack(expand = True)

        #   Options Subframe    #
        self.options_subframe = ctk.CTkScrollableFrame(options_frame, 
                                                  fg_color=COLORS['white'],
                                                  corner_radius=0, 
                                                  scrollbar_button_color=COLORS['neutral 3'],
                                                  scrollbar_button_hover_color=COLORS['neutral 3'])
        self.options_subframe.columnconfigure(0, weight=1, uniform=1)


        #   Subsection 1    #
        options_subsection_1 = ctk.CTkFrame(self.options_subframe, fg_color='transparent')

        ctk.CTkLabel(options_subsection_1,
                     text='Filtering Options',
                     text_color= COLORS['primary 1'],
                     font = ('Poppins Medium', 20)).pack(side = 'left')
        
        self.auto_filter_button = ctk.CTkButton(options_subsection_1,
                                            text='Auto Filter',
                                            font=('Poppins Semibold', 16),
                                            corner_radius=16,
                                            fg_color=COLORS['white'],
                                            hover_color=COLORS['background'],
                                            height=36,
                                            text_color= COLORS['primary 1'],
                                            border_width= 2, 
                                            border_color=COLORS['primary 1'],
                                            command=auto_filter_command)
        
        self.auto_filter_button.pack(side = 'right', padx = (0, 12))
        
        options_subsection_1.grid(row = 0, column = 0, sticky = 'we', pady = (0, 40))

        #   Subsection 2    #
        #   Filter Genes by Spots
        
        ctk.CTkLabel(self.options_subframe,
                     text='Minimum spots',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 1, column = 0, sticky = 'w', pady = (0, 16))
        
        options_subsection_2 = ctk.CTkFrame(self.options_subframe, fg_color='transparent')
        options_subsection_2.columnconfigure(0, weight=1, uniform=1)
        options_subsection_2.columnconfigure(1, weight=4, uniform=1)
        options_subsection_2.columnconfigure(2, weight=1, uniform=1)

        ctk.CTkButton(options_subsection_2,
            text= '-',
            text_color=COLORS['primary 1'],
            height=56,
            width=56,
            fg_color=COLORS['primary 2'],
            corner_radius=16,
            hover_color=COLORS['primary 2 hover'],
            font=('Poppins', 20),
            command=lambda: update_min_spots_entry('dec')
            ).grid(row = 0, column = 0, sticky = 'w', padx = (0, 16))
        
        self.min_spots_entry = ctk.CTkEntry(options_subsection_2,
                        placeholder_text='30',
                        border_width=1,
                        border_color=COLORS['border'],
                        corner_radius=16,
                        height=56,
                        font = ('Poppins Medium', 18),
                        text_color= COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        placeholder_text_color=COLORS['neutral 2'],
                        justify='center',
                        validatecommand=(self.vcmdInt, '%P'),
                        validate='all')
        
        self.min_spots_entry.grid(row = 0, column = 1, sticky = 'we')

        ctk.CTkButton(options_subsection_2,
            text= '+',
            text_color=COLORS['primary 1'],
            height=56,
            width=56,
            fg_color=COLORS['primary 2'],
            corner_radius=16,
            hover_color=COLORS['primary 2 hover'],
            font=('Poppins', 20),
            command=lambda: update_min_spots_entry('inc')
            ).grid(row = 0, column = 2, sticky = 'e', padx = (16, 0))

        options_subsection_2.grid(row = 2, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))

        #   Subsection 3    #
        #   Filter by Expressed Genes (n_genes_by_counts)

        ctk.CTkLabel(self.options_subframe,
                     text='Filter by Expressed Genes',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 3, column = 0, sticky = 'w', pady = (0, 16))
        
        options_subsection_3 = ctk.CTkFrame(self.options_subframe, fg_color='transparent')
        options_subsection_3.columnconfigure(0, weight=3, uniform=1)
        options_subsection_3.columnconfigure(1, weight=1, uniform=1)
        options_subsection_3.columnconfigure(2, weight=3, uniform=1)

        self.slider_exp_genes = CTkRangeSlider(options_subsection_3, 
                                               variables=[self.n_genes_by_counts_min, self.n_genes_by_counts_max], 
                                               from_=0, 
                                               to=10,
                                               height = 16,
                                               fg_color	 = COLORS['neutral 4'],
                                               progress_color = COLORS['primary 1'],
                                               button_color= COLORS['primary 1'],
                                               button_hover_color= COLORS['primary 1 hover'],
                                               command=update_eg_entry)
        
        self.slider_exp_genes.grid(row = 0, column = 0, columnspan = 3, pady = (0, 16), sticky = 'we')

        self.min_exp_genes_entry = ctk.CTkEntry(options_subsection_3,
                        border_width=1,
                        border_color=COLORS['border'],
                        corner_radius=16,
                        height=56,
                        font = ('Poppins Medium', 18),
                        text_color= COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        justify='center',
                        validatecommand=(self.vcmdInt, '%P'),
                        validate='all')
        
        self.min_exp_genes_entry.grid(row = 1, column = 0, sticky = 'we')

        ctk.CTkLabel(options_subsection_3,
                     text='to',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 18)).grid(row = 1, column = 1, sticky = 'we', padx = 16)

        self.max_exp_genes_entry = ctk.CTkEntry(options_subsection_3,
                        border_width=1,
                        border_color=COLORS['border'],
                        corner_radius=16,
                        height=56,
                        font = ('Poppins Medium', 18),
                        text_color= COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        justify='center',
                        validatecommand=(self.vcmdInt, '%P'),
                        validate='all')
        
        self.max_exp_genes_entry.grid(row = 1, column = 2, sticky = 'we')             
        
        options_subsection_3.grid(row = 4, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))

        #   Subsection 4    #
        #   Filter by Total Counts (total_counts)

        ctk.CTkLabel(self.options_subframe,
                     text='Filter by Total Counts',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 5, column = 0, sticky = 'w', pady = (0, 16))
        
        options_subsection_4 = ctk.CTkFrame(self.options_subframe, fg_color='transparent')
        options_subsection_4.columnconfigure(0, weight=3, uniform=1)
        options_subsection_4.columnconfigure(1, weight=1, uniform=1)
        options_subsection_4.columnconfigure(2, weight=3, uniform=1)

        self.slider_total_counts = CTkRangeSlider(options_subsection_4, 
                                               variables=[self.total_counts_min, self.total_counts_max], 
                                               from_=0, 
                                               to=10,
                                               height = 16,
                                               fg_color	 = COLORS['neutral 4'],
                                               progress_color = COLORS['primary 1'],
                                               button_color= COLORS['primary 1'],
                                               button_hover_color= COLORS['primary 1 hover'],
                                               command=update_tc_entry)
        
        self.slider_total_counts.grid(row = 0, column = 0, columnspan = 3, pady = (0, 16), sticky = 'we')

        self.min_total_counts_entry = ctk.CTkEntry(options_subsection_4,
                        border_width=1,
                        border_color=COLORS['border'],
                        corner_radius=16,
                        height=56,
                        font = ('Poppins Medium', 18),
                        text_color= COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        justify='center',
                        validatecommand=(self.vcmdInt, '%P'),
                        validate='all')
        
        self.min_total_counts_entry.grid(row = 1, column = 0, sticky = 'we')

        ctk.CTkLabel(options_subsection_4,
                     text='to',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 18)).grid(row = 1, column = 1, sticky = 'we', padx = 16)

        self.max_total_counts_entry = ctk.CTkEntry(options_subsection_4,
                        border_width=1,
                        border_color=COLORS['border'],
                        corner_radius=16,
                        height=56,
                        font = ('Poppins Medium', 18),
                        text_color= COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        justify='center',
                        validatecommand=(self.vcmdInt, '%P'),
                        validate='all')
        
        self.max_total_counts_entry.grid(row = 1, column = 2, sticky = 'we')  

        options_subsection_4.grid(row = 6, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))

        #   Subsection 5    #
        #   Filter by Percent Mitochondrial Counts (pct_counts_mt)

        ctk.CTkLabel(self.options_subframe,
                     text='Filter by % Mitochondrial Counts',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 16)).grid(row = 7, column = 0, sticky = 'w', pady = (0, 16))
        
        options_subsection_5 = ctk.CTkFrame(self.options_subframe, fg_color='transparent')
        options_subsection_5.columnconfigure(0, weight=3, uniform=1)
        options_subsection_5.columnconfigure(1, weight=1, uniform=1)
        options_subsection_5.columnconfigure(2, weight=3, uniform=1)

        self.slider_mit_counts = CTkRangeSlider(options_subsection_5, 
                                               variables=[self.mit_counts_min, self.mit_counts_max], 
                                               from_=0, 
                                               to=10,
                                               height = 16,
                                               fg_color	 = COLORS['neutral 4'],
                                               progress_color = COLORS['primary 1'],
                                               button_color= COLORS['primary 1'],
                                               button_hover_color= COLORS['primary 1 hover'],
                                               command=update_mc_entry)
        
        self.slider_mit_counts.grid(row = 0, column = 0, columnspan = 3, pady = (0, 16), sticky = 'we')

        self.min_mit_counts_entry = ctk.CTkEntry(options_subsection_5,
                        border_width=1,
                        border_color=COLORS['border'],
                        corner_radius=16,
                        height=56,
                        font = ('Poppins Medium', 18),
                        text_color= COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        justify='center',
                        validatecommand=(self.vcmdInt, '%P'),
                        validate='all')
        
        self.min_mit_counts_entry.grid(row = 1, column = 0, sticky = 'we')

        ctk.CTkLabel(options_subsection_5,
                     text='to',
                     text_color= COLORS['neutral 2'],
                     font = ('Poppins Medium', 18)).grid(row = 1, column = 1, sticky = 'we', padx = 16)

        self.max_mit_counts_entry = ctk.CTkEntry(options_subsection_5,
                        border_width=1,
                        border_color=COLORS['border'],
                        corner_radius=16,
                        height=56,
                        font = ('Poppins Medium', 18),
                        text_color= COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        justify='center',
                        validatecommand=(self.vcmdInt, '%P'),
                        validate='all')
        
        self.max_mit_counts_entry.grid(row = 1, column = 2, sticky = 'we')  

        options_subsection_5.grid(row = 8, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))

        ctk.CTkButton(self.options_subframe,
                        text='Apply Filters',
                        font=('Poppins Semibold', 18),
                        corner_radius=16,
                        fg_color=COLORS['primary 1'],
                        hover_color=COLORS['primary 1 hover'],
                        height=64,
                        text_color= COLORS['white'],
                        command=self.apply_filter).grid(row = 9, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))

        #   Subsection 6    #
        #   Plotting Options    #
        options_subsection_6 = ctk.CTkFrame(self.options_subframe, fg_color='transparent')
        options_subsection_6.columnconfigure(0, weight=1, uniform=1)

        ctk.CTkLabel(options_subsection_6,
                    text='Plotting Options',
                    text_color= COLORS['primary 1'],
                    font = ('Poppins Medium', 20)).grid(row = 0, column = 0, sticky = 'w', pady= (0, 40))
        
        ctk.CTkLabel(options_subsection_6,
                    text='Metric',
                    text_color= COLORS['neutral 2'],
                    font = ('Poppins Medium', 16)).grid(row = 1, column = 0, sticky = 'w', pady = (0, 16))
        
        self.metric_optionmenu = customCTkOptionMenu(options_subsection_6,
                            variable = self.metric,
                            height=56, image_padx = 16,
                            values = ['Expressed Genes', 'Total Counts', '% Mitochondrial Counts'])
        
        self.metric_optionmenu.grid(row = 2, column = 0, sticky = 'we', pady = (0, 40))

        ctk.CTkLabel(options_subsection_6,
                    text='Plot type',
                    text_color= COLORS['neutral 2'],
                    font = ('Poppins Medium', 16)).grid(row = 3, column = 0, sticky = 'w', pady = (0, 16))        
        
        self.plot_type_optionmenu = customCTkOptionMenu(options_subsection_6,
                            variable = self.plot_type,
                            height=56, image_padx = 16,
                            values = ['Spatial Plot', 'Violin Plot', 'Bar Plot'])
        
        self.plot_type_optionmenu.grid(row = 4, column = 0, sticky = 'we', pady = (0, 40))

        options_subsection_6.grid(row = 10, column = 0, sticky = 'we', pady = (0, 40), padx = (0, 8))

        #   Place Frames    #
        self.no_plot_frame.grid(row = 2, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')
        options_frame.grid(row = 1, column = 1, pady = (24, 0), sticky='NSWE', rowspan = 2)
        buttons_frame.pack(side = 'bottom', fill = 'x', padx = 18, pady = (0, 10))
        line_canvas.pack(side = 'bottom', fill = 'x', padx = 2)
        self.options_subframe.pack(side = 'top', fill = 'both', expand = True, padx = (24,12), pady = (40, 2)) 
        self.filter_stats_frame.grid(row = 1, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')

        #   Bindings    #
        self.min_spots_entry.bind("<Return>", lambda _: self.root.focus_set())

        self.parent.ngenesbycounts.trace("w", lambda *args: (self.slider_exp_genes.configure(to=self.parent.ngenesbycounts.get()), self.n_genes_by_counts_max.set(self.parent.ngenesbycounts.get()), update_eg_entry()))    #Whenever the data is imported, set the slider values.
        self.parent.total_counts.trace("w", lambda *args: (self.slider_total_counts.configure(to=self.parent.total_counts.get()), self.total_counts_max.set(self.parent.total_counts.get()), update_tc_entry()))
        self.parent.mit_counts.trace("w", trace_mc)

        self.max_exp_genes_entry.bind('<Return>', eg_max_focusout)
        self.min_exp_genes_entry.bind('<Return>', eg_min_focusout)
        self.max_total_counts_entry.bind('<Return>', tc_max_focusout)
        self.min_total_counts_entry.bind('<Return>', tc_min_focusout)
        self.max_mit_counts_entry.bind('<Return>', mc_max_focusout)
        self.min_mit_counts_entry.bind('<Return>', mc_min_focusout)
        
    def init_plot(self):

        color = {
                'Expressed Genes' : "n_genes_by_counts",
                'Total Counts' : "total_counts",
                '% Mitochondrial Counts' : "pct_counts_mt"
            }
        
        limits = {
                'Expressed Genes' : (self.max_exp_genes_entry.get(),self.min_exp_genes_entry.get()),
                'Total Counts' : (self.max_total_counts_entry.get(),self.min_total_counts_entry.get()),
                '% Mitochondrial Counts' : (self.max_mit_counts_entry.get(),self.min_mit_counts_entry.get())
            }
        
        def generate_args():

            if self.plot_type.get() == "Spatial Plot":
                upper_lim_vmax, lower_lim_vmin = limits[self.metric.get()]
                args = {
                    'color': color[self.metric.get()],
                    'colorbar' : False,
                    'vmin' : lower_lim_vmin,
                    'vmax' : upper_lim_vmax,
                }

            elif self.plot_type.get() == 'Bar Plot':
                args = {
                    'data' : self.parent.adata.obs[color[self.metric.get()]],
                    'kde' : False
                }
                
            elif self.plot_type.get() == 'Violin Plot':
                args = {
                    'color' : color[self.metric.get()],
                    'jitter' : 0.4
                }

            return args

        # lock app window
        
        self.root.lockAppWindow()

        if self.metric == '% Mitochondrial Counts' and self.parent.mit_counts.get() == 0:
            response = self.root.build_msgmodal(msgbox_args={
                'no_button': False,
                'cancel_button': False,
                'yes_text': 'Okay',
                'title': 'Error',
                'icon': 'fail',
                'title_color': COLORS['error'],
                'message': "Zero Mitochondrial Counts in Dataset.",
                'auto_size': True
            })

            self.root.releaseAppWindow()
            return

        self.listofcanvases = []
        plotting_args_list = [] 
        if self.plot_type.get() == 'Spatial Plot':
            self.listofCB_canvases = []

        if hasattr(self, 'container_parent'):                                                            #Check and destroy if the container frame already exists
            self.container_parent.destroy()

        if self.no_plot_frame.winfo_exists():                                                            #Check and remove the no plots frame
            self.no_plot_frame.grid_forget()

        self.container_parent = ctk.CTkFrame(self, fg_color=COLORS['background'])                         #Make the parent container frame
        self.container_parent.grid(row = 2, column = 0, padx = (0, 24), pady = (24, 0), sticky='NSWE')
        self.container_parent.columnconfigure(0 , weight=1, uniform=1)
        self.container_parent.rowconfigure(0 , weight=1, uniform=1)

        spatial_container = True if self.plot_type.get() == "Spatial Plot" else False
        
        make_container(self, title=self.metric.get(), row= 0, col = 0, index = 0, spatial_container=spatial_container)

        self.root.update_idletasks()                                                                                            #Update the containers on screen. 
        y = self.listofcanvases[0].winfo_height()                                                                               #Get a container height. 
        x = self.listofcanvases[0].winfo_width()                                                                                #Get a container width. 

        plotting_args_list.append(generate_args())

        self.root.build_IndProgressModal()

        if self.plot_type.get() == "Spatial Plot":
            secondary_thread = threading.Thread(target=make_squidpy_plots, args=[self, x, y, plotting_args_list, self.parent.imgsize, 'filterdata', 'bottom', self.cb_canvas_height])
            secondary_thread.start()

        if self.plot_type.get() == 'Violin Plot':                                              # For Bar plots
            data = self.parent.adata.obs[color[self.metric.get()]]
            upper_lim , lower_lim = limits[self.metric.get()]
            upper_lim , lower_lim = int(upper_lim) , int(lower_lim)

            secondary_thread = threading.Thread(target=make_violin_plots, args=[self, data, x, y, upper_lim, lower_lim, self.metric.get()])
            secondary_thread.start()   
        
        if self.plot_type.get() == 'Bar Plot':                                              # For Bar plots
            data = self.parent.adata.obs[color[self.metric.get()]]
            upper_lim , lower_lim = limits[self.metric.get()]
            upper_lim , lower_lim = int(upper_lim) , int(lower_lim)

            secondary_thread = threading.Thread(target=make_bar_plots, args=[self, data, x, y, upper_lim, lower_lim, self.metric.get(), 'Spots'])
            secondary_thread.start()   

    def finalize_plotting(self, plots, fullsize_plots, colorbars = None):
        self.plot_images = []
        self.fullsize_plots = fullsize_plots

        if self.plot_type.get() == 'Spatial Plot':

            if colorbars is None:
               raise ValueError('No colorbars provided.') 
            
            self.cb_images = []
            # For colorbar

            cb_image = ImageTk.PhotoImage(colorbars[0])
            self.cb_images.append(cb_image)
            cb_canvas = self.listofCB_canvases[0]

            cb_canvas_width = cb_canvas.winfo_width()
            cb_canvas_height = cb_canvas.winfo_height()

            cb_x = (cb_canvas_width - cb_image.width()) // 2
            cb_y = (cb_canvas_height - cb_image.height()) // 2   

            cb_canvas.create_image(cb_x, cb_y, image=self.cb_images[0], anchor='nw')

        # For the main plot image 

        plot_image = ImageTk.PhotoImage(plots[0])           # Convert the first image to a PhotoImage
        self.plot_images.append(plot_image)                 # Keep a reference to the PhotoImage
        canvas = self.listofcanvases[0]

        # Get canvas width and height
        canvas_width = canvas.winfo_width()
        canvas_height = canvas.winfo_height()

        # Calculate coordinates to place the image at the center of the canvas
        x = (canvas_width - plot_image.width()) // 2
        y = (canvas_height - plot_image.height()) // 2

        canvas.create_image(x, y, image=self.plot_images[0], anchor='nw')

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
            
    def apply_filter(self):
        self.root.lockAppWindow()
        min_spots = self.min_spots_entry.get() if self.min_spots_entry.get() is not '' else 'None'
        confirmation_string = f"Proceeding will apply the following filters:\n\n1. Minimum Spots per Gene: {min_spots}\n2. Expressed Genes: {self.min_exp_genes_entry.get()} to {self.max_exp_genes_entry.get()}\n3. Total Counts: {self.min_total_counts_entry.get()} to {self.max_total_counts_entry.get()}\n4. % Mitochondrial Counts: {self.min_mit_counts_entry.get()} to {self.max_mit_counts_entry.get()}\n\nDo you want to proceed?"

        response = self.root.build_msgmodal(msgbox_args={
                'no_button': False,
                'cancel_button': True,
                'yes_text': 'Continue',
                'title': 'Confirmation',
                'icon': 'info',
                'title_color': COLORS['primary 1'],
                'message': confirmation_string,
                'auto_size': True
            })

        if response == 'Cancel':
            self.root.releaseAppWindow()
            return

        adata_copy = self.parent.adata.copy()

        if min_spots != 'None':
            sc.pp.filter_genes(adata_copy, min_cells=int(min_spots))
        adata_copy = adata_copy[(self.n_genes_by_counts_min.get() <= adata_copy.obs.n_genes_by_counts) & (adata_copy.obs.n_genes_by_counts <= self.n_genes_by_counts_max.get())]
        adata_copy = adata_copy[(self.total_counts_min.get() <= adata_copy.obs.total_counts) & (adata_copy.obs.total_counts <= self.total_counts_max.get())]
        adata_copy = adata_copy[(self.mit_counts_min.get() <= adata_copy.obs.pct_counts_mt) & (adata_copy.obs.pct_counts_mt <= self.mit_counts_max.get())]

        if len(adata_copy.obs) == 0:
            error = "The applied filters have resulted in the exclusion of all spots from the dataset.\nTo ensure meaningful analysis, please adjust the filtering parameters to a less stringent level.\nThe dataset will be reverted to its previous state."
            response = self.root.build_msgmodal(msgbox_args={
                    'no_button': False,
                    'cancel_button': False,
                    'yes_text': 'Okay',
                    'title': 'Error',
                    'icon': 'error',
                    'title_color': COLORS['error'],
                    'message': error,
                    'auto_size': True
                })

            if response:
                self.root.releaseAppWindow()           
                del adata_copy
                return
        
        if len(adata_copy.var) == 0:
            error = "The applied filters have resulted in the exclusion of all genes from the dataset.\nTo ensure meaningful analysis, please adjust the filtering parameters to a less stringent level.\nThe dataset will be reverted to its previous state."
            response = self.root.build_msgmodal(msgbox_args={
                    'no_button': False,
                    'cancel_button': False,
                    'yes_text': 'Okay',
                    'title': 'Error',
                    'icon': 'error',
                    'title_color': COLORS['error'],
                    'message': error,
                    'auto_size': True
                })

            if response:
                self.root.releaseAppWindow()           
                del adata_copy
                return
        
        self.parent.adata = adata_copy.copy()
        del adata_copy

        ngenesbycounts_max = self.parent.adata.obs["n_genes_by_counts"].max()
        totalcounts_max = self.parent.adata.obs["total_counts"].max()
        pctcountsmt_max = self.parent.adata.obs["pct_counts_mt"].max()

        self.parent.ngenesbycounts.set(ngenesbycounts_max)
        self.parent.total_counts.set(totalcounts_max)
        self.parent.mit_counts.set(0 if math.isnan(pctcountsmt_max) else pctcountsmt_max)

        self.parent.gene_names = self.parent.adata.var.index.tolist()
        self.parent.total_spots.set(len(self.parent.adata.obs))
        self.parent.total_genes.set(len(self.parent.adata.var))

        response = self.root.build_msgmodal(msgbox_args={
                'no_button': False,
                'cancel_button': False,
                'yes_text': 'Okay',
                'title': 'Successful',
                'icon': 'info',
                'title_color': COLORS['primary 1'],
                'message': 'Selected filters successfully applied!',
                'auto_size': True
            })
        
        self.root.releaseAppWindow()