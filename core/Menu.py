#       Import Packages     #
import customtkinter as ctk
from PIL import Image, ImageTk
from core.assets.themecolors import COLORS


#       Menu Class     #
class Menu(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent, corner_radius=0, width=349, fg_color=COLORS['white'])
        self.parent = parent
        self.pack(side = 'left', fill = 'y')
        self.current_selection = 'dashboard'            #A value to track the current selected button

        #Load icons to use in menu buttons
        self.button_icons = {
            'dashboard' : self.load_menubutton_icon('dashboard'),
            'manage_data' : self.load_menubutton_icon('manage_data'),
            'view_data' : self.load_menubutton_icon('view_data'),
            'filter_data' : self.load_menubutton_icon('filter_data'),
            'normalization' : self.load_menubutton_icon('normalization'),
            'integration' : self.load_menubutton_icon('integration'),
            'variable_genes' : self.load_menubutton_icon('variable_genes'),
            'clustering' : self.load_menubutton_icon('clustering'),
            'diff_exp' : self.load_menubutton_icon('diff_exp'),
            'pathway' : self.load_menubutton_icon('pathway'),
            'decon' : self.load_menubutton_icon('decon'),
            'trajectory' : self.load_menubutton_icon('trajectory')
        }

        #These are the arguments that are used to return a menu button to its "selected" state. 
        self.button_selected_args = {
            'font' : ('Poppins SemiBold', 18), 
            'text_color' : 'white', 
            'fg_color' : COLORS['primary 1'], 
            'hover' : False
        }

        #Initialize the Side Menu
        self.initialize()

    #A function to load the icons for the menu buttons
    def load_menubutton_icon(self, name):
        ico_b = Image.open(f'app\\core\\assets\\menu_icons\\{name}.png')
        ico_b = ico_b.resize((24, 24), Image.Resampling.LANCZOS)
        ico_b = ImageTk.PhotoImage(ico_b)

        ico_w = Image.open(f'app\\core\\assets\\menu_icons\\white\\{name}.png')
        ico_w = ico_w.resize((24, 24), Image.Resampling.LANCZOS)
        ico_w = ImageTk.PhotoImage(ico_w)

        return [ico_b, ico_w]

    #A function to return a menu button to the "clear" state.
    def clear_button_selection(self):
        current_but = getattr(self, f'{self.current_selection}_button')
        current_but.configure(font = ('Poppins Medium', 18), text_color = COLORS['primary 1'], fg_color = COLORS['white'], hover = True, image =  self.button_icons[self.current_selection][0])

    #Functions when button is pressed
    def call_dashboard(self):
        if self.current_selection != 'dashboard':
            self.clear_button_selection()
        self.dashboard_button.configure(image =  self.button_icons['dashboard'][1], **self.button_selected_args)
        self.current_selection = 'dashboard'
        self.parent.mainwindow.select(0)

    def call_manage_data(self):
        if self.current_selection != 'manage_data':
            self.clear_button_selection()
        self.manage_data_button.configure(image =  self.button_icons['manage_data'][1], **self.button_selected_args)
        self.current_selection = 'manage_data'
        self.parent.mainwindow.select(1)

    def call_view_data(self):
        if self.current_selection != 'view_data':
            self.clear_button_selection()
        self.view_data_button.configure(image =  self.button_icons['view_data'][1], **self.button_selected_args)
        self.current_selection = 'view_data'
        self.parent.mainwindow.select(2)

    def call_filter_data(self):
        if self.current_selection != 'filter_data':
            self.clear_button_selection()
        self.filter_data_button.configure(image =  self.button_icons['filter_data'][1], **self.button_selected_args)
        self.current_selection = 'filter_data'
        self.parent.mainwindow.select(3)

    def call_normalization(self):
        if self.current_selection != 'normalization':
            self.clear_button_selection()
        self.normalization_button.configure(image =  self.button_icons['normalization'][1], **self.button_selected_args)
        self.current_selection = 'normalization'
        self.parent.mainwindow.select(4)

    def call_integration(self):
        if self.current_selection != 'integration':
            self.clear_button_selection()
        self.integration_button.configure(image =  self.button_icons['integration'][1], **self.button_selected_args)
        self.current_selection = 'integration'
        #self.parent.mainwindow.select(1)

    def call_variable_genes(self):
        if self.current_selection != 'variable_genes':
            self.clear_button_selection()
        self.variable_genes_button.configure(image =  self.button_icons['variable_genes'][1], **self.button_selected_args)
        self.current_selection = 'variable_genes'
        #self.parent.mainwindow.select(1)

    def call_clustering(self):
        if self.current_selection != 'clustering':
            self.clear_button_selection()
        self.clustering_button.configure(image =  self.button_icons['clustering'][1], **self.button_selected_args)
        self.current_selection = 'clustering'
        #self.parent.mainwindow.select(1)

    def call_diff_exp(self):
        if self.current_selection != 'diff_exp':
            self.clear_button_selection()
        self.diff_exp_button.configure(image =  self.button_icons['diff_exp'][1], **self.button_selected_args)
        self.current_selection = 'diff_exp'
        #self.parent.mainwindow.select(1)    

    def call_pathway(self):
        if self.current_selection != 'pathway':
            self.clear_button_selection()
        self.pathway_button.configure(image =  self.button_icons['pathway'][1], **self.button_selected_args)
        self.current_selection = 'pathway'
        #self.parent.mainwindow.select(1)

    def call_decon(self):
        if self.current_selection != 'decon':
            self.clear_button_selection()
        self.decon_button.configure(image =  self.button_icons['decon'][1], **self.button_selected_args)
        self.current_selection = 'decon'
        #self.parent.mainwindow.select(1)

    #Initialize Function
    def initialize(self):

        #Pack app logo
        logo = Image.open('app\\core\\assets\\app_logo.png').resize((260, 80), Image.Resampling.LANCZOS)
        ctk.CTkLabel(self, image=ImageTk.PhotoImage(logo), text='').pack(pady = 44.5, padx = 44.5)

        #Build the main frame
        self.main_frame = ctk.CTkFrame(self, fg_color='transparent')

        #These are the arguments that are used to make a menu button in its "clear" state. 
        self.button_clear_args = {
            'master' : self.main_frame,
            'text_color' : COLORS['primary 1'],
            'font' : ('Poppins Medium', 18),
            'fg_color' : COLORS['white'],
            'hover_color' : COLORS['background'],
            'anchor' : 'w', 
            'width' : 301, 
            'height' : 56, 
            'corner_radius' : 12
        }

        #Build the buttons
        self.dashboard_button = ctk.CTkButton(text='Dashboard', image= self.button_icons['dashboard'][0], command=self.call_dashboard, **self.button_clear_args)
        self.manage_data_button = ctk.CTkButton(text='Manage Data', image= self.button_icons['manage_data'][0], command=self.call_manage_data, **self.button_clear_args)
        self.view_data_button = ctk.CTkButton(text='View Data', image= self.button_icons['view_data'][0], command=self.call_view_data, **self.button_clear_args)
        self.filter_data_button = ctk.CTkButton(text='Preprocessing & Filtering', image= self.button_icons['filter_data'][0], command=self.call_filter_data, **self.button_clear_args)
        self.normalization_button = ctk.CTkButton(text='Normalization', image= self.button_icons['normalization'][0], command=self.call_normalization, **self.button_clear_args)
        self.integration_button = ctk.CTkButton(text='Spatial Integration', image= self.button_icons['integration'][0], command=self.call_integration, **self.button_clear_args)
        self.variable_genes_button = ctk.CTkButton(text='Spatially Variable Genes', image= self.button_icons['variable_genes'][0], command=self.call_variable_genes, **self.button_clear_args)
        self.clustering_button = ctk.CTkButton(text='Clustering', image= self.button_icons['clustering'][0], command=self.call_clustering, **self.button_clear_args)
        self.diff_exp_button = ctk.CTkButton(text='Differential Expression', image= self.button_icons['diff_exp'][0], command=self.call_diff_exp, **self.button_clear_args)
        self.pathway_button = ctk.CTkButton(text='Pathway Analysis', image= self.button_icons['pathway'][0], command=self.call_pathway, **self.button_clear_args)
        self.decon_button = ctk.CTkButton(text='Spatial Deconvolution', image= self.button_icons['decon'][0], command=self.call_decon, **self.button_clear_args)

        #Pack Widgets
        self.main_frame.pack()
        self.dashboard_button.pack()
        self.manage_data_button.pack(pady = 16)
        self.view_data_button.pack()
        self.filter_data_button.pack(pady = 16)
        self.normalization_button.pack()
        self.integration_button.pack(pady = 16)
        self.variable_genes_button.pack()
        self.clustering_button.pack(pady = 16)
        self.diff_exp_button.pack()
        self.pathway_button.pack(pady = 16)
        self.decon_button.pack()