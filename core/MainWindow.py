#       Import Packages     #
import customtkinter as ctk
from tkinter import ttk
from core.assets.themecolors import COLORS
from core.Dashboard import Dashboard
from core.ManageData import ManageData
from core.ViewData import ViewData
from core.FilterData import FilterData
from core.Normalization import Normalization


#       MainWindow Class     #

class MainWindow(ttk.Notebook):
    def __init__(self, parent, menu, command_queue):

        #Configure style of notebook
        style=ttk.Style()
        style.layout('TNotebook', [])
        style.layout('TNotebook.Tab', [])
        style.configure('TNotebook', background=COLORS['background'])

        super().__init__(parent)
        self.menu = menu
        self.parent = parent
        self.command_queue = command_queue
        self.pack(side = 'left', fill = 'both', expand = True)

        #   Variables   #
        self.adata = None
        self.bdata = None
        self.gene_names = []
        self.imgsize = None
        self.dataname = None
        self.dataset_displayname = ctk.StringVar(value='No dataset imported')
        self.total_spots = ctk.StringVar(value='0')
        self.total_genes = ctk.StringVar(value='0')
        self.default_total_spots = ctk.StringVar(value='0')
        self.default_total_genes = ctk.StringVar(value='0')
        self.removed_total_genes = ctk.StringVar(value='0')
        self.removed_total_spots = ctk.StringVar(value='0')
        self.ngenesbycounts = ctk.IntVar(value=10)
        self.total_counts = ctk.IntVar(value=10)
        self.mit_counts = ctk.IntVar(value=10)
        self.isNormalized = ctk.BooleanVar(value=False)
        self.isScaled = ctk.BooleanVar(value=False)
        self.isLog1p = ctk.BooleanVar(value=False)
        self.isIntegrated = False
        self.isMerged = False
        self.isClustered = None
        self.isVariable = None
        self.isXenium = None
        self.xenium_imgpath = None

        #   Bindings
        def update_removed_spots(*args):
            new_spots = int(self.default_total_spots.get()) - int(self.total_spots.get())
            self.removed_total_spots.set(str(new_spots))

        def update_removed_genes(*args):
            new_genes = int(self.default_total_genes.get()) - int(self.total_genes.get())
            self.removed_total_genes.set(str(new_genes))

        self.total_spots.trace_add("write", update_removed_spots)
        self.total_genes.trace_add("write", update_removed_genes)

        #Initialize the Main Window Notebook
        self.initialize()

    def initialize(self):
        
        #Define Tabs
        self.Dashboard = Dashboard(self)
        self.ManageData = ManageData(self)
        self.ViewData = ViewData(self)
        self.FilterData = FilterData(self)
        self.Normalization = Normalization(self)
        self.add(self.Dashboard)
        self.add(self.ManageData, padding=32, sticky='nsew')
        self.add(self.ViewData, padding=32, sticky='nsew')
        self.add(self.FilterData, padding=32, sticky='nsew')
        self.add(self.Normalization, padding=32, sticky='nsew')