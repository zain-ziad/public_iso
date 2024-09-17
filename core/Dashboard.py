#       Import Packages     #
import customtkinter as ctk
from core.assets.themecolors import COLORS

#       Dashboard Class     #
class Dashboard(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent, fg_color=COLORS['background'])
        ctk.CTkLabel(self, text='Dashboard!').pack()