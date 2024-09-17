#       Import Packages     #
import customtkinter as ctk
from queue import Queue
import time
from PIL import ImageGrab, ImageTk, Image, ImageEnhance
import tkinter as tk
from core.Menu import Menu
from core.MainWindow import MainWindow
from core.assets.font import font_init
from core.assets.themecolors import COLORS
from core.customwidgets import msgbox
import warnings

#       App     #
class App(ctk.CTk):
    def __init__ (self, size):

        #main setup
        super().__init__()
        self.title('Interface for Spatial Omics')           #Name of the window                       
        self.resizable(False, False)                          #Resizable Window
        self.size = size                                    #Size of the window

        #queue variables
        self.command_queue = Queue()                        #Used for multiprocessing

        #Init Font
        font_init(self)

        #core widgets
        self.menu = Menu(self)
        self.mainwindow = MainWindow(self, self.menu, self.command_queue)

        #run commands
        self.center_window(self.size[0], self.size[1])      #Center the window

        #scaling widgets according to resolution
        self.scaling_factor = self.size[0]/1920
        ctk.set_widget_scaling(self.scaling_factor)

        # Suppress all UserWarnings
        warnings.filterwarnings("ignore", category=UserWarning)

        #self.attributes("-fullscreen", True)
        self.run()

    #methods
    def check_queue(self):
        while not self.command_queue.empty():
            command, args = self.command_queue.get()
            if command == 'destroy_progress_window':
                self.destroy_ProgressModal(**args)
            if command == 'managedata_refreshstats':
                self.mainwindow.ManageData.refresh_statistics()
            if command == 'viewdata_finalize_plotting':
                self.mainwindow.ViewData.finalize_plotting(*args)
            if command == 'normalization_finalize_plotting':
                self.mainwindow.Normalization.finalize_plotting(*args)
            if command == 'filterdata_finalize_plotting':
                self.mainwindow.FilterData.finalize_plotting(*args)
            if command == 'normalization_on_completion':
                self.mainwindow.Normalization.on_completion(*args)

        self.after(100, self.check_queue)
    
    def lockAppWindow(self):
        # Capture a screenshot of the application window
        x = self.winfo_rootx()
        y = self.winfo_rooty()
        w = x + self.winfo_width()
        h = y + self.winfo_height()
        screenshot = ImageGrab.grab().crop((x, y, w, h)).convert('RGBA')
        self.screenshot_size = screenshot.size
        overlayed_screenshot = Image.alpha_composite(screenshot, Image.new("RGBA", screenshot.size, (0, 0, 0, int(255 * 0.6))))

        # Convert to a format suitable for Tkinter
        self.screenshot_image = ImageTk.PhotoImage(overlayed_screenshot)

        # Create a full-size label to display the image
        self.overlay_label = tk.Label(self, image=self.screenshot_image)
        self.overlay_label.place(x=0, y=0, relwidth=1, relheight=1)

        # Disable main window
        self.attributes("-disabled", True)

    def releaseAppWindow(self):
        # To bypass the black artifacts, first we make a new window with the screenshot label. Destroy the original label, and then destroy the window.
    
        self.overlay_window = tk.Toplevel(self)
        self.overlay_window.geometry(f"{self.winfo_width()}x{self.winfo_height()}+{self.winfo_rootx()}+{self.winfo_rooty()}")
        self.overlay_window.overrideredirect(True)
        label = tk.Label(self.overlay_window, image=self.screenshot_image)
        label.place(x=0, y=0, relwidth=1, relheight=1)
        self.overlay_label.destroy()
        self.after(200, self.overlay_window.destroy)
        self.attributes("-disabled", False)
        self.attributes('-topmost', True)
        self.attributes('-topmost', False)

    def checkAppWindow(self):
        if self.cget('state') == 'disabled':
            return True
        else:
            return False

    def build_msgmodal(self, msgbox_args = {}):
        response = ctk.StringVar('')
        msg = msgbox(parent=self, callback = lambda i : response.set(i), **msgbox_args)
        self.wait_window(msg)
        self.attributes('-topmost', True)
        self.attributes('-topmost', False)
        return response.get()

    def build_IndProgressModal(self, msg = 'Loading...'):
        self.indprogressmodal = msgbox(parent=self, callback = None, title=msg, cancel_button=False, yes_button=False, no_button=False, icon='loading', message='')

    def build_DetProgressModal(self, msg = 'Please wait...'):
        self.detprogressmodal = msgbox(parent=self, callback = None, title=msg, cancel_button=False, yes_button=False, no_button=False, icon='det_progress_bar', message='det_loading')

    def destroy_ProgressModal(self, value = 'ind', title = None, msg = None, show = 'none', releaseAppWindow = True):

        def cont_func():
            if show == 'success':
                success = msgbox(parent=self, callback = None, 
                                title=title, 
                                cancel_button=False, 
                                yes_button=True, 
                                no_button=False,
                                icon='success', 
                                message=msg,
                                yes_text = 'Continue')
                self.wait_window(success)

            if show == 'fail':
                fail = msgbox(parent=self, callback = None, 
                                title=title, 
                                cancel_button=False, 
                                yes_button=True, 
                                no_button=False,
                                icon='fail', 
                                message=msg,
                                yes_text = 'Try Again',
                                title_color=COLORS['error'])
                self.wait_window(fail)

            self.attributes('-topmost', True)
            self.attributes('-topmost', False)
            
            if releaseAppWindow:
                self.releaseAppWindow()

        if value == 'ind':
            self.after(1000, self.indprogressmodal.destroy_ind_modal)
            self.after(1010, cont_func)
        elif value == 'det':
            self.detprogressmodal.destroy()
            cont_func()

    def center_window(self, width, height):
        # get screen width and height
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        # calculate position x and y coordinates
        x = (screen_width/2) - (width/2)
        y = (screen_height/2) - (height/2)
        self.geometry('%dx%d+%d+%d' % (width, height, x, y))

    def run(self):
        if __name__ == "__main__":
            self.check_queue()
            self.mainloop()
        
App(size = (1920, 1080))