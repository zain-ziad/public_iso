import tkinter as tk
from typing import Tuple
import customtkinter as ctk
from core.assets.themecolors import COLORS
from PIL import Image, ImageTk
import sys
import _thread
import time
from math import floor

class customTkLabel(ctk.CTkCanvas):
    #text : The Text to display.
    #height : Height of the text box in pixels.
    #textMoveY : A variable that moves the text in the textbox.
    #font : Font of the text. 
    #text_color : Color of the text.
    #bg_color : Text background color. 
    def __init__(self, parent, text="", height = 100, textMoveY = -15, font=("Calibri", 140), text_color = 'black', bg_color = 'red', **kwargs):
        super().__init__(parent, height=height, bg=bg_color, highlightthickness=0, **kwargs)
        self.grid()
        self.txtid = self.create_text(0, textMoveY, text=text, fill=text_color, font=font, anchor='nw')

class customCTkOptionMenu(ctk.CTkOptionMenu):
    def __init__(self, parent, height = 64,
                        font=('Poppins Medium', 18),
                        text_color=COLORS['neutral 1'],
                        fg_color=COLORS['white'],
                        button_color=COLORS['white'],
                        button_hover_color=COLORS['white'],
                        button_fg_color = COLORS['white'], 
                        corner_radius= 16,
                        dynamic_resizing=False,
                        dropdown_fg_color='#f0f0f0',
                        dropdown_font= ('Poppins Regular', 14),
                        dropdown_hover_color= COLORS['neutral 3'],
                        dropdown_text_color= COLORS['neutral 1'],
                        add_border = True,
                        text = None,
                        image_padx = 20, 
                        **kwargs):
        
        super().__init__(parent,
                        height = height,
                        font=font,
                        text_color=text_color,
                        fg_color=fg_color,
                        button_color=button_color,
                        button_hover_color=button_hover_color,
                        corner_radius= corner_radius,
                        dynamic_resizing=dynamic_resizing,
                        dropdown_fg_color=dropdown_fg_color,
                        dropdown_font= dropdown_font,
                        dropdown_hover_color= dropdown_hover_color,
                        dropdown_text_color= dropdown_text_color,
                        **kwargs)
        #self.pack()
        self.label_text = text
        ico = Image.open(f'app\\core\\assets\\arrow_drop_down.png')
        self.dropdown_image = self.dropdown_image = ctk.CTkImage(light_image=ico, size=(25,15))
        self.image_label = ctk.CTkLabel(self, text="", image=self.dropdown_image, fg_color=button_fg_color)
        self._canvas.delete("dropdown_arrow")

        # Getting information to place the image
        grid_info = self._text_label.grid_info()
        grid_info["padx"], grid_info["sticky"] = image_padx, "e"
        self.image_label.grid(**grid_info)

        # If text exists:
        if self.label_text:
            self._text_label.configure(text=self.label_text)

        # Adding events for clicks and hovers
        self.image_label.bind("<Button-1>", self._clicked)
        self.image_label.bind("<Enter>", self._on_enter)
        self.image_label.bind("<Leave>", self._on_leave)

        #Add border
        if add_border:
            self._canvas.configure(highlightthickness=1, highlightcolor=COLORS['border'])

    def _on_enter(self, event=0):
        super()._on_enter(event)
        
        if self.image_label:
            color = self._apply_appearance_mode(self._button_hover_color)
            self.image_label.configure(fg_color=color, bg_color=color)

    def _on_leave(self, event=0):
        super()._on_leave(event)

        if self.image_label:
            color = self._apply_appearance_mode(self._button_color)
            self.image_label.configure(fg_color=color, bg_color=color)

    def _dropdown_callback(self, value: str):
        self._current_value = value
        if not self.label_text:
            self._text_label.configure(text=self._current_value)

        if self._variable is not None:
            self._variable_callback_blocked = True
            self._variable.set(self._current_value)
            self._variable_callback_blocked = False

        if self._command is not None:
            self._command(self._current_value)

class msgbox(ctk.CTkToplevel):
    def __init__(self,
                 parent,
                 callback,
                 width: int = 600,
                 height: int = 300,
                 title: str = "Message Box",
                 message: str = "This is a CTkMessagebox!",
                 border_color: str = COLORS['white'],
                 fg_color: str = COLORS['white'],
                 icon: str = 'upload',
                 corner_radius: int = 15,
                 cancel_button: bool = True,
                 no_button: bool = True,
                 yes_button: bool = True,
                 yes_text: str = 'Yes',
                 title_color: str = COLORS['primary 1'],
                 auto_size: bool = True
                 ):
        
        super().__init__()
        self.parent = parent
        self.callback = callback

        #animation vars
        self.animation_steps = 30
        self.current_step = 0

        self.title(title)
        self.resizable(width=False, height=False)
        self.overrideredirect(1)
        self.attributes("-topmost", True)
        self.attributes("-toolwindow", True)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
        if sys.platform.startswith("win"):
            self.transparent_color = self._apply_appearance_mode(self.cget("fg_color"))
            self.attributes("-transparentcolor", self.transparent_color)
        elif sys.platform.startswith("darwin"):
            self.transparent_color = 'systemTransparent'
            self.attributes("-transparent", True)
        else:
            self.transparent_color = '#000001'
            corner_radius = 0

        self.config(background=self.transparent_color)
        self.protocol("WM_DELETE_WINDOW", None)
        self.round_corners = corner_radius

        self.frame_top = ctk.CTkFrame(self, corner_radius=self.round_corners, width=width, border_width=0,
                                            bg_color=self.transparent_color, fg_color=COLORS['white'])
        self.frame_top.grid(sticky="nswe")

        self.frame_top.bind("<B1-Motion>", self.move_window)
        self.frame_top.bind("<ButtonPress-1>", self.oldxyset)

        #   CONTENT     #
        self.midframe = ctk.CTkFrame(self.frame_top, fg_color=COLORS['white'])

        if icon == 'loading':
            self.gif_label = ctk.CTkLabel(self.midframe, text='', corner_radius=0)
            self.gif = GifPlayer(self.gif_label, "app\\core\\assets\\spinner.gif", 0.025)
            self.gif_label.grid(row=0, column=0, columnspan=3)
            self.gif.play()

        elif icon == 'det_progress_bar':
            self.progress_bar = ctk.CTkProgressBar(self.midframe, 
                                                   mode='determinate', 
                                                   determinate_speed=0.5, 
                                                   width=280, 
                                                   height=12, 
                                                   progress_color=COLORS['primary 1'],
                                                   border_color=COLORS['border'], 
                                                   fg_color=COLORS['background'])
            self.progress_bar.set(0)
            self.progress_bar.grid(row=0, column=0, columnspan=3, padx = 20)
            
        elif icon != 'no icon':
            width = int(80*self.parent.scaling_factor)
            height = int(80*self.parent.scaling_factor)

            if icon == 'upload':
                icon_image = ImageTk.PhotoImage(Image.open('app\\core\\assets\\upload_msgbox.png').resize((width, height), Image.Resampling.LANCZOS))
            if icon == 'success':
                icon_image = ImageTk.PhotoImage(Image.open('app\\core\\assets\\success_msgbox.png').resize((width, height), Image.Resampling.LANCZOS))
            if icon == 'fail' or icon == 'error':
                icon_image = ImageTk.PhotoImage(Image.open('app\\core\\assets\\failed_msgbox.png').resize((width, height), Image.Resampling.LANCZOS))
            if icon == 'info':
                icon_image = ImageTk.PhotoImage(Image.open('app\\core\\assets\\info_msgbox.png').resize((width, height), Image.Resampling.LANCZOS))
            if icon == 'warning':
                icon_image = ImageTk.PhotoImage(Image.open('app\\core\\assets\\warning_msgbox.png').resize((width, height), Image.Resampling.LANCZOS))

            ctk.CTkLabel(self.midframe, text="", image=icon_image).grid(row=0, column=0, columnspan=3)

        ctk.CTkLabel(self.midframe, text=title, font = ('Poppins Semibold', 24), text_color=title_color).grid(row = 1, column = 0, columnspan = 3, pady = (28, 0))
        if message != "":
            if message == 'det_loading':
                self.text_box = ctk.CTkLabel(self.midframe, font = ('Poppins Medium', 20), text_color=COLORS['neutral 2'], text="Initiating...")
            else:
                self.text_box = ctk.CTkLabel(self.midframe, font = ('Poppins', 16), text_color=COLORS['neutral 1'], text= message)
            self.text_box.grid(row = 2, column = 0, columnspan = 3, pady = (0, 32), sticky = 'nswe')

        # Calculate the number of buttons
        button_count = sum([cancel_button, no_button, yes_button])

        # Determine column positions based on button count
        positions = {
            1: [1],
            2: [0, 2],
            3: [0, 1, 2]
        }

        current_col = 0

        # Cancel Button
        if cancel_button:
            ctk.CTkButton(self.midframe, 
                        fg_color=COLORS['white'], 
                        width=128, 
                        height=48, 
                        corner_radius=16, 
                        font=('Poppins Semibold', 18), 
                        text_color=COLORS['neutral 2'], 
                        text='Cancel', 
                        border_color=COLORS['neutral 2'],
                        hover_color=COLORS['background'], 
                        border_width=2,
                        command=lambda: self.button_event(choice='Cancel')).grid(row=3, column=positions[button_count][current_col], sticky='we', padx=12)
            current_col += 1

        # No Button
        if no_button:
            ctk.CTkButton(self.midframe, 
                        fg_color=COLORS['white'], 
                        width=128, 
                        height=48, 
                        corner_radius=16, 
                        font=('Poppins Semibold', 18), 
                        text_color=COLORS['neutral 2'], 
                        text='No', 
                        border_color=COLORS['neutral 2'], 
                        hover_color=COLORS['background'], 
                        border_width=2,
                        command=lambda: self.button_event(choice='No')).grid(row=3, column=positions[button_count][current_col], sticky='we', padx=12)
            current_col += 1

        # Yes Button
        if yes_button:
            ctk.CTkButton(self.midframe, 
                        fg_color=COLORS['primary 1'], 
                        width=128, 
                        height=48, 
                        corner_radius=16, 
                        font=('Poppins Semibold', 18), 
                        text_color=COLORS['white'], 
                        hover_color=COLORS['primary 1 hover'],
                        text=yes_text,
                        command=lambda: self.button_event(choice='Yes')).grid(row=3, column=positions[button_count][current_col], sticky='we', padx=12)
            
        # Define the desired pad values for the mainframe
        desired_padx = 42
        desired_pady = 42

        self.midframe.pack(pady=desired_pady, padx=desired_padx, anchor= 'center')

        if auto_size:
            self.update_idletasks()
            
            if self.parent.scaling_factor <= 1:
                width = int(self.winfo_reqwidth() + (self.winfo_reqwidth() * (1-self.parent.scaling_factor)))
                height = int(self.winfo_reqheight() + (self.winfo_reqheight() * (1-self.parent.scaling_factor)))
            else:
                width = self.winfo_reqwidth()
                height = self.winfo_reqheight()                

            # Ensure width is not less than the minimum
            if width < 400 * self.parent.scaling_factor:
                width = int(400 * self.parent.scaling_factor)

            # Calculate spawn positions
            self.spawn_x = int(self.parent.winfo_width() * 0.5 + self.parent.winfo_x() - 0.5 * width + 7)
            self.spawn_y = int(self.parent.winfo_height() * 0.5 + self.parent.winfo_y() - 0.5 * height + 20)
            
            # Define the starting size
            self.starting_width = width * 0.95
            self.starting_height = height * 0.95

            # Define the transition animation
            def animate():
                self.current_step += 1

                # Calculate the current size and position
                step_width = self.starting_width + (width - self.starting_width) * (self.current_step / self.animation_steps)
                step_height = self.starting_height + (height - self.starting_height) * (self.current_step / self.animation_steps)
                current_x = self.spawn_x + (width - step_width) / 2
                current_y = self.spawn_y + (height - step_height) / 2

                # Update the window geometry
                self.geometry(f"{int(step_width)}x{int(step_height)}+{int(current_x)}+{int(current_y)}")

                # Continue animation or finalize
                if self.current_step < self.animation_steps:
                    self.after(1, animate)
                else:
                    self.geometry(f"{width}x{height}+{self.spawn_x}+{self.spawn_y}")
                    self.midframe.pack_configure(padx=desired_padx, pady=desired_pady)

            # Set geometry
            self.midframe.pack_forget()
            self.geometry(f"{self.starting_width}x{self.starting_height}+{self.spawn_x + width // 2}+{self.spawn_y + height // 2}")
            self.after(1, animate)

    def button_event(self, choice):
        if self.callback:
            self.callback(choice)
        self.destroy()

    def update_det_progress(self, mins, percent, elapsed, *args):
        self.progress_bar.step()
        self.text_box.configure(text = f"{percent}% complete. {mins} minutes left.\nElapsed time: {elapsed} mins.")

    def oldxyset(self, event):
        self.oldx = event.x
        self.oldy = event.y
    
    def move_window(self, event):
        self.y = event.y_root - self.oldy
        self.x = event.x_root - self.oldx
        self.geometry(f'+{self.x}+{self.y}')

    def destroy_ind_modal(self, *args):
        self.gif.stop_thread()
        self.destroy()
        
class GifPlayer:
    def __init__(self,canvas,giffile,delay):
        self.frame=[]
        self.thread_working = True
        i=0
        while 1:
            try:
                image=tk.PhotoImage(file = giffile, format="gif -index "+str(i))
                self.frame.append(image)
                i=i+1
            except:
                break
        self.totalFrames=i-1
        self.delay=delay
        self.canvas = canvas
        # Create an image item on the canvas
        self.canvas.configure(image = self.frame[0])

    def play(self):
        """
        Plays the gif
        """
        self.gif_thread = _thread.start_new_thread(self.infinite,())

    def infinite(self):
        i = 0
        try:
            while self.thread_working == True:
                self.canvas.configure(image=self.frame[i])
                i = (i + 1) % self.totalFrames
                time.sleep(self.delay)
        except Exception:
            return

    def stop_thread(self):
        self.thread_working = False

class customCTkCheckBox(ctk.CTkCheckBox):
    def __init__(self, parent, text, **kwargs):
        super().__init__(parent , text=text,
                        checkbox_width = 20,
                        checkbox_height = 20, 
                        corner_radius = 4,
                        border_width = 2,
                        border_color = COLORS['neutral 2'],
                        hover_color = COLORS['neutral 3'],
                        text_color = COLORS['neutral 1'],
                        fg_color = COLORS['primary 1'],
                        font = ('Poppins Medium', 18),
                        **kwargs)
        #self._text_label.grid_forget()
        self._text_label.grid(row=0, column=2, sticky="w", pady=(2,0))

class ColorMapWidget(ctk.CTkFrame):
    def __init__(self, parent, colormap, root, **kwargs):
        super().__init__(parent,
                        fg_color='white',
                        border_width=0,
                        height=40, **kwargs)
        
        self.cmap = colormap.get()
        self.root = root
        self.update_cmap()

        self.left_label = ctk.CTkLabel(self, image=self.left_photo, text='')
        self.middle_label = ctk.CTkLabel(self, image=self.middle_photo, text='')
        self.right_label = ctk.CTkLabel(self, image=self.right_photo, text='')

        self.left_label.grid(row=0, column=0, sticky="nw", padx= 0)
        self.middle_label.grid(row=0, column=1, sticky="we", padx= 0)
        self.right_label.grid(row=0, column=2, sticky="ne", padx= 0)

        self.grid_columnconfigure(0, weight=0)  # Left part does not stretch
        self.grid_columnconfigure(1, weight=1)  # Middle part stretches
        self.grid_columnconfigure(2, weight=0)  # Right part does not stretch

        self.root.bind("<Configure>", self.resize, add='+')

    
    def update_cmap(self, refresh = False):
        # Load images according to the cmap thats selected
        directory = f"app\\core\\assets\\viewdata\\colormaps\\{self.cmap}\\"

        self.left_image = Image.open(directory+"left.png")
        self.middle_image = Image.open(directory+"middle.png")
        self.right_image = Image.open(directory+"right.png")

        # Convert images to PhotoImage for use in Label widgets
        self.left_photo = ctk.CTkImage(light_image = self.left_image, size=(15,41))
        self.middle_photo = ctk.CTkImage(light_image = self.middle_image, size=(423,41))
        self.right_photo = ctk.CTkImage(light_image = self.right_image, size=(15,41))

        if refresh:
            self.resize()

    def resize(self, *args):
        self.fixed_height = self.left_label.winfo_height()
        resized = self.middle_image.resize((self.middle_label.winfo_width()+2, self.fixed_height), Image.Resampling.BICUBIC)
        self.middle_photo = ImageTk.PhotoImage(resized)
        self.middle_label.configure(image=self.middle_photo)

class PlotTopLevel(ctk.CTkToplevel):
    def __init__(self, parent, name, geometry, img, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.geometry(geometry)
        self.title(name)
        self.attributes("-topmost", True)
        self.after(100, lambda: self.attributes("-topmost", False))

        self.name = name
        self.parent = parent
        self.img = img                                                                          #img contains the original image of the plot.
        self.tk_img = ImageTk.PhotoImage(img)
        self.canvas = ctk.CTkCanvas(self)
        self.zoom_options = ctk.StringVar(value = 'Fit to Frame')
        self.canvas_image = None

        self.canvas.pack(fill = 'both', expand=True)
        self.protocol("WM_DELETE_WINDOW", self.on_close)

        # Mouse Bindings
        self.canvas.bind('<ButtonPress-1>', lambda event: self.canvas.scan_mark(event.x, event.y))              #Pan function
        self.canvas.bind("<B1-Motion>", lambda event: self.canvas.scan_dragto(event.x, event.y, gain=1))
        self.canvas.bind("<Configure>", lambda event: self.zoom_event())

        self.construct()
        
    def on_close(self):
        for toplevel in self.parent.current_plot_toplevels:
            if toplevel == self:
                self.parent.current_plot_toplevels.remove(toplevel)
        self.destroy()

    def zoom_event(self, *args):
        if self.img:                                                                #Check if the original image exists            
            ow, oh = self.img.size                                                  #Original width and height of the plotted image
            aspect_ratio = ow / oh                                                  #The aspect ratio of the original image
            rw, rh = (self.canvas.winfo_width(), self.canvas.winfo_height())        #The frame size we need to fit the resized image into

            if self.zoom_options.get() == 'Fit to Frame':
                new_size = (int(rw), int(rw / aspect_ratio))                        #If fit to frame, calculate the appropriate size
            else:
                scale = int(self.zoom_options.get()[:-1])/100                       #Calculate the scale (1.5x, 2x, 2.5x....)
                rw, rh = (rw*scale, rh*scale)                                       #Now this is our new resized frame
                new_size = (int(rw), int(rw / aspect_ratio)) 
                    
            #Updating the canvas
            if 0 not in new_size:
                img = self.img.resize(new_size, Image.BOX)                                          #resizing the image
                self.tk_img = ImageTk.PhotoImage(img)                                               #set the new image as at tkinter image object
                self.canvas.itemconfig(self.canvas_image, image=self.tk_img, anchor='nw')           #Place it in the canvas
                self.canvas.configure(scrollregion=self.canvas.bbox('all'))                         #Set bounds so user can't pan outside the frame

    def show_imgplot(self):
        self.canvas_image = self.canvas.create_image(0, 0, image=self.tk_img, anchor='nw')
        self.canvas.configure(scrollregion=self.canvas.bbox('all'))
        self.zoom_event()

    def construct(self):
        self.show_imgplot()
        
class OtherButton(ctk.CTkButton):
    def __init__(self, parent, text, var, **kwargs):
        super().__init__(parent , text=text,
                        font = ('Poppins Medium', 18),
                        border_color=COLORS['border'],
                        fg_color=COLORS['white'],
                        corner_radius=16,
                        border_width=2,
                        text_color=COLORS['neutral 2'],
                        height=56,
                        command=self.toggle_boolvar,
                        hover=False,
                        anchor="w",
                        compound='right',
                        image = None,
                        **kwargs)
        
        self.toggle_var = var
        self.bind("<Enter>", self.on_enter)
        self.bind("<Leave>", self.on_leave)

    def toggle_boolvar(self, *args):
        current_value = self.toggle_var.get()
        self.toggle_var.set(not current_value)
        if self.toggle_var.get() is True:
            self.configure(border_color = COLORS['primary 1'],
                            text_color=COLORS['neutral 1'],
            )

        else:
            self.configure(border_color = COLORS['border'],
                            text_color=COLORS['neutral 2'],
            )

    def action_completed(self, *args):
            self.configure(image = ctk.CTkImage(light_image = Image.open('app\\core\\assets\\normalization\\check.png'),size=(18,18)), state = 'disabled')
            self.configure(state = 'normal')
            self.grid_columnconfigure(1, weight=1)
            self.grid_columnconfigure(2, weight=1)
            self._image_label.grid_configure(sticky='e', columnspan=2, padx=(0,20))
            self.toggle_boolvar()

    def on_enter(self, *args):
        if self.toggle_var.get() is False:
            self.configure(border_color = COLORS['primary 1'])

    def on_leave(self, *args):
        if self.toggle_var.get() is False:
            self.configure(border_color = COLORS['border'])
        

        
        

        
