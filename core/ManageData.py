#       Import Packages     #
import customtkinter as ctk
from customtkinter import filedialog
from core.assets.themecolors import COLORS
from core.customwidgets import *
from core.backend.managedata_backend import *
import threading
import numpy as np

#       ManageData Class     #
class ManageData(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent, fg_color=COLORS['background'])
        self.columnconfigure((0,1), weight=1, uniform=1)
        self.rowconfigure(0, weight=0)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=2)
        self.parent = parent
        self.root = self.parent.parent
        self.xenium_imgpath = None
        self.export_selections = set()

        #Variables to be used
        self.directory = ctk.StringVar("")
        self.data_type = ctk.StringVar(value='10X Visium')
                
        #Initialize the ManageData Class
        self.initialize()

        #   Bind Events     #
        self.upload_button.bind("<Configure>", self.upload_button_on_configure)

    def initialize(self):
        #   Top Bar  #
        topbar_frame = ctk.CTkFrame(self, fg_color='transparent')
        customTkLabel(topbar_frame, text='Manage Data', font=('Poppins Medium', 40), height=54, textMoveY= -22, text_color=COLORS['primary 1'], bg_color=COLORS['background']).pack(anchor = 'w', side = 'top')
        ctk.CTkLabel(topbar_frame, textvariable = self.parent.dataset_displayname, font = ('Poppins Medium', 18), text_color=COLORS['neutral 2']).pack(anchor = 'w', side = 'bottom', pady = 6)
        topbar_frame.grid(row = 0, column = 0, columnspan = 2, sticky = 'NSWE')

        #   Import your dataset     #
        midframe_1 = ctk.CTkFrame(self, fg_color='transparent')
        ctk.CTkLabel(midframe_1, text='Import your dataset', font = ('Poppins Medium', 20), text_color=COLORS['primary 1']).pack(anchor = 'w', side = 'top')        
        import_frame = ctk.CTkFrame(midframe_1, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])
        import_frame_1 = ctk.CTkFrame(import_frame, fg_color='transparent')
        ctk.CTkLabel(import_frame_1, text='Select the dataset type', font = ('Poppins Medium', 16), text_color=COLORS['neutral 2']).pack(anchor = 'w')
        customCTkOptionMenu(import_frame_1, values=['10X Visium', '10X Xenium', 'MERFISH', 'SeqFISH', 'SlideSeq', 'Adata'], variable = self.data_type).pack(anchor = 'w', fill = 'x', expand = True)
        ctk.CTkLabel(import_frame_1, text='Import your dataset', font = ('Poppins Medium', 16), text_color=COLORS['neutral 2']).pack(anchor = 'w', pady = (28, 0))
        self.upload_ico = ImageTk.PhotoImage(Image.open('app\\core\\assets\\upload.png').resize((16, 16), Image.Resampling.LANCZOS))
        self.upload_button = ctk.CTkButton(import_frame_1, height=56, fg_color=COLORS['white'], border_color=COLORS['border'], border_width=1, corner_radius=16, font=('Poppins Medium', 18), text_color=COLORS['primary 1'], text='Browse', image=self.upload_ico, hover_color=COLORS['background'], command=self.open_dialog)
        self.upload_button.pack(anchor = 'w', fill = 'x', expand = True)
        ctk.CTkButton(import_frame_1, width=662, height=64, fg_color=COLORS['primary 1'], corner_radius=16, font=('Poppins Semibold', 20), text_color=COLORS['white'], text='Import', hover_color=COLORS['primary 1 hover'], command=self.load_spatial_data).pack(anchor = 'center', pady = (40, 0))

        #   Export your dataset     #
        midframe_2 = ctk.CTkFrame(self, fg_color='transparent')
        ctk.CTkLabel(midframe_2, text='Export your dataset', font = ('Poppins Medium', 20), text_color=COLORS['primary 1']).pack(anchor = 'w', side = 'top')
        export_frame = ctk.CTkFrame(midframe_2, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'])
        self.export_frame_1 = ctk.CTkFrame(export_frame, fg_color='transparent')
        ctk.CTkLabel(self.export_frame_1, text='Select the files you want to export', font = ('Poppins Medium', 16), text_color=COLORS['neutral 2']).pack(anchor = 'w')
        self.export_datasets_frame = ctk.CTkScrollableFrame(self.export_frame_1, height=143, fg_color=COLORS['white'], border_color=COLORS['border'], border_width=1, corner_radius=16, scrollbar_button_color=COLORS['neutral 3'])
        self.export_datasets_frame._scrollbar.configure(height=0)
        self.export_datasets_frame._scrollbar.grid_configure(padx = 2)
        self.export_datasets_frame.pack(fill = 'both', expand = True)
        ctk.CTkButton(self.export_frame_1, width=662, height=64, fg_color=COLORS['primary 1'], corner_radius=16, font=('Poppins Semibold', 20), text_color=COLORS['white'], text='Export', hover_color=COLORS['primary 1 hover'], command=self.export).pack(pady = (39, 0))

        #   Disabled State  #
        self.dsb_export_frame = ctk.CTkFrame(export_frame, fg_color='transparent')
        dsb_ico = ImageTk.PhotoImage(Image.open('app\\core\\assets\\managedata\\waiting_for_import.png').resize((40, 36), Image.Resampling.LANCZOS))
        dsb_export_frame_centered = ctk.CTkFrame(self.dsb_export_frame, fg_color='transparent')
        ctk.CTkLabel(dsb_export_frame_centered, text="", image=dsb_ico).pack(pady = (0, 16))
        ctk.CTkLabel(dsb_export_frame_centered, text="Export options will become available after you've imported the dataset.", font = ('Poppins Medium', 16), text_color=COLORS['neutral 2']).pack()
        dsb_export_frame_centered.pack(expand = True)

        #   Summary     #
        midframe_3 = ctk.CTkFrame(self, fg_color='transparent')
        self.summary_frame = ctk.CTkFrame(midframe_3, fg_color=COLORS['white'], corner_radius=20, border_width=1, border_color=COLORS['border'], height = 342)
        ctk.CTkLabel(midframe_3, text='Summary', font = ('Poppins Medium', 20), text_color=COLORS['primary 1']).pack(anchor = 'w', side = 'top')

        #   Disabled State  #
        self.dsb_summary_frame = ctk.CTkFrame(self.summary_frame, fg_color='transparent')
        dsb_ico = ImageTk.PhotoImage(Image.open('app\\core\\assets\\managedata\\waiting_for_import.png').resize((40, 36), Image.Resampling.LANCZOS))
        ctk.CTkLabel(self.dsb_summary_frame, text="", image=dsb_ico).pack(pady = (0, 16))
        ctk.CTkLabel(self.dsb_summary_frame, text="Dataset Summary will be available following the import.", font = ('Poppins Medium', 16), text_color=COLORS['neutral 2']).pack()
        self.dsb_summary_frame.pack(expand = True)

        midframe_1.grid(row = 1, column = 0, pady = 24, padx = (0, 12), sticky='NSEW')
        midframe_2.grid(row = 1, column = 1, pady = 24, padx = (12, 0), sticky='NSEW')
        midframe_3.grid(row = 2, column = 0, columnspan = 2, sticky='NSWE')
        import_frame.pack(pady = (16, 0), fill = 'both', expand = True)
        import_frame_1.pack(pady = 40, padx = 40, fill = 'x', expand = True)
        export_frame.pack(pady = (16, 0), fill = 'both', expand = True)
        self.dsb_export_frame.pack(pady = 40, padx = 40, fill = 'both', expand = True)
        self.summary_frame.pack(pady = (16, 12), fill = 'both', expand = True)
        self.summary_frame.pack_propagate(0)

    def open_dialog(self):                                                                                  #File dialog asking the user to select the data folder.
        if self.data_type.get() in ['10X Visium', '10X Xenium', 'MERFISH', 'SlideSeq', 'SeqFISH']:
            path = filedialog.askdirectory(title='Select Output Folder')
            if path != '':
                self.directory.set(path)
                self.upload_button.configure(image= None, text='', font = ('Poppins', 16), text_color = COLORS['neutral 1'])
                self.upload_button_on_configure()
            else:
                self.directory.set('')
                self.upload_button.configure(image= self.upload_ico, text='Browse', font=('Poppins Medium', 18), text_color=COLORS['primary 1'])
        else:
            path = filedialog.askopenfilename(title="Select an Anndata file in .H5AD format",filetypes=[("H5AD files", "*.h5ad")])
            if path != '':
                self.directory.set(path)
                self.upload_button.configure(image= None, text='', font = ('Poppins', 16), text_color = COLORS['neutral 1'])
                self.upload_button_on_configure()
            else:
                self.directory.set('')
                self.upload_button.configure(image= self.upload_ico, text='Browse', font=('Poppins Medium', 18), text_color=COLORS['primary 1'])

    def upload_button_on_configure(self, *args):
        width = self.upload_button.winfo_width()
        if self.upload_button.cget('text') != 'Browse':
            # 0.1 characters fill in every pixel width of the button
            total_len = round((width * 0.1))
            if total_len >= len(self.directory.get()):
                self.upload_button.configure(text = self.directory.get())
            else:
                self.upload_button.configure(text = self.directory.get()[:total_len] + '...')

    def load_spatial_data(self):
        if self.parent.adata is None:
            self.root.lockAppWindow()
            if self.data_type.get() == '10X Xenium':
                response = self.root.build_msgmodal(msgbox_args = {'height' : 342, 
                                                                    'title':'Upload Image', 
                                                                    'message' : 'Would you like to import the Xenium Tissue image?'})
                if response == 'Cancel':
                    self.root.releaseAppWindow()
                    return
                elif response == 'Yes':
                    self.xenium_imgpath = filedialog.askopenfilename(filetypes=[("PNG files", "*.png")])
                    if self.xenium_imgpath is None:
                        self.root.releaseAppWindow()
                        return
                    
            self.root.build_IndProgressModal()
            secondary_thread = threading.Thread(target=read_10x, args=[self, self.data_type.get(), self.directory.get(), self.xenium_imgpath])
            secondary_thread.start()

    def refresh_statistics(self):
        self.refresh_export()
        if self.dsb_summary_frame.winfo_exists():
            self.dsb_summary_frame.destroy()
            self.adata_stats = ctk.CTkLabel(self.summary_frame, 
                                            text=f'AnnData object with n_obs × n_vars = {self.parent.total_spots.get()} × {self.parent.total_genes.get()}',
                                            font=('Poppins Medium', 18),
                                            text_color=COLORS['neutral 1'])
            self.adata_stats.pack(anchor = 'nw', pady = (16, 8), padx = 40)
            self.summary_frame_scroll = ctk.CTkScrollableFrame(self.summary_frame, 
                                                               fg_color=COLORS['white'], 
                                                               corner_radius=0, 
                                                               border_width=0, 
                                                               orientation='horizontal', 
                                                               scrollbar_button_color=COLORS['neutral 3'])
            self.summary_frame_scroll.pack(fill= 'both', expand = True, pady = (0, 16), padx = 34)

            #obs
            self.obs = ctk.CTkFrame(self.summary_frame_scroll, fg_color='transparent')
            ctk.CTkLabel(self.obs, text='obs:', font = ('Poppins Medium', 18), text_color=COLORS['neutral 1']).pack(side = 'left', padx = (0, 6), pady = 0)
            for col in self.parent.adata.obs:
                ctk.CTkLabel(self.obs, text=col, 
                            corner_radius=4, 
                            fg_color=COLORS['primary 2'], 
                            font=('Poppins Medium', 16), 
                            text_color=COLORS['primary 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            self.obs.pack(side = 'top', pady = 0, padx = 0, fill = 'x', anchor = 'nw')

            #var
            self.var = ctk.CTkFrame(self.summary_frame_scroll, fg_color='transparent')
            ctk.CTkLabel(self.var, text='var:', font = ('Poppins Medium', 18), text_color=COLORS['neutral 1']).pack(side = 'left', padx = (0, 6), pady = 0)
            for col in self.parent.adata.var:
                ctk.CTkLabel(self.var, text=col, 
                            corner_radius=4, 
                            fg_color=COLORS['primary 2'], 
                            font=('Poppins Medium', 16), 
                            text_color=COLORS['primary 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            self.var.pack(side = 'top', pady = (8, 0), padx = 0, fill = 'x', anchor = 'nw')

            #uns
            self.uns = ctk.CTkFrame(self.summary_frame_scroll, fg_color='transparent')
            ctk.CTkLabel(self.uns, text='uns:', font = ('Poppins Medium', 18), text_color=COLORS['neutral 1']).pack(side = 'left', padx = (0, 6), pady = 0)
            for col in self.parent.adata.uns:
                ctk.CTkLabel(self.uns, text=col, 
                             corner_radius=4, 
                             fg_color=COLORS['primary 2'], 
                             font=('Poppins Medium', 16), 
                             text_color=COLORS['primary 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            if len(self.parent.adata.uns) == 0:
                ctk.CTkLabel(self.uns, text='None', 
                            corner_radius=4, 
                            fg_color=COLORS['background'], 
                            font=('Poppins Medium', 16), 
                            text_color=COLORS['neutral 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            self.uns.pack(side = 'top', pady = (8, 0), padx = 0, fill = 'x', anchor = 'nw')

            #obsm
            self.obsm = ctk.CTkFrame(self.summary_frame_scroll, fg_color='transparent')
            ctk.CTkLabel(self.obsm, text='obsm:', font = ('Poppins Medium', 18), text_color=COLORS['neutral 1']).pack(side = 'left', padx = (0, 6), pady = 0)
            for col in self.parent.adata.obsm:
                ctk.CTkLabel(self.obsm, text=col, 
                             corner_radius=4, 
                             fg_color=COLORS['primary 2'], 
                             font=('Poppins Medium', 16), 
                             text_color=COLORS['primary 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            if len(self.parent.adata.obsm) == 0:
                ctk.CTkLabel(self.obsm, text='None', 
                            corner_radius=4, 
                            fg_color=COLORS['background'], 
                            font=('Poppins Medium', 16), 
                            text_color=COLORS['neutral 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            self.obsm.pack(side = 'top', pady = (8, 0), padx = 0, fill = 'x', anchor = 'nw')

            #varm
            self.varm = ctk.CTkFrame(self.summary_frame_scroll, fg_color='transparent')
            ctk.CTkLabel(self.varm, text='varm:', font = ('Poppins Medium', 18), text_color=COLORS['neutral 1']).pack(side = 'left', padx = (0, 6), pady = 0)
            for col in self.parent.adata.varm:
                ctk.CTkLabel(self.varm, text=col, 
                             corner_radius=4, 
                             fg_color=COLORS['primary 2'], 
                             font=('Poppins Medium', 16), 
                             text_color=COLORS['primary 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            if len(self.parent.adata.varm) == 0:
                ctk.CTkLabel(self.varm, text='None', 
                            corner_radius=4, 
                            fg_color=COLORS['background'], 
                            font=('Poppins Medium', 16), 
                            text_color=COLORS['neutral 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            self.varm.pack(side = 'top', pady = (8, 0), padx = 0, fill = 'x', anchor = 'nw')

            #obsp
            self.obsp = ctk.CTkFrame(self.summary_frame_scroll, fg_color='transparent')
            ctk.CTkLabel(self.obsp, text='obsp:', font = ('Poppins Medium', 18), text_color=COLORS['neutral 1']).pack(side = 'left', padx = (0, 6), pady = 0)
            for col in self.parent.adata.obsp:
                ctk.CTkLabel(self.obsp, text=col, 
                             corner_radius=4, 
                             fg_color=COLORS['primary 2'], 
                             font=('Poppins Medium', 16), 
                             text_color=COLORS['primary 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            if len(self.parent.adata.obsp) == 0:
                ctk.CTkLabel(self.obsp, text='None', 
                            corner_radius=4, 
                            fg_color=COLORS['background'], 
                            font=('Poppins Medium', 16), 
                            text_color=COLORS['neutral 1']).pack(side = 'left', ipadx = 4, ipady = 4, padx = 6, pady = 1)
            self.obsp.pack(side = 'top', pady = (8, 0), padx = 0, fill = 'x', anchor = 'nw')

    def refresh_export(self):
        def update_list(name, var):
            if var.get():
                self.export_selections.add(name)
            else:
                self.export_selections.discard(name)

        var = ctk.BooleanVar()
        if self.dsb_export_frame.winfo_exists():
            self.dsb_export_frame.destroy()
            self.export_frame_1.pack(pady = 40, padx = 40, fill = 'both', expand = True)
        customCTkCheckBox(self.export_datasets_frame, text='Anndata', variable=var, command=lambda n='Anndata', v=var: update_list(n, v)).pack(side = 'top', anchor = 'nw', padx = (16, 0), pady = (6, 0))
        if 'hires' in self.parent.adata.uns['spatial'][self.parent.dataname]['images']:
            var = ctk.BooleanVar()
            customCTkCheckBox(self.export_datasets_frame, text='Hires Tissue Image', variable=var, command=lambda n='Hires Tissue Image', v=var: update_list(n, v)).pack(side = 'top', anchor = 'nw', padx = (16, 0), pady = (6, 0))

    def export(self):
        self.root.lockAppWindow()
        if len(self.export_selections) == 0:
            info = msgbox(parent=self.root, callback = None, 
                            title='No Files Selected', 
                            cancel_button=False, 
                            yes_button=True, 
                            no_button=False, 
                            width=440, 
                            height=342, 
                            icon='info', 
                            message='Please select at least one file to proceed with the export.',
                            yes_text = 'Continue')
            self.root.wait_window(info)
            return
        
        path = filedialog.askdirectory(title='Select directory')
        if path == '':
            self.root.releaseAppWindow()
            return
        try:
            exception_occured = False
            process_anndata = False
            for select in self.export_selections:
                if select == 'Hires Tissue Image':
                    image = Image.fromarray((self.parent.adata.uns['spatial'][self.parent.dataname]['images']['hires']).astype(np.uint8), 'RGB')
                    image.save(f"{path}\\tissue_hires_image.png")

                ### ADD MORE FILES HERE             

                if select == 'Anndata':
                    self.root.build_IndProgressModal()
                    secondary_thread = threading.Thread(target=saveas, args=[self, path])
                    secondary_thread.start()
                    process_anndata = True

        except Exception as e:
            print(e)
            exception_occured = True
            fail = msgbox(parent=self.root, callback = None, 
                            title='Failed to export', 
                            cancel_button=False, 
                            yes_button=True, 
                            no_button=False, 
                            width=642, 
                            height=382, 
                            icon='fail', 
                            message='Failed to export one or more files. \n Please check the log for more information.',
                            yes_text = 'Try Again',
                            title_color=COLORS['error'])
            self.wait_window(fail)

            self.root.attributes('-topmost', True)
            self.root.attributes('-topmost', False)
            self.root.releaseAppWindow()
        
        finally:
            if exception_occured == False and process_anndata == False:
                success = msgbox(parent=self.root, callback = None, 
                title='Export Successful', 
                cancel_button=False, 
                yes_button=True, 
                no_button=False, 
                width=440, 
                height=342, 
                icon='success', 
                message='Successfully exported the selected files.',
                yes_text = 'Continue')
                self.wait_window(success)
                
                self.root.attributes('-topmost', True)
                self.root.attributes('-topmost', False)
                self.root.releaseAppWindow()
