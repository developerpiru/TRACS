import tkinter as tk
from tkinter import font as tkfont
from tkinter import *
from tkinter import filedialog
from tkinter.ttk import Frame, Label, Button, Style
import os
import subprocess
import sys
import pandas as pd
import numpy as np
from scipy import stats
import webbrowser
from datetime import datetime
import global_vars_v1 as global_vars


# version
app_version = "1.1.4"

# TRACS: Toolset for the Ranked Analysis of CRISPR Screens
# Created by Pirunthan Perampalam
# Please go to the official TRACS GitHub page for documentation and the latest updates
github_url = "https://github.com/developerpiru/TRACS"
# If you experience bugs/errors/other hurdles in your analyses, please leave a bug report on GitHub

# Data visualization and exploration for TRACS is enabled by VisualizeTRACS, a companion R shiny app
# You can launch VisualizeTRACS from its GitHub page
viz_url = "https://github.com/developerpiru/VisualizeTRACS"
# Help and documentation for VisualizeTRACS is available there

# Fixed lib log2FC calcs by addressing columns by NAME
# Fixed addressing columns by name for reg log2FC
# validated with Excel TRACS
# dynamic column naming and addressing working for all functions in read count processing
# implemented fixed replicate naming of -A, -B, -C, etc for replicates
# removed replicate naming stuff from window frame
# library trimming, alignment, read counts working
# adding input for library read file
# update summary working for all variables
# trimming with cutadapt working
# converting to fasta library working
# building bowtie2 index working
# implemented FILE_FLAGS[] global vars for all file attributes
# trimmed reads save in trimmed-reads folder
# bowtie2 alignments save in alignments folder
# processing for all initial and final replicate files working
# multicore method for cutadapt function working
# implemented replicate naming schema
# horizontal and vertical scrollbars for listboxes so filenames don't have to be split
# replicate naming and number working for cutadapt and bowtie2 functions based on replicate naming schema
# fixed bowtie2 alignment not working anymore
# working mageck count function
# added Cas9- integration
# loading Cas9- files into mageck working
# TRACS algorithm working parts:
# - DONE - get total read counts per sample
# - DONE - get lib log2FC for library
# - DONE - determine mean of replicates for log2FC
# - DONE - get log2FC per sample
# - DONE - get sgw per sample
# - DONE - get rank per sample
# - DONE - get Library, Initial, and Final ES per gene
# - DONE - get ES ratio (ER) (log2 fold change of Initial and Final ES) per gene
# TRACS algorithm functions work for n replicates
# improved statistics calculations for faster speed for p value and q values calculations
# fixed layout and UI in 1.0.2
# removed AnalyzeCounts function and bridged with TrimmingReads function
# cleaned summary box outputs version 1.1.1
# added option to set number of cpu cores
# version 1.1.3: use bowtie2 alignments to improve read depth for library reads only


class TRACSApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # load global variables
        global_vars.init_vars()

        # widget styles
        self.style = Style()
        #self.style.theme_use("default")
        self.frame_relief = FLAT
        self.frame_borderwidth = 0
        
        # app title
        self.winfo_toplevel().title("TRACS v" + app_version)

        # set fonts
        self.title_font = tkfont.Font(family='Arial', size=18, weight="bold")
        self.label_font = tkfont.Font(family='Arial', size=12, weight="bold")
        self.label_font_small = tkfont.Font(family='Arial', size=8, weight="bold")
        self.body_font = tkfont.Font(family='Arial', size=12, weight="normal")
        self.version_font = tkfont.Font(family='Arial', size=9, weight="normal")

        # main container
        container = tk.Frame(self, bg="#3d6fdb")
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # define page frames
        self.frames = {}
        # create frames for each page
        for F in (MainPage, NewExpSettings, LoadFiles, InitialCas9Pos, FinalCas9Pos, InitialCas9Neg, FinalCas9Neg, SummaryPage, TrimmingReads):
            page_name = F.__name__
            frame = F(parent=container, controller=self)
            self.frames[page_name] = frame

            # draw all frames in same location
            frame.grid(row=0, column=0, sticky="nsew")

        # load main page\frame at startup
        self.show_frame("MainPage")

    def show_frame(self, page_name):
        # brings selected frame to the top
        frame = self.frames[page_name]
        frame.tkraise()

    # function to select experiment directory
    def browse_folder(self, frame_name, target_name):
        target_obj = getattr(self.frames[frame_name], target_name)
        folder = filedialog.askdirectory(initialdir=os.getcwd(), title='Select folder...')
        target_obj.delete(0, tk.END)
        target_obj.insert(tk.END, folder)

    # function to select files
    def browse_file(self, multi_files, frame_name, target_name, file_type):
        target_obj = getattr(self.frames[frame_name], target_name)

        if file_type == "csv":
            file_options = (("CSV", "*.csv"),
                            ("All files", "*"))
        elif file_type == "fastq":
            file_options = (("FASTQ files", "*.fastq"),
                            ("GZipped Fastq", "*.fastq.gz"),
                            ("All files", "*"))

        if multi_files == "TRUE":
            files = filedialog.askopenfilenames(title='Select files...',
                                                filetypes=(file_options))
            for item in files:
                target_obj.insert(END, item)
        else:
            file = filedialog.askopenfilename(title='Select file...',
                                              filetypes=(file_options))
            target_obj.delete(0, END)
            target_obj.insert(END, file)

    # function to delete selected file from file list
    def delete_selected(self, frame_name, target_name):
        target_obj = getattr(self.frames[frame_name], target_name)
        selected = target_obj.curselection()
        for i in selected[::-1]:
            target_obj.delete(i)

    # function to get next replicate alphabetical character (e.g. A, B, C, etc.)
    # requires first replicate letter ('A') and incremental value (n)
    # returns correct replicate letter for the A+nth replicate
    def get_rep_char(self, n):
        start_char = "A"
        next_char = chr(ord(start_char) + int(n)).upper()
        return str(next_char)


class MainPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        frame = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        frame.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # info
        fr_1 = tk.Frame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = tk.Label(fr_1, text="TRACS: Toolset for the Ranked Analysis of CRISPR Screens",
                         wraplength="800", justify="left", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_info.pack(padx=5, pady=5)

        fr_2 = tk.Frame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_2.pack(fill=X, side=TOP, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = tk.Label(fr_2, text="For help and documentation, please visit our GitHub page at:"
                                    "\n" + github_url +
                                    " or click the button below to open it in your browser."
                                    "\n\n\nData visualization and exploration of data generated by TRACS is enabled by "
                                    "its companion R shiny app, VisualizeTRACS."
                                    "\nYou can access VisualizeTRACS at its "
                                    "own GitHub page at " + viz_url +
                                    "\n\n\nYou can find help and documentation for VisualizeTRACS there as well.",
                         wraplength="900", justify="left", font=controller.body_font, bg="#3d6fdb", fg="white")
        lbl_info.pack(side=TOP, padx=5, pady=5)
        lbl_version = tk.Label(fr_2, text="\n\n TRACS version: " + app_version,
                            wraplength="800", justify="left", font=controller.version_font, bg="#3d6fdb", fg="white")
        lbl_version.pack(side=TOP, padx=5, pady=5)

        btnNewExp = Button(frame, text="Start New Experiment", command=lambda: controller.show_frame("NewExpSettings"))
        btnNewExp.pack(padx=5, pady=5)

        btnLaunchGithub = Button(frame, text="Goto TRACS @ GitHub", command=lambda: self.openURL(github_url))
        btnLaunchGithub.pack(padx=5, pady=5)

        btnLaunchVisualizeTRACS = Button(frame, text="Goto VisualizeTRACS @ GitHub", command=lambda: self.openURL(viz_url))
        btnLaunchVisualizeTRACS.pack(padx=5, pady=5)

        self.pack(fill=BOTH, expand=True)

        btnBack = Button(self, text="Exit", command=lambda: app.destroy())
        btnBack.pack(side=RIGHT, padx=5, pady=5)

    def openURL(self, url):
        webbrowser.open(url, new=1)


class NewExpSettings(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="New experiment settings", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_exp_settings_Cas9_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # experiment name
        fr_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_exp_name = tk.Label(fr_1, text="Experiment name", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_exp_name.pack(side=LEFT, padx=5, pady=5)
        self.txt_exp_name = Entry(fr_1)
        self.txt_exp_name.pack(fill=X, padx=5, pady=5, expand=True)

        # output directory
        fr_2 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_2.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_expdir = tk.Label(fr_2, text="Experiment directory", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_expdir.pack(side=LEFT, anchor=N, padx=5, pady=5)
        self.txt_expdir = Entry(fr_2, width="50")
        self.txt_expdir.pack(fill=X, padx=5, pady=5, expand=True)
        btn_browse = Button(fr_2, text="Browse", command=lambda: self.controller.browse_folder(
            frame_name="NewExpSettings", target_name="txt_expdir"))
        btn_browse.pack(side=RIGHT, padx=5, pady=5)

        # library type
        fr_3 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_3.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)

        # name of initial group
        fr_4 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_4.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_name_initial = tk.Label(fr_4, text="Name of initial condition (T0)", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_name_initial.pack(side=LEFT, padx=5, pady=5)
        self.txt_name_initial = Entry(fr_4)
        self.txt_name_initial.pack(fill=X, padx=5, pady=5, expand=True)

        # name of final group
        fr_5 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_5.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_name_final = tk.Label(fr_5, text="Name of final condition (Tf)", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_name_final.pack(side=LEFT, padx=5, pady=5)
        self.txt_name_final = Entry(fr_5)
        self.txt_name_final.pack(fill=X, padx=5, pady=5, expand=True)

        # false discovery rate (FDR) for stats
        fr_6 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_6.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_fdr = tk.Label(fr_6, text="False discovery rate (e.g. 0.1 for 10%)", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_fdr.pack(side=LEFT, padx=5, pady=5)
        self.txt_fdr = Entry(fr_6)
        self.txt_fdr.pack(fill=X, padx=5, pady=5, expand=True)

        # number of CPU cores to use
        fr_7 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_7.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_cores = tk.Label(fr_7, text="Number of CPU cores to use", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_cores.pack(side=LEFT, padx=5, pady=5)
        self.txt_cores = Entry(fr_7)
        self.txt_cores.pack(fill=X, padx=5, pady=5, expand=True)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: self.controller.show_frame("LoadFiles"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("MainPage"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


class LoadFiles(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="Initial pooled sgRNA library files", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_load_files_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # reference library
        fr_2 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_2.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_reflibrary = tk.Label(fr_2, text="Library reference file (CSV)", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_reflibrary.pack(side=LEFT, anchor=N, padx=5, pady=5)
        self.txt_reflibrary = tk.Entry(fr_2)
        self.txt_reflibrary.pack(fill=X, padx=5, pady=5, expand=True)
        btn_browse = Button(fr_2, text="Browse", command=lambda: self.controller.browse_file(multi_files="FALSE",
                                                                                             frame_name="LoadFiles",
                                                                                             target_name="txt_reflibrary",
                                                                                             file_type="csv"))
        btn_browse.pack(side=RIGHT, padx=5, pady=5)

        # library reads file
        fr_6 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_6.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_libreads = tk.Label(fr_6, text="Initial library (L0) read file (FASTQ)", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_libreads.pack(side=LEFT, anchor=N, padx=5, pady=5)
        self.txt_libreads = tk.Entry(fr_6)
        self.txt_libreads.pack(fill=X, padx=5, pady=5, expand=True)
        btn_browse = Button(fr_6, text="Browse", command=lambda: self.controller.browse_file(multi_files="FALSE",
                                                                                             frame_name="LoadFiles",
                                                                                             target_name="txt_libreads",
                                                                                             file_type="fastq"))
        btn_browse.pack(side=RIGHT, padx=5, pady=5)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("InitialCas9Pos"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("NewExpSettings"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


class InitialCas9Pos(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="Load initial Cas9-positive files", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_load_files_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # container frame for listbox labels
        fr_3 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_3.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial label
        fr_4 = tk.Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_4.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_i_info = tk.Label(fr_4, text="Please select your initial condition (T0) Cas9-positive read files (FASTQ)", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_i_info.pack(padx=5, pady=5, side=LEFT)

        # container frame for list boxes
        list_container = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        list_container.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial listbox
        initial_list_container = tk.Frame(list_container, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#FFFFFF")
        initial_list_container.pack(fill=BOTH, padx=0, pady=0, side=LEFT, expand=True)

        # frame for initial listbox buttons
        initial_btn_container = tk.Frame(initial_list_container, bg="#3d6fdb")
        initial_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # initial listbox
        self.lst_initial_files = Listbox(initial_list_container, height="18", width="900", selectmode=EXTENDED)

        # y scroll for initial listbox
        yscrl_lst_initial_files = Scrollbar(initial_list_container, orient=VERTICAL)
        yscrl_lst_initial_files.pack(side=RIGHT, fill=Y)
        yscrl_lst_initial_files.config(command=self.lst_initial_files.yview)
        self.lst_initial_files.config(yscrollcommand=yscrl_lst_initial_files.set)

        # x scroll for initial listbox
        xscrl_lst_initial_files = Scrollbar(initial_list_container, orient=HORIZONTAL)
        xscrl_lst_initial_files.pack(side=BOTTOM, fill=X)
        xscrl_lst_initial_files.config(command=self.lst_initial_files.xview)
        self.lst_initial_files.config(xscrollcommand=xscrl_lst_initial_files.set)

        self.lst_initial_files.pack(side=LEFT, padx=0, pady=0)

        # add
        btn_i_add = Button(initial_btn_container, text="Add", command=lambda: self.controller.browse_file(multi_files="TRUE",
                                                                                                          frame_name="InitialCas9Pos",
                                                                                                          target_name="lst_initial_files",
                                                                                                          file_type="fastq"))
        btn_i_add.pack(side=TOP, padx=0, pady=0)
        # remove
        btn_i_rem = Button(initial_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="InitialCas9Pos",
                                                                           target_name="lst_initial_files"))
        btn_i_rem.pack(side=TOP, padx=0, pady=0)
        # clear
        btn_i_clear = Button(initial_btn_container, text="Clear", command=lambda: self.lst_initial_files.delete(0, END))
        btn_i_clear.pack(side=TOP, padx=0, pady=0)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("FinalCas9Pos"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("LoadFiles"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


class FinalCas9Pos(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="Load final Cas9-positive files", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_load_files_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # container frame for listbox labels
        fr_3 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_3.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for final label
        fr_5 = tk.Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_5.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_f_info = tk.Label(fr_5, text="Please select your final condition (T0) Cas9-positive read files (FASTQ)", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_f_info.pack(padx=5, pady=5, side=LEFT)

        # container frame for final list boxes
        list_container = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        list_container.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for final listbox
        final_list_container = tk.Frame(list_container, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#FFFFFF")
        final_list_container.pack(fill=BOTH, padx=0, pady=0, side=LEFT, expand=True)

        # frame for final listbox buttons
        final_btn_container = tk.Frame(final_list_container, bg="#3d6fdb")
        final_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # final listbox
        self.lst_final_files = Listbox(final_list_container, height="18", width="900", selectmode=EXTENDED)

        # y scroll for final listbox
        yscrl_lst_final_files = Scrollbar(final_list_container, orient=VERTICAL)
        yscrl_lst_final_files.pack(side=RIGHT, fill=Y)
        yscrl_lst_final_files.config(command=self.lst_final_files.yview)
        self.lst_final_files.config(yscrollcommand=yscrl_lst_final_files.set)
        # x scroll for final listbox
        xscrl_lst_final_files = Scrollbar(final_list_container, orient=HORIZONTAL)
        xscrl_lst_final_files.pack(side=BOTTOM, fill=X)
        xscrl_lst_final_files.config(command=self.lst_final_files.xview)
        self.lst_final_files.config(xscrollcommand=xscrl_lst_final_files.set)

        self.lst_final_files.pack(side=LEFT, padx=0, pady=0)

        # add
        btn_f_add = Button(final_btn_container, text="Add",
                           command=lambda: self.controller.browse_file(multi_files="TRUE",
                                                                       frame_name="FinalCas9Pos",
                                                                       target_name="lst_final_files",
                                                                       file_type="fastq"))
        btn_f_add.pack(side=TOP, padx=0, pady=0)
        # remove
        btn_f_rem = Button(final_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="FinalCas9Pos",
                                                                           target_name="lst_final_files"))
        btn_f_rem.pack(side=TOP, padx=0, pady=0)
        # clear
        btn_f_clear = Button(final_btn_container, text="Clear", command=lambda: self.lst_final_files.delete(0, END))
        btn_f_clear.pack(side=TOP, padx=0, pady=0)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("InitialCas9Neg"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("InitialCas9Pos"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


# function to load Cas9- files
class InitialCas9Neg(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="Load initial Cas9-negative files", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_load_files_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # container frame for listbox labels
        fr_3 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_3.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial label
        fr_4 = tk.Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_4.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_i_info = tk.Label(fr_4, text="Please select your initial condition (T0) Cas9-negative read files (FASTQ)", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_i_info.pack(padx=5, pady=5, side=LEFT)

        # container frame for list boxes
        list_container = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        list_container.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial listbox
        initial_list_container = tk.Frame(list_container, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#FFFFFF")
        initial_list_container.pack(fill=BOTH, padx=0, pady=0, side=LEFT, expand=True)

        # frame for initial listbox buttons
        initial_btn_container = tk.Frame(initial_list_container, bg="#3d6fdb")
        initial_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # initial listbox
        self.lst_initial_files = Listbox(initial_list_container, height="18", width="900", selectmode=EXTENDED)

        # y scroll for initial listbox
        yscrl_lst_initial_files = Scrollbar(initial_list_container, orient=VERTICAL)
        yscrl_lst_initial_files.pack(side=RIGHT, fill=Y)
        yscrl_lst_initial_files.config(command=self.lst_initial_files.yview)
        self.lst_initial_files.config(yscrollcommand=yscrl_lst_initial_files.set)

        # x scroll for initial listbox
        xscrl_lst_initial_files = Scrollbar(initial_list_container, orient=HORIZONTAL)
        xscrl_lst_initial_files.pack(side=BOTTOM, fill=X)
        xscrl_lst_initial_files.config(command=self.lst_initial_files.xview)
        self.lst_initial_files.config(xscrollcommand=xscrl_lst_initial_files.set)

        self.lst_initial_files.pack(side=LEFT, padx=0, pady=0)

        # add
        btn_i_add = Button(initial_btn_container, text="Add",
                           command=lambda: self.controller.browse_file(multi_files="TRUE",
                                                                       frame_name="InitialCas9Neg",
                                                                       target_name="lst_initial_files",
                                                                       file_type="fastq"))
        btn_i_add.pack(side=TOP, padx=0, pady=0)
        # remove
        btn_i_rem = Button(initial_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="InitialCas9Neg",
                                                                           target_name="lst_initial_files"))
        btn_i_rem.pack(side=TOP, padx=0, pady=0)
        # clear
        btn_i_clear = Button(initial_btn_container, text="Clear", command=lambda: self.lst_initial_files.delete(0, END))
        btn_i_clear.pack(side=TOP, padx=0, pady=0)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("FinalCas9Neg"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("FinalCas9Pos"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


class FinalCas9Neg(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="Load final Cas9-negative files", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_load_files_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # container frame for listbox labels
        fr_3 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_3.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for final label
        fr_5 = tk.Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_5.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_f_info = tk.Label(fr_5, text="Please select your final condition (T0) Cas9-negative read files (FASTQ)", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_f_info.pack(padx=5, pady=5, side=LEFT)

        # container frame for final list boxes
        list_container = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        list_container.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for final listbox
        final_list_container = tk.Frame(list_container, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#FFFFFF")
        final_list_container.pack(fill=BOTH, padx=0, pady=0, side=LEFT, expand=True)

        # frame for final listbox buttons
        final_btn_container = tk.Frame(final_list_container, bg="#3d6fdb")
        final_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # final listbox
        self.lst_final_files = Listbox(final_list_container, height="18", width="900", selectmode=EXTENDED)

        # y scroll for final listbox
        yscrl_lst_final_files = Scrollbar(final_list_container, orient=VERTICAL)
        yscrl_lst_final_files.pack(side=RIGHT, fill=Y)
        yscrl_lst_final_files.config(command=self.lst_final_files.yview)
        self.lst_final_files.config(yscrollcommand=yscrl_lst_final_files.set)
        # x scroll for final listbox
        xscrl_lst_final_files = Scrollbar(final_list_container, orient=HORIZONTAL)
        xscrl_lst_final_files.pack(side=BOTTOM, fill=X)
        xscrl_lst_final_files.config(command=self.lst_final_files.xview)
        self.lst_final_files.config(xscrollcommand=xscrl_lst_final_files.set)

        self.lst_final_files.pack(side=LEFT, padx=0, pady=0)

        # add
        btn_f_add = Button(final_btn_container, text="Add",
                           command=lambda: self.controller.browse_file(multi_files="TRUE",
                                                                       frame_name="FinalCas9Neg",
                                                                       target_name="lst_final_files",
                                                                       file_type="fastq"))
        btn_f_add.pack(side=TOP, padx=0, pady=0)
        # remove
        btn_f_rem = Button(final_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="FinalCas9Neg",
                                                                           target_name="lst_final_files"))
        btn_f_rem.pack(side=TOP, padx=0, pady=0)
        # clear
        btn_f_clear = Button(final_btn_container, text="Clear", command=lambda: self.lst_final_files.delete(0, END))
        btn_f_clear.pack(side=TOP, padx=0, pady=0)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: self.process_inputs())
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("InitialCas9Neg"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    def process_inputs(self):

        # load all settings into global variables
        global_vars.EXPERIMENT_SETTINGS['Experiment name'] = \
            self.controller.frames['NewExpSettings'].txt_exp_name.get()  # experiment name
        global_vars.EXPERIMENT_SETTINGS['Experiment directory'] = \
            self.controller.frames['NewExpSettings'].txt_expdir.get()  # experiment directory
        global_vars.EXPERIMENT_SETTINGS['Initial condition name'] = \
            self.controller.frames['NewExpSettings'].txt_name_initial.get()  # initial condition name (day 0)
        global_vars.EXPERIMENT_SETTINGS['Final condition name'] = \
            self.controller.frames['NewExpSettings'].txt_name_final.get()  # final condition name (day X)
        global_vars.EXPERIMENT_SETTINGS['FDR'] = \
            self.controller.frames['NewExpSettings'].txt_fdr.get()  # FDR
        global_vars.EXPERIMENT_SETTINGS['CPU cores'] = \
            self.controller.frames['NewExpSettings'].txt_cores.get()  # CPU cores
        global_vars.EXPERIMENT_SETTINGS['Library reference file'] = \
            self.controller.frames['LoadFiles'].txt_reflibrary.get()  # library reference file

        # get library read fastq file
        global_vars.LIBRARY_READ_PATH = self.controller.frames['LoadFiles'].txt_libreads.get()

        print("\nLibrary read path: %s" % global_vars.LIBRARY_READ_PATH)

        # get Initial condition Cas9+ fastq files
        global_vars.LIST_INITIAL_PATHS = self.controller.frames['InitialCas9Pos'].lst_initial_files.get(0, END)

        # get Final condition Cas9+ fastq files
        global_vars.LIST_FINAL_PATHS = self.controller.frames['FinalCas9Pos'].lst_final_files.get(0, END)

        # get Initial condition Cas9- fastq files
        global_vars.LIST_CAS9_NEG_I_PATHS = self.controller.frames['InitialCas9Neg'].lst_initial_files.get(0, END)

        # get Final condition Cas9- fastq files
        global_vars.LIST_CAS9_NEG_F_PATHS = self.controller.frames['FinalCas9Neg'].lst_final_files.get(0, END)

        # reset status box
        summary_box = self.controller.frames['SummaryPage'].txt_summary
        summary_box.config(state=NORMAL)
        summary_box.delete(1.0, END)

        # print experiment settings into status text box
        for keys, values in global_vars.EXPERIMENT_SETTINGS.items():
            tmp_text = "%s: %s\n" % (keys, values)
            summary_box.insert(END, tmp_text)

        # print library read file manually since it's only a single file
        summary_box.insert(END, "\n\nInitial Library read file:")
        summary_box.insert(END, "\n\n" + global_vars.LIBRARY_READ_PATH)

        # call function to print files into status box
        self.print_files(tar_name=global_vars.LIST_INITIAL_PATHS,
                         condition=global_vars.EXPERIMENT_SETTINGS['Initial condition name'], cas9_type="-positive")
        self.print_files(tar_name=global_vars.LIST_FINAL_PATHS,
                         condition=global_vars.EXPERIMENT_SETTINGS['Final condition name'], cas9_type="-positive")
        self.print_files(tar_name=global_vars.LIST_CAS9_NEG_I_PATHS,
                         condition=global_vars.EXPERIMENT_SETTINGS['Initial condition name'], cas9_type="-negative")
        self.print_files(tar_name=global_vars.LIST_CAS9_NEG_F_PATHS,
                         condition=global_vars.EXPERIMENT_SETTINGS['Final condition name'], cas9_type="-negative")

        summary_box.config(state=DISABLED)

        # show summary page after processing all inputs
        self.controller.show_frame("SummaryPage")

    # function to print selected files into status box
    def print_files(self, tar_name, condition, cas9_type):
        summary_box = self.controller.frames['SummaryPage'].txt_summary

        tmp_text = "\n\nSample files for Cas9%s '%s' condition:\n" % (cas9_type, condition)

        summary_box.insert(END, tmp_text)

        for i in tar_name:
            tmp_text = "\nFile %s: %s" % (tar_name.index(i) + 1, i)
            summary_box.insert(END, tmp_text)


class SummaryPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="Summary", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_exp_settings_Cas9_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # experiment name
        fr_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = tk.Label(fr_1, text="Please confirm all details below before continuing:", wraplength="900", justify="left", font=controller.label_font, bg="#3d6fdb", fg="white")
        lbl_info.pack(fill=X, side=LEFT, padx=5, pady=5)

        # summary text box
        fr_2 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#FFFFFF")
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        self.txt_summary = Text(fr_2, height="18")
        self.txt_summary.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        scrl_txt_summary = Scrollbar(fr_2)
        scrl_txt_summary.pack(side=RIGHT, fill=Y)
        scrl_txt_summary.config(command=self.txt_summary.yview)
        self.txt_summary.config(yscrollcommand=scrl_txt_summary.set)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("TrimmingReads"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("FinalCas9Neg"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

        btn_Back.pack(side=RIGHT, padx=5, pady=5)


class TrimmingReads(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="Run analysis", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_exp_settings_Cas9_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # experiment name
        fr_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = tk.Label(fr_1, text="You are now ready to start analysis! You have two options:"
                                    "\n\n1. If you are starting from raw read files to get read counts, select"
                                    " option 1. TRACS will automatically do the following:"
                                    "\n\n-Trim your reads\n-Build a reference library\n-Perform alignments"
                                    "\n-Generate a read counts file\n-Run the TRACS algorithm to determine"
                                    " gene essentiality."
                                    "\n\n2. If you have already generated a read counts file in the correct format"
                                    " and have it saved in the 'readcounts' folder in your experiment directory, select"
                                    " option 2. Using your read counts file, TRACS will automatically:"
                                    "\n\n-Run the TRACS algorithm to determine gene essentiality.",
                         wraplength="900", justify="left", font=controller.body_font, bg="#3d6fdb", fg="white")
        lbl_info.pack(fill=X, side=LEFT, padx=5, pady=5)

        # summary text box
        fr_2 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)

        # option 1 button: run the entire pipline
        btnTrim = Button(fr_2, text="1. I have raw reads (FASTQ files). Do everything for me!",
                         command=lambda: self.do_steps_option_1())
        btnTrim.pack(side=TOP)

        # option 2 button: only run the TRACS algorithm
        btnTest = Button(fr_2, text="2. I already have a read counts file (TXT). Just give me the essential genes!",
                         command=lambda: self.do_steps_option_2())
        btnTest.pack(side=TOP)

        # navigation buttons
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("SummaryPage"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    # function to save experiment settings to file
    def save_exp_settngs(self):

        summary_box = self.controller.frames['SummaryPage'].txt_summary

        file_name = global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "-run_settings.txt"

        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        f = open(file_name, 'w')
        f.write("TRACS version: %s\n" % app_version)
        f.write("\nDate: %s\n\n\n" % str(datetime.now()))
        f.write(summary_box.get(1.0, END))
        f.close()

    def do_steps_option_1(self):

        # save experiment settings to file
        self.save_exp_settngs()

        # call function to convert to fasta library
        self.convert_to_fasta()

        # trim reads
        print("\nTrimming reads...")

        # call do_trimming_v2 function to trim reads for all samples
        # Library reads
        self.do_trimming_v2(file_list=global_vars.LIBRARY_READ_PATH,
                            condition_name="Library",
                            cas9_type="")

        # all Cas9+ samples
        self.do_trimming_v2(file_list=global_vars.LIST_INITIAL_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[0])
        self.do_trimming_v2(file_list=global_vars.LIST_FINAL_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[0])

        # all Cas9- samples
        self.do_trimming_v2(file_list=global_vars.LIST_CAS9_NEG_I_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[1])
        self.do_trimming_v2(file_list=global_vars.LIST_CAS9_NEG_F_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[1])

        # build index
        self.build_index()

        # align reads
        print("\nAligning reads with library index...")

        # call do_bowtie2_alignment function to align the library reads using bowtie2
        # Library reads
        self.do_bowtie2_alignment(file_list=global_vars.LIBRARY_READ_PATH,
                                  condition_name="Library",
                                  cas9_type="")

        # run mageck_count function to get raw read counts
        self.do_mageck_count(file_list_1=global_vars.LIST_INITIAL_PATHS,
                             file_list_2=global_vars.LIST_FINAL_PATHS,
                             condition_name_1=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                             condition_name_2=global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        print("\n\nDone alignments and read count generation!")
        print("\n\nStarting TRACS algorithm...")

        # start do_steps_option_2 function
        self.do_steps_option_2()

    # function to call cutadapt to trim provided read sample files
    def do_trimming_v2(self, file_list, condition_name, cas9_type):

        # save experiment settings to file
        self.save_exp_settngs()

        # change directory
        trimmed_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                           global_vars.FILE_FLAGS['Trimmed dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(trimmed_target_path) is False:
            print("\nCreating directory for trimmed files...")
            print("\n%s" % trimmed_target_path)
            os.mkdir(global_vars.FILE_FLAGS['Trimmed dir'])

        os.chdir(trimmed_target_path)

        for i in file_list:

            if condition_name == "Library":
                # process for library read file
                this_filename = condition_name + global_vars.FILE_FLAGS['Trimmed read file']
                readfile = global_vars.LIBRARY_READ_PATH
            else:
                # process for initial and final reads files
                # get replicate name
                next_char = self.controller.get_rep_char(file_list.index(i))
                this_filename = condition_name + cas9_type + next_char + global_vars.FILE_FLAGS['Trimmed read file']
                readfile = i

            # determine line count in file to overcome multicore ban in cutadapt
            line_count = (os.popen('wc -l \'%s\'' % readfile).read()).split(" ", 1)

            # build command pieces
            cmd_prefix = "head \'%s\' -n %s| " % (readfile, line_count[0])
            cmd_main = "cutadapt -g NAAGGACGAAACANNN...GTTTTAGAGCTAGAAAT -> "
            cmd_f_out = "\'" + this_filename + "\' "
            cmd_args = "-e 0.1 --cores " + global_vars.EXPERIMENT_SETTINGS['CPU cores'] + " --trimmed-only"

            cmd = cmd_prefix + cmd_main + cmd_f_out + cmd_args

            print("\n\nSaving as filename: %s" % this_filename)

            # run command in shell
            cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            show_in_console = True

            # output
            while True:
                out = cmd_output.stderr.read(1)
                if out == b'' and cmd_output.poll() != None:
                    break
                if out != '':
                    output = out.decode('utf-8')
                    if show_in_console:
                        sys.stdout.write(output)
                        sys.stdout.flush()

            # check if we are doing the steps for the Library read file, in which case only one file is required
            # if true, then break the for loop from continuing
            # otherwise for loop continues to iterate over other replicates (for initial and final conditions)
            if condition_name == "Library":
                break

    # function to convert CSV reference library file to fasta format for use with bowtie2
    def convert_to_fasta(self):

        print("\n>Creating directory for reference library...")

        # change directory
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        # build command pieces
        cmd_prefix = "awk -F \',\' \'{print \">\"$1\"\\n\"$2}\' "
        cmd_f_in = "\'" + global_vars.EXPERIMENT_SETTINGS['Library reference file'] + "\'"
        cmd_f_out = " > \'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + \
                    global_vars.FILE_FLAGS['Fasta library file'] + "\'"
        cmd = cmd_prefix + cmd_f_in + cmd_f_out

        print("\nConverting reference library file from CSV to FASTA...")

        # run command
        cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        show_in_console = True

        # output
        while True:
            out = cmd_output.stderr.read(1)
            if out == b'' and cmd_output.poll() != None:
                break
            if out != '':
                output = out.decode('utf-8')
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

        print("\nFinished converting reference library file!")

    # function to generate bowtie2 index for reference library
    def build_index(self):

        print("\nBuilding library index...")

        # change directory
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        global_vars.FILE_FLAGS['Bowtie2 index full path'] = \
            os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                         global_vars.FILE_FLAGS['Bowtie2 index name'] + "/")

        if os.path.exists(global_vars.FILE_FLAGS['Bowtie2 index full path']) is False:
            print("\n\'bowtie2-index\' directory not found; creating it now...")
            os.mkdir(global_vars.FILE_FLAGS['Bowtie2 index name'])

        os.chdir(global_vars.FILE_FLAGS['Bowtie2 index name'])

        # build command pieces
        # structure:    bowtie2-build LibraryAB_fasta.fa bowtie2_index_LibraryAB
        cmd_prefix = "bowtie2-build "
        cmd_f_in = "\'../" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + \
                   global_vars.FILE_FLAGS['Fasta library file'] + "\' "
        cmd_f_out = "\'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "-" + \
                    global_vars.FILE_FLAGS['Bowtie2 index name'] + "\'"
        cmd = cmd_prefix + cmd_f_in + cmd_f_out

        # run command
        cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        show_in_console = True

        # output
        while True:
            out = cmd_output.stderr.read(1)
            if out == b'' and cmd_output.poll() != None:
                break
            if out != '':
                output = out.decode('utf-8')
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

    # function to perform alignments for provided sample read files using bowtie2
    def do_bowtie2_alignment(self, file_list, condition_name, cas9_type):

        # change directory
        alignment_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                             global_vars.FILE_FLAGS['Alignments dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(alignment_target_path) is False:
            print("\nCreating directory for aligned files...")
            print("\n%s" % alignment_target_path)
            os.mkdir(global_vars.FILE_FLAGS['Alignments dir'])

        trimmed_read_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                         global_vars.FILE_FLAGS['Trimmed dir'])

        os.chdir(alignment_target_path)

        # loop through file list for current condition (Library, Initial, Final conditions)
        for i in file_list:
            # update status box
            print("\nGetting alignment data for initial library reads %s: %s" % (file_list.index(i) + 1, i + global_vars.FILE_FLAGS['Trimmed read file']))

            # check if current specified condition is the Library condition - special treatment if true
            if condition_name == "Library":
                # process for library read file
                input_file = os.path.join(trimmed_read_path, condition_name + global_vars.FILE_FLAGS['Trimmed read file'])
                output_file = condition_name + global_vars.FILE_FLAGS['Aligned bam file']
            else:
                # process for initial and final sample read files
                # get replicate name
                next_char = self.controller.get_rep_char(file_list.index(i))
                input_file = os.path.join(trimmed_read_path, condition_name + cas9_type + next_char + global_vars.FILE_FLAGS['Trimmed read file'])
                output_file = condition_name + cas9_type + next_char + global_vars.FILE_FLAGS['Aligned bam file']

            # build command pieces
            cmd_prefix = "bowtie2 -x \'%s" % global_vars.FILE_FLAGS['Bowtie2 index full path'] + \
                         global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "-" + global_vars.FILE_FLAGS[
                             'Bowtie2 index name'] + "\' -U "
            cmd_f_in = "\'" + input_file + "\'"
            cmd_f_out = "\'" + output_file + "\'"
            cmd_suffix = " --norc -N 1 -p %i | samtools view -bS - > %s" % (4*int(global_vars.EXPERIMENT_SETTINGS['CPU cores']), cmd_f_out)
            cmd = cmd_prefix + cmd_f_in + cmd_suffix

            # run command
            cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            show_in_console = True

            # output
            while True:
                out = cmd_output.stderr.read(1)
                if out == b'' and cmd_output.poll() != None:
                    break
                if out != '':
                    output = out.decode('utf-8')
                    if show_in_console:
                        sys.stdout.write(output)
                        sys.stdout.flush()

            # check if we are doing the steps for the Library read file,
            # if true, then break the for loop from continuing because it only has one replicate; no looping required
            if condition_name == "Library":
                break

    # function to run mageck's count function to get read counts from alignments for all samples
    def do_mageck_count(self, file_list_1, file_list_2, condition_name_1, condition_name_2):

        # change directory
        counts_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                             global_vars.FILE_FLAGS['Readcounts dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(counts_target_path) is False:
            print("\nCreating directory for read count files...")
            print("\n%s" % counts_target_path)
            os.mkdir(global_vars.FILE_FLAGS['Readcounts dir'])

        aligned_read_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                         global_vars.FILE_FLAGS['Alignments dir'])  # + "/"

        os.chdir(counts_target_path)

        cmd_sample_list = ""
        cmd_input_file_list = ""

        # prepare library read file first
        input_file = os.path.join(aligned_read_path, "Library" + global_vars.FILE_FLAGS['Aligned bam file'])
        # add Library name to sample list
        cmd_sample_list = cmd_sample_list + "Library" + ","
        # add trimmed bam-aligned Library read file to input file list
        cmd_input_file_list = cmd_input_file_list + "\'" + input_file + "\' "

        # loop for all Cas9+ and Cas9- sample read files
        for j in global_vars.CAS9_TYPE_NAMING:

            # initial condition labels and files
            for i in file_list_1:
                print("\nNow processing file %s: %s" % (file_list_1.index(i) + 1, i + global_vars.FILE_FLAGS['Trimmed read file']))

                # get replicate name
                next_char = self.controller.get_rep_char(file_list_1.index(i))

                input_file = os.path.join(aligned_read_path, condition_name_1 + j + next_char + global_vars.FILE_FLAGS['Trimmed read file'])

                cmd_sample_list = cmd_sample_list + condition_name_1 + j + next_char + ","

                cmd_input_file_list = cmd_input_file_list + "\'" + input_file + "\' "

            # final condition labels and files
            for i in file_list_2:
                print("\nNow processing file %s: %s" % (file_list_2.index(i) + 1, i + global_vars.FILE_FLAGS['Trimmed read file']))

                # get replicate name
                next_char = self.controller.get_rep_char(file_list_2.index(i))

                input_file = os.path.join(aligned_read_path, condition_name_2 + j + next_char + global_vars.FILE_FLAGS['Trimmed read file'])

                cmd_sample_list = cmd_sample_list + condition_name_2 + j + next_char + ","

                cmd_input_file_list = cmd_input_file_list + "\'" + input_file + "\' "

        # build command pieces
        cmd_prefix = "mageck count -l \'" + global_vars.EXPERIMENT_SETTINGS['Library reference file'] + "\'"
        cmd_namtag = " -n \'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "\'"
        cmd_sampletags = " --sample-label \'" + cmd_sample_list[:-1] + "\'"
        cmd_readfiles = " --fastq " + cmd_input_file_list
        cmd = cmd_prefix + cmd_namtag + cmd_sampletags + cmd_readfiles

        # run command
        cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
        show_in_console = True

        # output
        while True:
            out = cmd_output.stderr.read(1)
            if out == b'' and cmd_output.poll() != None:
                break
            if out != '':
                output = out.decode('utf-8')
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

    def do_steps_option_2(self):

        # run function to analyze reads with TRACS algorithm
        self.process_reads()

    # function to run TRACS algorithm to calculate Enrichment Scores and Enrichment Ratio from read counts table
    def process_reads(self):

        # change directory
        reads_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                         global_vars.FILE_FLAGS['Readcounts dir'] + "/")

        os.chdir(reads_target_path)

        # get name of read counts file
        readsfile = global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".count.txt"

        # load counts into pandas dataframe
        self.raw_counts = pd.read_csv(readsfile, delimiter='\t')

        # set gene names as index
        self.raw_counts.index = self.raw_counts['Gene']

        # remove gene and sgRNA columns
        del self.raw_counts['Gene'], self.raw_counts['sgRNA']

        # sort by gene name
        self.raw_counts = self.raw_counts.sort_index()

        # determine total number of replicates
        self.num_replicates = (len(self.raw_counts.columns) - 1) / 4

        # print total number of replicates
        print("\nNumber of replicates: %i" % self.num_replicates)

        # # # NORMALIZING READS # # #
        print("\nNormalizing read counts...")

        # first add 1 to each read count to prevent divide by 0 errors when calculating log2FC
        self.raw_counts += 1

        # get sum of each column (total reads for each sample)
        total_reads = self.raw_counts.sum(axis=0)

        # normalize each column/sample by the total reads for each column
        i = 0
        for col in self.raw_counts.columns:
            # create new column name
            name = col + global_vars.COLUMN_NAMES['Normalized']

            # normalize
            self.raw_counts[name] = (self.raw_counts[col] * 1.0) / total_reads[i]

            # format for correct digits
            self.raw_counts[name] = self.raw_counts[name].map(lambda x: '%.16f' % x)
            self.raw_counts[name] = self.raw_counts[name].astype(float)
            i += 1

        # # # CALCULATE LOG2 FOLD CHANGE # # #
        print("\nCalculating log2 fold changes...")

        # save original column names
        orig_column_names = self.raw_counts.columns

        # log2FC structure: log2FC = [final reads / initial reads]
        # call function to calculate log2FC for library-Cas9-neg-initial [initial Cas9 neg / library]
        self.lib_log2FC_calculator()

        # call function to calculate log2FC for initial condition [initial Cas9 pos / initial Cas9 neg ]
        self.log2FC_calculator(name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'])

        # call function to calculate log2FC for final condition [final Cas9 pos / final Cas9 neg ]
        self.log2FC_calculator(name = global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        # # # GENE SCORE # # #
        print("\nCalculating gene scores...")

        # determine number of sgRNAs for each gene
        self.num_sgRNAs = pd.DataFrame()
        self.num_sgRNAs['num sgRNA'] = self.raw_counts.groupby(self.raw_counts.index)[orig_column_names[0]].count()

        print("\nDetermining average...")

        # call function to calculate average of log2FC for every gene for each log2FC sample column
        # calculates a single log2FC from all the guides for each gene
        self.avg = pd.DataFrame()  # new dataframe to hold averages
        self.avg_calculator(name = global_vars.COLUMN_NAMES['Library'])
        self.avg_calculator(name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'])
        self.avg_calculator(name = global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        # # # ENRICHMENT SCORE # # #
        print("\nCalculating sgw...")

        # determine sgw for each gene
        self.sgw = pd.DataFrame()  # new dataframe to hold sgw
        self.get_sgw()

        print("\nCalculating rank...")
        # determine rank for each sample
        self.get_rank()

        # save to file
        print("\nSaving to weighted ranks (sgw) to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.sgw.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".sgw.csv", index=True,
                        sep=',',
                        float_format='%.16f')

        print("\nCalculating enrichment scores...")
        # determine enrichment score for each gene
        self.ES = pd.DataFrame()

        # copy sgRNA number into ES dataframe
        self.ES['Num.sgRNA'] = self.num_sgRNAs['num sgRNA']

        # new ES method
        # call function to get Library.ES
        self.get_ES(name = global_vars.COLUMN_NAMES['Library'])

        # call function to get Initial.ES
        self.get_ES(name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'])

        # call function to get Final.ES
        self.get_ES(name = global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        # save to file
        print("\nSaving enrichment scores to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.ES.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".ES.csv", index=True,
                               sep=',',
                               float_format='%.16f')

        # # # ENRICHMENT RATIO # # #
        print("\nCalculating enrichment ratio...")
        # determine enrichment ratio (ER); log2 ratio of ES ( log2(final ES / initial ES)
        initial_col_name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'] + global_vars.COLUMN_NAMES['ES']
        final_col_name = global_vars.EXPERIMENT_SETTINGS['Final condition name'] + global_vars.COLUMN_NAMES['ES']
        self.ES[global_vars.COLUMN_NAMES['ER']] = np.log2(self.ES[final_col_name].astype(float) / self.ES[initial_col_name].astype(float))
        self.ES[global_vars.COLUMN_NAMES['ER']] = self.ES[global_vars.COLUMN_NAMES['ER']].map(lambda x: '%.6f' % x)

        # # # CALCULATE STATISTICS # # #
        print("\nCalculating statistics...")
        # create temp_stats dataframe
        self.temp_stats = pd.DataFrame()
        # call function to calculate statistics
        self.get_stats(initial_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                       final_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        # copy stats to end of ES dataframe
        self.ES['pval'] = self.temp_stats['pvalue']
        self.ES['qval'] = self.temp_stats['qvalue']

        # rename columns in final ES dataframe
        self.ES = self.ES.rename(columns={initial_col_name: global_vars.COLUMN_NAMES['Initial'], final_col_name: global_vars.COLUMN_NAMES['Final']})
        
        # save to file
        print("\nSaving final table to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.ES.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".csv",
                       index=True, sep=',', float_format='%.6f')

        print("\nAll done!")

    # function to calculate log2FC for Library. Log2FC = log2[Initial Cas9- / Library]
    def lib_log2FC_calculator(self):

        # loop through all initial Cas9- replicates
        for col in range(0, int(self.num_replicates)):
            # get name of library column
            library_norm_name = global_vars.COLUMN_NAMES['Library'] + global_vars.COLUMN_NAMES['Normalized']

            # get name of current initial cas9- replicate
            current_col_name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'] + global_vars.COLUMN_NAMES['Cas9 neg'] + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['Normalized']

            # create new column name
            new_col_name = global_vars.COLUMN_NAMES['Library'] + "-" + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['log2FC']

            # get log2FC for current replicate
            self.raw_counts[new_col_name] = np.log2(self.raw_counts[current_col_name].astype(float) / self.raw_counts[library_norm_name].astype(float))

    # function to calculate log2FC for Initial and Final condition
    # condition is determined by 'name' variable
    # log2FC = log2[condition Cas9+ / condition Cas9-]
    def log2FC_calculator(self, name):

        # temp dataframe to hold mean values
        temp_mean = pd.DataFrame()

        # copy correct columns (denominator) into temp dataframe so mean function can be used across columns
        for col in range(0, int(self.num_replicates)):
            current_col_name = name + global_vars.COLUMN_NAMES['Cas9 neg'] + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['Normalized']
            temp_mean[current_col_name] = self.raw_counts[current_col_name]

        # calculate mean across columns in temp dataframe
        temp_mean['mean'] = temp_mean.mean(axis=1)

        # loop through condition's Cas9+ columns (either initial or final; determined by 'name' variable and find log2FC
        for col in range(0, int(self.num_replicates)):
            # get current column name
            current_col_name = name + global_vars.COLUMN_NAMES['Cas9 pos'] + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['Normalized']

            # create new column name
            new_col_name = name + "-" + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['log2FC']

            # calculate log2FC for current replicate
            self.raw_counts[new_col_name] = np.log2(self.raw_counts[current_col_name].astype(float) / temp_mean['mean'].astype(float))

    # function to group all guides for a single together and determine the average log2FC for those guides
    # e.g. 6 sgRNA ID log2FC --> reduced to single log2FC per gene
    # current condition determined by 'name' variable
    def avg_calculator(self, name):

        # determine average of all guides per gene for current condition
        for col in range(0, int(self.num_replicates)):
            # get current column name
            current_col_name = name + "-" + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['log2FC']

            # create new column name
            new_col_name = name + "-" + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['average']

            # determine average log2FC for all guides per gene
            self.avg[new_col_name] = self.raw_counts.groupby(self.raw_counts.index)[current_col_name].mean()

    # function to determine sgRNA weighted score (sgw)
    # works across all columns in self.avg dataframe
    def get_sgw(self):

        # get column names from self.avg dataframe
        column_names = self.avg.columns

        # loop through all columns in self.avg dataframe
        for col in range(0, len(self.avg.columns)):
            # get current column name
            current_name = column_names[col]

            # create new column name
            new_col_name = current_name[:-4] + global_vars.COLUMN_NAMES['sgw']  # remove 4 characters from the end (removes ".avg")

            # determine sgw: sgRNA weighted score (log2FC * # of sgRNAs)
            self.sgw[new_col_name] = self.avg[current_name].astype(float) * self.num_sgRNAs['num sgRNA'].astype(float)

    # function to determine rank of each gene
    # sorts self.sgw dataframe columns from smallest to largest (ascending order)
    # ranks genes from 1 to x, where x is the highest value
    def get_rank(self):

        # get column names from self.sgw dataframe
        column_names = self.sgw.columns

        # loop through all columns in self.sgw dataframe
        for col in range(0, len(self.sgw.columns)):
            # get current column name
            current_name = column_names[col]

            # create new column name
            new_col_name = current_name[:-4] + global_vars.COLUMN_NAMES['rank']  # remove 4 characters from the end (removes ".sgw")

            # sort values by ascending order in current_name column
            self.sgw = self.sgw.sort_values(current_name)

            # determine rank of gene
            # increment the rank from 1 to x
            self.sgw[new_col_name] = pd.RangeIndex(stop=self.sgw.shape[0]) + 1

    # function to calculate pvalues and qvalues for each gene using get_rank() integer values
    def get_stats(self, initial_name, final_name):

        # temp dataframe to hold mean values for initial and final conditions
        temp_mean_initial = pd.DataFrame()
        temp_mean_final = pd.DataFrame()

        # copy columns of stated condition (determined by 'name') into temp dataframe
        for col in range(0, int(self.num_replicates)):
            # get current column initial name
            current_initial_col_name = initial_name + "-" + self.controller.get_rep_char(col) + \
                                       global_vars.COLUMN_NAMES['rank']
            current_final_col_name = final_name + "-" + self.controller.get_rep_char(col) + \
                                       global_vars.COLUMN_NAMES['rank']

            # copy column into temp_mean so mean can easily be calculated
            temp_mean_initial[current_initial_col_name] = self.sgw[current_initial_col_name]
            temp_mean_final[current_final_col_name] = self.sgw[current_final_col_name]

        # calculate the mean for initial and final conditions
        self.temp_stats['i_mean'] = temp_mean_initial.mean(axis=1)
        self.temp_stats['f_mean'] = temp_mean_final.mean(axis=1)

        # calculate standard deviation for initial and final conditions
        self.temp_stats['i_std'] = temp_mean_initial.apply(lambda x: stats.tstd(x), axis=1)
        self.temp_stats['f_std'] = temp_mean_final.apply(lambda x: stats.tstd(x), axis=1)

        # calculate p values
        self.temp_stats['pvalue'] = self.temp_stats.apply(lambda x: stats.ttest_ind_from_stats(mean1=x['i_mean'],
                                                          std1=x['i_std'],
                                                          nobs1=int(self.num_replicates),
                                                          mean2=x['f_mean'],
                                                          std2=x['f_std'],
                                                          nobs2=int(self.num_replicates),
                                                          equal_var=False)[1], axis=1)

        # sort values from smallest to largest p value
        self.temp_stats = self.temp_stats.sort_values('pvalue')

        # determine rank of p value
        # increment the rank from 1 to x
        self.temp_stats['i'] = pd.RangeIndex(stop=self.temp_stats.shape[0]) + 1

        # calculate q values
        self.temp_stats['qvalue'] = self.temp_stats.apply(
            lambda x: (x['i'] / len(self.temp_stats.index)) * float(global_vars.EXPERIMENT_SETTINGS['FDR']), axis=1)

        # format decimals\digits
        self.temp_stats['pvalue'] = self.temp_stats['pvalue'].map(lambda x: '%.16f' % x)
        self.temp_stats['qvalue'] = self.temp_stats['qvalue'].map(lambda x: '%.16f' % x)

    # function to determine Enrichment Score (ES)
    # current condition determined by 'name' variable
    def get_ES(self, name):

        # temp dataframe to hold mean values
        temp_mean = pd.DataFrame()

        # copy columns of stated condition (determined by 'name') into temp dataframe
        for col in range(0, int(self.num_replicates)):
            # get current column name
            current_col_name = name + "-" + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['rank']

            # copy column into temp_mean so mean can easily be calculated
            temp_mean[current_col_name] = self.sgw[current_col_name]

        # calculate mean across columns in temp dataframe
        temp_mean['mean'] = temp_mean.mean(axis=1)

        # create new column name
        new_col_name = name + global_vars.COLUMN_NAMES['ES']

        # calculate ES for current condition
        self.ES[new_col_name] = temp_mean['mean'] / self.num_sgRNAs['num sgRNA'].astype(float)


if __name__ == "__main__":
    app = TRACSApp()
    app.geometry("940x680+500+50")
    app.configure(bg="red")
    app.resizable(False, False)
    app.mainloop()



