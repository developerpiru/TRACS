# version 1.1.19 2023


import tkinter as tk
from tkinter import font as tkfont
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter.ttk import Button, Progressbar, Style
import os
import subprocess
import sys
import pandas as pd
import numpy as np
from scipy import stats
import webbrowser
from datetime import datetime
import _thread
import global_vars_v1 as global_vars

# version
app_version = "1.1.8"

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
# version 1.1.4: new theme, bug fixes
# version 1.1.4: added error checking for all user inputs
# version 1.1.5: bug fixes
# version 1.1.6: added error checking for all analysis steps; TRACS stops on any error now
# added new frame to show once TRACS completes
# TrimmingReads function/frame renamed to StartAnalysis
# version 1.1.7: added logging
# version 1.1.8: added threaded progress bar and status updates


def openURL(url):
    webbrowser.open(url, new=1)

# function to write messages to log file and console
def write_to_log(message):

    # log file path
    file_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                             global_vars.EXPERIMENT_SETTINGS['Experiment name'] +
                             global_vars.FILE_FLAGS['Log file'])

    # write to console
    # print(message)

    # show message in status label
    # global_vars.status_label.config(text="%s\n%s" % (global_vars.status_label.cget("text"), message))

    # show message in status text box
    global_vars.status_label.config(state="normal")
    global_vars.status_label.insert(END, "\n%s" % message)

    # scroll down on the text box
    global_vars.status_label.see("end")

    # disable status text box
    global_vars.status_label.config(state=DISABLED)

    try:
        # write to log file
        file = open(file_path, 'a+')
        file.write("\n[%s] %s" % (str(datetime.now()), message))
        file.close()
    except IOError:
        print("\nWarning: Cannot write to log file!")


# function to error check multi-file input for all read files
# make sure at least 2 files are selected
def check_read_inputs(source_object, condition, controller, next_frame, self):

    if not source_object.get(1, END):  # the first file path is 0, second is 1
        messagebox.showerror(title="Error",
                             message="You must select at least 2 FASTQ read files for your %s" % condition)
    else:
        # if there were no errors, continue to the next frame/page/step
        # check to see if the passed frame contains PROCESSINPUTS message, in which case it is the final file-loading
        # frame and the FinalCas9Neg.process_inputs() should be run
        if next_frame == "PROCESSINPUTS":
            FinalCas9Neg.process_inputs(self)
        else:
            controller.show_frame(next_frame)


# function to check if the file located at file_path is accessible and be read
def check_file_access(file_path):
    # returns True only if file can be accessed
    try:
        file = open(file_path, 'r')
        file.close()
    except IOError:
        return False
    return True


# main app class
class TRACSApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # load global variables
        global_vars.init_vars()

        # widget styles
        self.style = Style()
        self.style.theme_use("clam")
        # self.style.configure("progressbar_style", foreground="#0DD9A2", background="#0DD9A2")
        self.frame_relief = FLAT
        self.frame_borderwidth = 0
        
        # app title
        self.winfo_toplevel().title("TRACS " + app_version)

        # set fonts
        self.title_font = tkfont.Font(family='Arial', size=18, weight="bold")
        self.label_font = tkfont.Font(family='Arial', size=12, weight="bold")
        self.label_font_small = tkfont.Font(family='Arial', size=8, weight="bold")
        self.intro_body_font = tkfont.Font(family='Arial', size=12, weight="normal")
        self.status_label_font = tkfont.Font(family='Arial', size=12, weight="bold")
        self.version_font = tkfont.Font(family='Arial', size=9, weight="normal")

        # main container
        container = tk.Frame(self, bg="#3d6fdb")
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # define page frames
        self.frames = {}
        # create frames for each page
        for F in (MainPage, NewExpSettings, LoadFiles, InitialCas9Pos, FinalCas9Pos, InitialCas9Neg, FinalCas9Neg,
                  SummaryPage, StartAnalysis, EndAnalysis):
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


# frame start page
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
                                    "\nYou can download VisualizeTRACS at its "
                                    "own GitHub page at " + viz_url +
                                    "\n\n\nYou can find help and documentation for VisualizeTRACS there as well.",
                         wraplength="900", justify="left", font=controller.intro_body_font, bg="#3d6fdb", fg="white")
        lbl_info.pack(side=TOP, padx=5, pady=5)
        lbl_version = tk.Label(fr_2, text="\n\n TRACS version: " + app_version,
                            wraplength="800", justify="left", font=controller.version_font, bg="#3d6fdb", fg="white")
        lbl_version.pack(side=TOP, padx=5, pady=5)

        btnNewExp = Button(frame, text="Start New Experiment", command=lambda: controller.show_frame("NewExpSettings"))
        btnNewExp.pack(padx=5, pady=5)

        # btnLaunchGithub = Button(frame, text="Goto TRACS @ GitHub", command=lambda: openURL(github_url))
        # btnLaunchGithub.pack(padx=5, pady=5)
        #
        # btnLaunchVisualizeTRACS = Button(frame, text="Goto VisualizeTRACS @ GitHub", command=lambda: openURL(viz_url))
        # btnLaunchVisualizeTRACS.pack(padx=5, pady=5)

        self.pack(fill=BOTH, expand=True)

        btnBack = Button(self, text="Exit", command=lambda: app.destroy())
        btnBack.pack(side=RIGHT, padx=5, pady=5)
    #
    # def openURL(self, url):
    #     webbrowser.open(url, new=1)


# frame for experiment settings
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
        lbl_expdir = tk.Label(fr_2, text="Experiment folder", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
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
        lbl_fdr = tk.Label(fr_6, text="False discovery rate (e.g. 0.05 for 5%)", width="40", font=controller.label_font, bg="#3d6fdb", fg="white")
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
        btn_Next = Button(self, text="Next", command=lambda: self.check_inputs())
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("MainPage"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    # function to error check inputs from NewExpSettings and save values to global variables
    def check_inputs(self):

        # variable to record user input errors
        error_log = ""

        # check each user input and save into global variables

        # experiment name
        if self.controller.frames['NewExpSettings'].txt_exp_name.get() != "":
            global_vars.EXPERIMENT_SETTINGS['Experiment name'] = \
                self.controller.frames['NewExpSettings'].txt_exp_name.get()  # experiment name
        else:
            error_log = error_log + "Experiment name cannot be blank"

        # experiment folder
        if self.controller.frames['NewExpSettings'].txt_expdir.get() != "":
            global_vars.EXPERIMENT_SETTINGS['Experiment directory'] = \
                self.controller.frames['NewExpSettings'].txt_expdir.get()  # experiment directory
        else:
            error_log = error_log + "\n\nExperiment folder cannot be blank"

        # initial condition name
        if self.controller.frames['NewExpSettings'].txt_name_initial.get() != "":
            global_vars.EXPERIMENT_SETTINGS['Initial condition name'] = \
                self.controller.frames['NewExpSettings'].txt_name_initial.get()  # initial condition name (day 0)
        else:
            error_log = error_log + "\n\nInitial condition name cannot be blank"

        # final condition name
        if self.controller.frames['NewExpSettings'].txt_name_final.get() != "":
            global_vars.EXPERIMENT_SETTINGS['Final condition name'] = \
                self.controller.frames['NewExpSettings'].txt_name_final.get()  # final condition name (day X)
        else:
            error_log = error_log + "\n\nFinal condition name cannot be blank"

        # FDR
        if self.controller.frames['NewExpSettings'].txt_fdr.get() != "":
            try:
                global_vars.EXPERIMENT_SETTINGS['FDR'] = \
                    float(int(self.controller.frames['NewExpSettings'].txt_fdr.get()) / 100)  # FDR
            except ValueError:
                global_vars.EXPERIMENT_SETTINGS['FDR'] = 0
        else:
            error_log = error_log + "\n\nFalse discovery rate cannot be blank"

        # CPU cores
        if self.controller.frames['NewExpSettings'].txt_cores.get() != "":
            try:
                global_vars.EXPERIMENT_SETTINGS['CPU cores'] = \
                    int(self.controller.frames['NewExpSettings'].txt_cores.get())  # CPU cores
            except ValueError:
                global_vars.EXPERIMENT_SETTINGS['CPU cores'] = 0
        else:
            error_log = error_log + "\n\nCPU cores cannot be blank"

        # additionally, check to make sure FDR and CPU cores have numerical values

        # FDR
        if self.controller.frames['NewExpSettings'].txt_fdr.get() != "":
            if not self.controller.frames['NewExpSettings'].txt_fdr.get().isdigit():
                error_log = error_log + "\n\nFalse discovery rate must be a number"
        # CPU cores
        if self.controller.frames['NewExpSettings'].txt_cores.get() != "":
            if not self.controller.frames['NewExpSettings'].txt_cores.get().isdigit():
                error_log = error_log + "\n\nNumber of CPU cores must be a number"

        # check condition names and make sure they are unique (can't have same name for both)
        if self.controller.frames['NewExpSettings'].txt_name_initial.get() != "" and \
                self.controller.frames['NewExpSettings'].txt_name_initial.get() == self.controller.frames['NewExpSettings'].txt_name_final.get():
            error_log = error_log + "\n\nInitial and final condition names must be different"

        # if error_log logged any errors, display a warning message with the errors
        if error_log != "":
            messagebox.showerror(title="Error", message="%s" % error_log)
        else:
            # if there were no errors, continue to the next frame/page/step
            self.controller.show_frame("LoadFiles")


# frame for loading library files
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
        btn_Next = Button(self, text="Next", command=lambda: self.check_inputs())
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("NewExpSettings"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    # function to error check inputs from LoadFiles and save values to global variables
    def check_inputs(self):

        error_log = ""
        # check user input and save into global variables
        # sgRNA library csv file
        if self.controller.frames['LoadFiles'].txt_reflibrary.get() != "":
            global_vars.EXPERIMENT_SETTINGS['Library reference file'] = \
                self.controller.frames['LoadFiles'].txt_reflibrary.get()  # library reference file
        else:
            error_log = error_log + "\nLibrary reference file (CSV)"

        # initial pooled sgRNA library reads fastq file
        if self.controller.frames['LoadFiles'].txt_libreads.get() != "":
            global_vars.LIBRARY_READ_PATH = self.controller.frames['LoadFiles'].txt_libreads.get()
        else:
            error_log = error_log + "\nInitial library (L0) read file (FASTQ)"

        # if error_log logged any errors, then display a warning message with the errors
        if error_log != "":
            messagebox.showerror(title="Error",
                                 message="You must select files for the following: \n%s" % error_log)
        else:
            # if there were no errors, continue to the next frame/page/step
            self.controller.show_frame("InitialCas9Pos")


# frame for loading initial Cas9 positive files
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
        btn_Next = Button(self, text="Next", command=lambda: check_read_inputs(
            source_object=self.controller.frames['InitialCas9Pos'].lst_initial_files,
            condition="initial Cas9-positive (T0) condition",
            controller=self.controller,
            next_frame="FinalCas9Pos",
            self=self))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("LoadFiles"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


# frame for loading final Cas9 positive files
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
        btn_Next = Button(self, text="Next", command=lambda: check_read_inputs(
            source_object=self.controller.frames['FinalCas9Pos'].lst_final_files,
            condition="final Cas9-positive (Tf) condition",
            controller=self.controller,
            next_frame="InitialCas9Neg",
            self=self))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("InitialCas9Pos"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


# frame for loading initial Cas9 negative files
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
        btn_Next = Button(self, text="Next", command=lambda: check_read_inputs(
            source_object=self.controller.frames['InitialCas9Neg'].lst_initial_files,
            condition="initial Cas9-negative (T0) condition",
            controller=self.controller,
            next_frame="FinalCas9Neg",
            self=self))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("FinalCas9Pos"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


# frame for loading final Cas9 negative files
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
        btn_Next = Button(self, text="Next", command=lambda: check_read_inputs(
            source_object=self.controller.frames['FinalCas9Neg'].lst_final_files,
            condition="final Cas9-negative (Tf) condition",
            controller=self.controller,
            next_frame="PROCESSINPUTS",
            self=self))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("InitialCas9Neg"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    # function to verify if the number of replicates provided for each condition matches
    def check_replicates(self):

        if len(global_vars.LIST_INITIAL_PATHS) and len(global_vars.LIST_FINAL_PATHS) and len(global_vars.LIST_CAS9_NEG_I_PATHS) == len(global_vars.LIST_CAS9_NEG_F_PATHS):
            return True
        else:
            return False

    def process_inputs(self):

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
        summary_box.insert(END, "\n\nInitial Library read file: %s" % global_vars.LIBRARY_READ_PATH)

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

        # check if the number of files supplied for each condition matches
        result = self.check_replicates()
        if result is False:
            msg = "\nYou have not provided the same number of samples/replicates for each condition. " \
                  "\nPlease ensure the number of samples provided for Initial and Final conditions for Cas9-positive " \
                  "and Cas9-negative sets are equal!" \
                  "\nStopping TRACS."
            print(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return
        else:
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

        # label
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
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("StartAnalysis"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("FinalCas9Neg"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

        btn_Back.pack(side=RIGHT, padx=5, pady=5)


class StartAnalysis(tk.Frame):

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

        # start button
        # fr_2 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        # fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        # global_vars.btnStart = Button(fr_2, text="Start",
        #                  command=lambda: _thread.start_new_thread(self.run_analysis_steps, ()))
        # global_vars.btnStart.pack(side=TOP)

        # visualize button
        # global_vars.btnVisualize = Button(fr_2, text="Open in VisualizeTRACS",
        #                               command=lambda: _thread.start_new_thread(self.run_analysis_steps, ()))
        # global_vars.btnVisualize.pack(side=LEFT)
        #
        # # disable visualize button
        # global_vars.btnVisualize.config(state=DISABLED)

        # progress bar
        fr_progress = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth,
                               bg="#3d6fdb")
        fr_progress.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        global_vars.progress_bar = Progressbar(fr_progress, orient="horizontal", mode="determinate", length="900")
        global_vars.progress_bar.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # initial values for progress bar
        global_vars.progress_bar['maximum'] = 100
        global_vars.progress_bar['value'] = 0

        # text box for status updates
        fr_2 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#FFFFFF")
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        global_vars.status_label = Text(fr_2, height="20")
        global_vars.status_label.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        scrl_txt_summary = Scrollbar(fr_2)
        scrl_txt_summary.pack(side=RIGHT, fill=Y)
        scrl_txt_summary.config(command=global_vars.status_label.yview)
        global_vars.status_label.config(yscrollcommand=scrl_txt_summary.set)

        # navigation buttons
        # global_vars.btnVisualize = Button(self, text="Open in VisualizeTRACS", command=lambda: controller.show_frame("StartAnalysis"))
        # global_vars.btnVisualize.pack(side=RIGHT, padx=5, pady=5)
        # global_vars.btnVisualize.config(state=DISABLED)

        global_vars.btnStart = Button(self, text="Start", command=lambda: _thread.start_new_thread(self.run_analysis_steps, ()))
        global_vars.btnStart.pack(side=RIGHT, padx=5, pady=5)

        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("SummaryPage"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

        # global_vars.btnStart = Button(self, text="Start",
        #                               command=lambda: _thread.start_new_thread(self.run_analysis_steps, ()))
        # global_vars.btnStart.pack(side=TOP)

        # _thread.start_new_thread(self.run_analysis_steps, ())

    # function to update progress bar
    def update_progressbar(self):

        global_vars.current_step += 1
        global_vars.progress_bar['value'] = (global_vars.current_step/global_vars.total_steps)*100

    # function to save experiment settings to file
    def save_exp_settngs(self):

        summary_box = self.controller.frames['SummaryPage'].txt_summary

        file_name = global_vars.EXPERIMENT_SETTINGS['Experiment name'] + global_vars.FILE_FLAGS['Run settings']

        # first chance to create experiment directory if it doesn't exist yet
        if os.path.exists(global_vars.EXPERIMENT_SETTINGS['Experiment directory']) is False:
            write_to_log("Experiment directory not found...")
            try:
                write_to_log("Creating it now...")
                os.mkdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
            except IOError or FileNotFoundError:
                write_to_log("Could not create experiment directory!")
                return False

        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        try:
            write_to_log("Saving experiment settings to file")
            f = open(file_name, 'w')
            f.write("TRACS version: %s\n" % app_version)
            f.write("\nDate: %s\n\n\n" % str(datetime.now()))
            f.write(summary_box.get(1.0, END))
            f.close()
        except IOError:
            write_to_log("Could not save experiment settings to file! "
                         "Please check the permissions for your experiment folder!")
            return False
        return True

    # function that calls individual functions for each step
    # function is called to run in a background process/thread
    def run_analysis_steps(self):

        # disable Start button
        global_vars.btnStart.config(state=DISABLED)

        # calculate the total steps for the experiment because it varies depending on the number of replicates
        # 16 fixed steps (including library trimming) + 4*number of replicates
        global_vars.total_steps = 16 + 4 * len(global_vars.LIST_INITIAL_PATHS)

        # save experiment settings to file
        result = self.save_exp_settngs()
        if result is False:
            msg = "Error encountered saving experiment settings. " \
                  "\nPlease make sure you have permission to write to the experiment folder!" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # update progress bar
        self.update_progressbar()

        # call function to convert to fasta library
        result = self.convert_to_fasta()
        if result is False:
            msg = "Error encountered during FASTA conversion." \
                  "\nPlease make sure you have permission to write to the experiment folder!" \
                  "\nCould not create FASTA file.\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # update progress bar
        self.update_progressbar()

        # trim reads
        write_to_log("Trimming reads...")

        # call do_trimming_v2 function to trim reads for all samples
        # Library reads
        result = self.do_trimming_v2(file_list=global_vars.LIBRARY_READ_PATH,
                            condition_name="Library",
                            cas9_type="")
        if result is False:
            msg = "Error encountered during library trimming." \
                  "\nCould not open initial library read file!" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # update progress bar
        self.update_progressbar()

        # all Cas9 positive samples
        result = self.do_trimming_v2(file_list=global_vars.LIST_INITIAL_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[0])
        if result is False:
            msg = "Error encountered during initial Cas9-positive sample trimming." \
                  "\nCould not open initial Cas9-positive sample(s)" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        result = self.do_trimming_v2(file_list=global_vars.LIST_FINAL_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[0])
        if result is False:
            msg = "Error encountered during final Cas9-positive sample trimming." \
                  "\nCould not open final Cas9-positive sample(s)" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # all Cas9 negative samples
        result = self.do_trimming_v2(file_list=global_vars.LIST_CAS9_NEG_I_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[1])
        if result is False:
            msg = "Error encountered during initial Cas9-negative sample trimming." \
                  "\nCould not open initial Cas9-negative sample(s)" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        result = self.do_trimming_v2(file_list=global_vars.LIST_CAS9_NEG_F_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[1])
        if result is False:
            msg = "Error encountered during final Cas9-negative sample trimming." \
                  "\nCould not open final Cas9-negative sample(s)" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # build index
        result = self.build_index()
        if result is False:
            msg = "Error encountered during index assembly. " \
                  "\nCould not create index files!" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # update progress bar
        self.update_progressbar()

        # align reads
        write_to_log("Aligning reads with library index...")

        # call do_bowtie2_alignment function to align the library reads using bowtie2
        # Library reads
        result = self.do_bowtie2_alignment(file_list=global_vars.LIBRARY_READ_PATH,
                                  condition_name="Library",
                                  cas9_type="")
        if result is False:
            msg = "Error encountered during library alignment." \
                  "\nCould not create alignment file!" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # update progress bar
        self.update_progressbar()

        # run mageck_count function to get raw read counts
        result = self.do_mageck_count(file_list_1=global_vars.LIST_INITIAL_PATHS,
                             file_list_2=global_vars.LIST_FINAL_PATHS,
                             condition_name_1=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                             condition_name_2=global_vars.EXPERIMENT_SETTINGS['Final condition name'])
        if result is False:
            msg = "Error encountered during read count calculations." \
                  "\nCould not create read counts file!" \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # update progress bar
        self.update_progressbar()

        write_to_log("Done alignments and read count generation!")
        write_to_log("Starting core TRACS calculations...")

        # run function to analyze reads with TRACS algorithm
        result = self.tracs_core()
        if result is False:
            msg = "Error encountered during core TRACS calculations." \
                  "\nStopping TRACS."
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return

        # self.controller.show_frame("EndAnalysis")
        write_to_log("All done!")

        # enable visualize button
        global_vars.btnVisualize.config(state="normal")

    # function to call cutadapt to trim provided read sample files
    def do_trimming_v2(self, file_list, condition_name, cas9_type):

        # change directory
        trimmed_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                           global_vars.FILE_FLAGS['Trimmed dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(trimmed_target_path) is False:
            try:
                write_to_log("Creating directory for trimmed files...")
                os.mkdir(global_vars.FILE_FLAGS['Trimmed dir'])
            except IOError:
                write_to_log("Could not create directory for trimmed files! "
                             "Please check the permissions for your experiment folder!")
                return False

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

            # check if input file is accessible
            result = check_file_access(readfile)

            if result is True:

                # determine line count in file to overcome multicore ban in cutadapt
                line_count = (os.popen('wc -l \'%s\'' % readfile).read()).split(" ", 1)

                # build command pieces
                cmd_prefix = "head \'%s\' -n %s| " % (readfile, line_count[0])
                cmd_main = "cutadapt -g NAAGGACGAAACANNN...GTTTTAGAGCTAGAAAT -> "
                cmd_f_out = "\'" + this_filename + "\' "
                cmd_args = "-e 0.1 --cores %i --trimmed-only" % global_vars.EXPERIMENT_SETTINGS['CPU cores']

                cmd = cmd_prefix + cmd_main + cmd_f_out + cmd_args

                write_to_log("Trimming file %s and saving as %s" % (readfile, this_filename))

                # run command in shell
                cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
                show_in_console = global_vars.show_in_console

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

                # update progress bar - not for library but for Cas9 positive and negative samples only
                # the one only updates for the initial and final condition samples
                self.update_progressbar()
            else:
                write_to_log("Could not trim file: %s!" % readfile)
                return False

        return True

    # function to convert CSV reference library file to fasta format for use with bowtie2
    def convert_to_fasta(self):

        result = check_file_access(global_vars.EXPERIMENT_SETTINGS['Library reference file'])

        if result is True:

            # change directory
            os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

            # build command pieces
            cmd_prefix = "awk -F \',\' \'{print \">\"$1\"\\n\"$2}\' "
            cmd_f_in = "\'" + global_vars.EXPERIMENT_SETTINGS['Library reference file'] + "\'"
            cmd_f_out = " > \'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + \
                        global_vars.FILE_FLAGS['Fasta library file'] + "\'"
            cmd = cmd_prefix + cmd_f_in + cmd_f_out

            write_to_log("Converting reference sgRNA library file from CSV to FASTA...")

            # run command
            cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            show_in_console = global_vars.show_in_console

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

            write_to_log("Finished converting reference sgRNA library file!")
            return True
        else:
            write_to_log("Failed to convert reference sgRNA library file to FASTA!")
            return False

    # function to generate bowtie2 index for reference library
    def build_index(self):

        # change directory
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        # check if library FASTA file is accessible
        result = check_file_access(global_vars.EXPERIMENT_SETTINGS['Experiment name'] +
                                   global_vars.FILE_FLAGS['Fasta library file'])

        if result is True:

            write_to_log("Building sgRNA library index...")

            global_vars.FILE_FLAGS['Bowtie2 index full path'] = \
                os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                             global_vars.FILE_FLAGS['Bowtie2 index name'] + "/")

            if os.path.exists(global_vars.FILE_FLAGS['Bowtie2 index full path']) is False:

                try:
                    write_to_log("Creating bowtie2-index directory...")
                    os.mkdir(global_vars.FILE_FLAGS['Bowtie2 index name'])
                except IOError:
                    write_to_log("Could not create bowtie2-index directory! "
                                 "Please check the permissions for your experiment folder!")
                    return False

            os.chdir(global_vars.FILE_FLAGS['Bowtie2 index name'])

            # build command pieces
            # structure: bowtie2-build LibraryAB_fasta.fa bowtie2_index_LibraryAB
            cmd_prefix = "bowtie2-build "
            cmd_f_in = "\'../" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + \
                       global_vars.FILE_FLAGS['Fasta library file'] + "\' "
            cmd_f_out = "\'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "-" + \
                        global_vars.FILE_FLAGS['Bowtie2 index name'] + "\'"
            cmd = cmd_prefix + cmd_f_in + cmd_f_out

            # run command
            cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            show_in_console = global_vars.show_in_console

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

            write_to_log("Finished building sgRNA index!")
            return True
        else:
            write_to_log("Could not open Library FASTA file!")
            return False

    # function to perform alignments for provided sample read files using bowtie2
    def do_bowtie2_alignment(self, file_list, condition_name, cas9_type):

        # path to alignments target directory
        alignment_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                             global_vars.FILE_FLAGS['Alignments dir'] + "/")

        # make alignment target directory if it doesn't exist
        if os.path.exists(alignment_target_path) is False:
            write_to_log("Creating directory for aligned files...")
            os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
            os.mkdir(global_vars.FILE_FLAGS['Alignments dir'])

        trimmed_read_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                         global_vars.FILE_FLAGS['Trimmed dir'])

        # change directory
        os.chdir(alignment_target_path)

        # loop through file list for current condition (Library, Initial, Final conditions)
        for i in file_list:

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

            # get index full path + index name
            index = global_vars.FILE_FLAGS['Bowtie2 index full path'] + \
                         global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "-" + \
                         global_vars.FILE_FLAGS['Bowtie2 index name']

            # check if each index file is accessible; six files in total
            # result_index = check_file_access([index+".1.bt2", index+".1.bt2")

            # check if trimmed read file is accessible
            result_input = check_file_access(input_file)

            if result_input is True:

                write_to_log("Getting alignment data for file: %s" % input_file)

                # build command pieces
                cmd_prefix = "bowtie2 -x \'%s\'" % index
                cmd_f_in = " -U \'%s\'" % input_file
                cmd_suffix = " --norc -N 1 -p %i | samtools view -bS - > \'%s\'" % (4 * int(global_vars.EXPERIMENT_SETTINGS['CPU cores']), output_file)
                cmd = cmd_prefix + cmd_f_in + cmd_suffix

                # run command
                cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
                show_in_console = global_vars.show_in_console

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
                    write_to_log("Successfully aligned initial library reads!")
                    break
                else:
                    write_to_log("Successfully aligned file: %s!" % input_file)
            else:
                write_to_log("Error opening trimmed read file: %s!" % input_file)
                return False
        return True

    # function to run mageck's count function to get read counts from alignments for all samples
    def do_mageck_count(self, file_list_1, file_list_2, condition_name_1, condition_name_2):

        write_to_log("Calculating read counts...")

        # change directory
        counts_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                             global_vars.FILE_FLAGS['Readcounts dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(counts_target_path) is False:
            try:
                write_to_log("Creating directory for read count files...")
                os.mkdir(global_vars.FILE_FLAGS['Readcounts dir'])
            except IOError:
                write_to_log("Could not create directory for read count files! "
                             "Please check the permissions for your experiment folder!")
                return False

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

        # check if trimmed library read file is accessible
        result = check_file_access(input_file)

        if result is True:
            # loop for all Cas9+ and Cas9- sample read files
            for j in global_vars.CAS9_TYPE_NAMING:

                # initial condition labels and files
                for i in file_list_1:

                    # get replicate name
                    next_char = self.controller.get_rep_char(file_list_1.index(i))

                    input_file = os.path.join(aligned_read_path, condition_name_1 + j + next_char + global_vars.FILE_FLAGS['Trimmed read file'])

                    write_to_log("Processing file: %s" % input_file)

                    cmd_sample_list = cmd_sample_list + condition_name_1 + j + next_char + ","

                    cmd_input_file_list = cmd_input_file_list + "\'" + input_file + "\' "

                    if check_file_access(input_file) is False:
                        write_to_log("Could not access trimmed read file: %s" % input_file)
                        return False

                # final condition labels and files
                for i in file_list_2:

                    # get replicate name
                    next_char = self.controller.get_rep_char(file_list_2.index(i))

                    input_file = os.path.join(aligned_read_path, condition_name_2 + j + next_char + global_vars.FILE_FLAGS['Trimmed read file'])

                    write_to_log("Processing file: %s" % input_file)

                    cmd_sample_list = cmd_sample_list + condition_name_2 + j + next_char + ","

                    cmd_input_file_list = cmd_input_file_list + "\'" + input_file + "\' "

                    if check_file_access(input_file) is False:
                        write_to_log("Could not access trimmed read file: %s" % input_file)
                        return False

            # build command pieces
            cmd_prefix = "mageck count -l \'" + global_vars.EXPERIMENT_SETTINGS['Library reference file'] + "\'"
            cmd_namtag = " -n \'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "\'"
            cmd_sampletags = " --sample-label \'" + cmd_sample_list[:-1] + "\'"
            cmd_readfiles = " --fastq " + cmd_input_file_list
            cmd = cmd_prefix + cmd_namtag + cmd_sampletags + cmd_readfiles

            # run command
            cmd_output = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
            show_in_console = global_vars.show_in_console

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
        else:
            write_to_log("Could not access trimmed library FASTQ file!")
            return False

        write_to_log("Finished calculating read counts!")
        return True

    # function to run TRACS algorithm to calculate Enrichment Scores and Enrichment Ratio from read counts table
    def tracs_core(self):

        # change directory
        reads_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                         global_vars.FILE_FLAGS['Readcounts dir'] + "/")

        os.chdir(reads_target_path)

        # get name of read counts file
        readsfile = global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".count.txt"

        result = check_file_access(readsfile)

        if result is False:
            write_to_log("Could not open read counts file!")
            return False

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
        write_to_log("Number of replicates: %i" % self.num_replicates)

        # # # NORMALIZING READS # # #
        write_to_log("Normalizing read counts...")

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

        # update progress bar
        self.update_progressbar()

        # # # CALCULATE LOG2 FOLD CHANGE # # #
        write_to_log("Calculating log2 fold changes...")

        # save original column names
        orig_column_names = self.raw_counts.columns

        # log2FC structure: log2FC = [final reads / initial reads]
        # call function to calculate log2FC for library-Cas9-neg-initial [initial Cas9 neg / library]
        try:
            self.lib_log2FC_calculator()
        except KeyError:
            msg = "Error encountered during log2 fold change calculations." \
                  "\nDid you enter the same file for the initial library reads and sample reads?"
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return False

        # call function to calculate log2FC for initial condition [initial Cas9 pos / initial Cas9 neg ]
        self.log2FC_calculator(name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'])

        # call function to calculate log2FC for final condition [final Cas9 pos / final Cas9 neg ]
        self.log2FC_calculator(name = global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        # update progress bar
        self.update_progressbar()

        # # # GENE SCORE # # #
        write_to_log("Calculating gene scores...")

        # determine number of sgRNAs for each gene
        self.num_sgRNAs = pd.DataFrame()
        self.num_sgRNAs['num sgRNA'] = self.raw_counts.groupby(self.raw_counts.index)[orig_column_names[0]].count()

        # update progress bar
        self.update_progressbar()

        write_to_log("Determining average...")

        # call function to calculate average of log2FC for every gene for each log2FC sample column
        # calculates a single log2FC from all the guides for each gene
        self.avg = pd.DataFrame()  # new dataframe to hold averages
        self.avg_calculator(name = global_vars.COLUMN_NAMES['Library'])
        self.avg_calculator(name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'])
        self.avg_calculator(name = global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        # update progress bar
        self.update_progressbar()

        # # # ENRICHMENT SCORE # # #
        write_to_log("Calculating sgw...")

        # determine sgw for each gene
        self.sgw = pd.DataFrame()  # new dataframe to hold sgw
        self.get_sgw()

        # update progress bar
        self.update_progressbar()

        write_to_log("Calculating rank...")
        # determine rank for each sample
        self.get_rank()

        # update progress bar
        self.update_progressbar()

        # save to file
        # write_to_log("Saving to weighted ranks (sgw) to file ...")
        # os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        # self.sgw.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".sgw.csv", index=True,
        #                 sep=',',
        #                 float_format='%.16f')

        write_to_log("Calculating enrichment scores...")
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

        # update progress bar
        self.update_progressbar()

        # save to file
        # write_to_log("Saving enrichment scores to file ...")
        # os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        # self.ES.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".ES.csv", index=True,
        #                        sep=',',
        #                        float_format='%.16f')

        # # # ENRICHMENT RATIO # # #
        write_to_log("Calculating enrichment ratio...")

        # determine enrichment ratio (ER); log2 ratio of ES ( log2(final ES / initial ES)
        initial_col_name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'] + global_vars.COLUMN_NAMES['ES']
        final_col_name = global_vars.EXPERIMENT_SETTINGS['Final condition name'] + global_vars.COLUMN_NAMES['ES']
        self.ES[global_vars.COLUMN_NAMES['ER']] = np.log2(self.ES[final_col_name].astype(float) / self.ES[initial_col_name].astype(float))
        self.ES[global_vars.COLUMN_NAMES['ER']] = self.ES[global_vars.COLUMN_NAMES['ER']].map(lambda x: '%.6f' % x)

        # update progress bar
        self.update_progressbar()

        # # # CALCULATE STATISTICS # # #
        write_to_log("Calculating statistics...")

        # create temp_stats dataframe
        self.temp_stats = pd.DataFrame()
        # call function to calculate statistics
        self.get_stats(initial_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                       final_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        # copy stats to end of ES dataframe
        self.ES['pval'] = self.temp_stats['pvalue']
        self.ES['qval'] = self.temp_stats['qvalue']

        # update progress bar
        self.update_progressbar()

        # rename columns in final ES dataframe
        self.ES = self.ES.rename(columns={initial_col_name: global_vars.COLUMN_NAMES['Initial'], final_col_name: global_vars.COLUMN_NAMES['Final']})
        
        # save to file
        try:
            write_to_log("Saving results to file ...")
            os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
            self.ES.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".csv",
                           index=True, sep=',', float_format='%.6f')
        except IOError:
            msg = "Could not save results to file! Please check the permissions for your experiment folder!"
            write_to_log(msg)
            messagebox.showerror(title="Error",
                                 message=msg)
            return False

        # update progress bar
        self.update_progressbar()

        return True

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

        # get column names from self.avg dataframecontroller.show_frame("StartAnalysis")
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


class EndAnalysis(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = tk.Label(fr_header_1, text="Analysis finished", font=controller.title_font, bg="#3d6fdb", fg="white")
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_exp_settings_Cas9_trans_800x175.png")
        lbl_exp_design_image = tk.Label(fr_banner, image=img_banner, bg="#3d6fdb")
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # label
        fr_1 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = tk.Label(fr_1, text="Analysis complete",
                         wraplength="900", justify="left", font=controller.intro_body_font, bg="#3d6fdb", fg="white")
        lbl_info.pack(fill=X, side=LEFT, padx=5, pady=5)

        # summary text box
        fr_2 = tk.Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)

        # navigation buttons
        btn_Next = Button(self, text="Exit TRACS", command=lambda: app.quit())
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("StartAnalysis"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

        btn_Back.pack(side=RIGHT, padx=5, pady=5)


if __name__ == "__main__":
    app = TRACSApp()
    app.geometry("940x680+500+50")
    app.resizable(False, False)
    _thread.start_new_thread(app.mainloop(), ())



