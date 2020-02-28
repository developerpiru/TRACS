import tkinter as tk
from tkinter import font as tkfont
from tkinter import *
from tkinter import filedialog
#from tkinter import ttk
from tkinter.ttk import Frame, Label, Button, Style
import os
import subprocess, sys
import pandas as pd
import numpy as np
import global_vars_v1 as global_vars
from scipy import stats
import webbrowser

# Fixed lib log2FC calcs by addressing columns by NAME
# Fixed addressing columns by name for reg log2FC
# validated with Excel TRACS
# dynamic column naming and addressing working for all functions in read count processing
# implemented fixed replicate naming of -A, -B, -C, etc for replicates
# removed replicate naming stuff from window frame

# # # DONE # # #
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

# # # TO DO # # #
# output R file for web-based viewing tool
# move summary_box outputs to separate function call ???? --- stuck in infinite loop with bowtie2 function
# implement error catching

# version
app_version = "1.1.1"

# please go to the official TRACS GitHub page for documentation and the latest updates
github_url = "https://github.com/developerpiru/TRACS"
# if you experience bugs/errors/other hurdles in your analyses, please leave a bug report on GitHub

# Data visualization and exploration for TRACS is enabled by VisualizeTRACS, a companion R shiny app
# You can launch VisualizeTRACS from its GitHub page
viz_url = "https://github.com/developerpiru/VisualizeTRACS"
# Help and documentation for VisualizeTRACS is available there


class TRACSApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # load global variables
        global_vars.init_vars()

        # widget styles
        self.style = Style()
        # self.style.theme_use("default")
        
        self.frame_relief = FLAT
        self.frame_borderwidth = 0
        
        # app title
        self.winfo_toplevel().title("TRACS v" + app_version)

        # set fonts
        self.title_font = tkfont.Font(family='Arial', size=18, weight="bold")
        self.label_font = tkfont.Font(family='Arial', size=12, weight="bold")
        self.label_font_small = tkfont.Font(family='Arial', size=8, weight="bold")
        self.body_font = tkfont.Font(family='Arial', size=12, weight="normal")

        # main container
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # define page frames
        self.frames = {}
        # create frames for each page
        for F in (MainPage, NewExpSettings, LoadFiles, LoadCas9NegFiles_v2, SummaryPage, TrimmingReads, FileExtractor):
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
    def browse_folder(self, frame_name, tar_name):
        target_obj = getattr(self.frames[frame_name], tar_name)
        folder = filedialog.askdirectory(initialdir=os.getcwd(), title='Select folder...')
        target_obj.delete(0, tk.END)
        target_obj.insert(tk.END, folder)

    # function to select files
    def browse_file(self, multi_files, frame_name, tar_name):
        target_obj = getattr(self.frames[frame_name], tar_name)
        if multi_files == "TRUE":
            files = filedialog.askopenfilenames(title='Select files...',
                                                filetypes=(("FASTQ files", "*.fastq"),
                                                           ("GZipped Fastq", "*fastq.gz"), ("All files", "*")))
            for item in files:
                target_obj.insert(END, item)
        else:
            file = filedialog.askopenfilename(title='Select file...',
                                              filetypes=(("FASTQ files", "*.fastq"),
                                                         ("GZipped Fastq", "*fastq.gz"), ("CSV", "*.csv"),
                                                         ("All files", "*")))
            target_obj.delete(0, END)
            target_obj.insert(END, file)

    # function to delete selected file from file list
    def delete_selected(self, frame_name, tar_name):
        target_obj = getattr(self.frames[frame_name], tar_name)
        selected = target_obj.curselection()
        for i in selected[::-1]:
            target_obj.delete(i)

    # function to move up file in list ---- not required anymore
    def move_up(self, frame_name, tar_name):
        target_obj = getattr(self.frames[frame_name], tar_name)
        selected = target_obj.curselection()

        if not selected:
            return
        for i in selected:
            if i != 0:
                saved_text = target_obj.get(i)
                target_obj.delete(i)
                target_obj.insert(i-1, saved_text)
                target_obj.select_set(i-1)

    # function to move down file in list ---- not required anymore
    def move_down(self, frame_name, tar_name):
        target_obj = getattr(self.frames[frame_name], tar_name)
        selected = target_obj.curselection()

        if not selected:
            return
        for i in selected:
            if i != target_obj.size()-1:
                saved_text = target_obj.get(i)
                target_obj.delete(i)
                target_obj.insert(i+1, saved_text)
                target_obj.select_set(i+1)

    # function to get next replicate alphabetical character (ie. A, B, C, etc.)
    # requires first replicate letter ('A') and incremental value (n)
    # returns correct replicate letter for the A+nth replicate
    def get_rep_char(self, n):
        start_char = "A"
        next_char = chr(ord(start_char) + int(n)).upper()
        return str(next_char)

    def extract_gzip(self):
        FILE_PATH_LIST = self.frames['FileExtractor'].file_list.get(0, END)
        summary_box = self.frames['FileExtractor'].txt_summary
        cmd = "gunzip -k"
        for i in FILE_PATH_LIST:
            cmd = cmd + " %s" % i

        summary_box.insert(END, "\n>" + cmd)
        cmd_output = os.popen(cmd).read()
        summary_box.insert(END, cmd_output)

    def new_window_browse(self, file_list):
        target_obj = getattr(self.new_window, file_list)
        files = filedialog.askopenfilenames(title='Select files...', filetypes=(("GZipped Fastq", "*fastq.gz")))
        for item in files:
            target_obj.insert(END, item)


class MainPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        frame = Frame(self, relief=GROOVE, borderwidth=0)
        frame.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # info
        fr_1 = Frame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = Label(fr_1, text="TRACS: Toolset for the Ranked Analysis of CRISPR Screens",
                         wraplength="800", justify="left", font=controller.title_font)
        lbl_info.pack(padx=5, pady=5)

        fr_2 = Frame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=X, side=TOP, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = Label(fr_2, text="For help and documentation, please visit our GitHub page at:"
                                    "\n" + github_url +
                                    " or click the button below to open it in your browser."
                                    "\n\nData visualization and exploration of data generated by TRACS is enabled by "
                                    "its companion R shiny app, VisualizeTRACS."
                                    "\nYou can access VisualizeTRACS at its "
                                    "own GitHub page at " + viz_url +
                                    "\n\nYou can find help and documentation for VisualizeTRACS there as well."
                                    "\n\nPlease report all bugs for TRACS and VisualizeTRACS on their GitHub pages!",
                         wraplength="900", justify="left", font=controller.body_font)
        lbl_info.pack(side=TOP, padx=5, pady=5)
        lbl_version = Label(fr_2, text="\n\n This is version " + app_version,
                            wraplength="800", justify="left", font=controller.body_font)
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
        fr_main_1 = Frame(self, relief=GROOVE, borderwidth=0)
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = Label(fr_header_1, text="Step 1: New experiment settings", font=controller.title_font)
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_exp_settings_Cas9_trans_800x175.png")
        lbl_exp_design_image = Label(fr_banner, image=img_banner)
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # experiment name
        fr_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_exp_name = Label(fr_1, text="Experiment name", width="40", font=controller.label_font)
        lbl_exp_name.pack(side=LEFT, padx=5, pady=5)
        self.txt_exp_name = Entry(fr_1)
        self.txt_exp_name.pack(fill=X, padx=5, pady=5, expand=True)

        # output directory
        fr_2 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_expdir = Label(fr_2, text="Experiment directory", width="40", font=controller.label_font)
        lbl_expdir.pack(side=LEFT, anchor=N, padx=5, pady=5)
        self.txt_expdir = Entry(fr_2, width="50")
        self.txt_expdir.pack(fill=X, padx=5, pady=5, expand=True)
        btn_browse = Button(fr_2, text="Browse", command=lambda: self.controller.browse_folder(
            frame_name="NewExpSettings", tar_name="txt_expdir"))
        btn_browse.pack(side=RIGHT, padx=5, pady=5)

        # library type
        fr_3 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_3.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)

        # name of initial group
        fr_4 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_4.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_name_initial = Label(fr_4, text="Name of initial condition (T0)", width="40", font=controller.label_font)
        lbl_name_initial.pack(side=LEFT, padx=5, pady=5)
        self.txt_name_initial = Entry(fr_4)
        self.txt_name_initial.pack(fill=X, padx=5, pady=5, expand=True)

        # name of final group
        fr_5 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_5.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_name_final = Label(fr_5, text="Name of final condition (Tf)", width="40", font=controller.label_font)
        lbl_name_final.pack(side=LEFT, padx=5, pady=5)
        self.txt_name_final = Entry(fr_5)
        self.txt_name_final.pack(fill=X, padx=5, pady=5, expand=True)

        # false discovery rate (FDR) for stats
        fr_6 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_6.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_fdr = Label(fr_6, text="False discovery rate (FDR, default is 0.05)", width="40", font=controller.label_font)
        lbl_fdr.pack(side=LEFT, padx=5, pady=5)
        self.txt_fdr = Entry(fr_6)
        self.txt_fdr.pack(fill=X, padx=5, pady=5, expand=True)

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
        fr_main_1 = Frame(self, relief=GROOVE, borderwidth=0)
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = Label(fr_header_1, text="Step 2: Load Cas9-positive files", font=controller.title_font)
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_load_files_trans_800x175.png")
        lbl_exp_design_image = Label(fr_banner, image=img_banner)
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # reference library
        fr_2 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_reflibrary = Label(fr_2, text="Library reference file (CSV)", width="40", font=controller.label_font)
        lbl_reflibrary.pack(side=LEFT, anchor=N, padx=5, pady=5)
        self.txt_reflibrary = tk.Entry(fr_2)
        self.txt_reflibrary.pack(fill=X, padx=5, pady=5, expand=True)
        btn_browse = Button(fr_2, text="Browse", command=lambda: self.controller.browse_file(multi_files="FALSE",
                                                                                             frame_name="LoadFiles",
                                                                                             tar_name="txt_reflibrary"))
        btn_browse.pack(side=RIGHT, padx=5, pady=5)

        # library reads file
        fr_6 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_6.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_libreads = Label(fr_6, text="Initial library (L0) read file (FASTQ)", width="40", font=controller.label_font)
        lbl_libreads.pack(side=LEFT, anchor=N, padx=5, pady=5)
        self.txt_libreads = tk.Entry(fr_6)
        self.txt_libreads.pack(fill=X, padx=5, pady=5, expand=True)
        btn_browse = Button(fr_6, text="Browse", command=lambda: self.controller.browse_file(multi_files="FALSE",
                                                                                             frame_name="LoadFiles",
                                                                                             tar_name="txt_libreads"))
        btn_browse.pack(side=RIGHT, padx=5, pady=5)

        # container frame for listbox labels
        fr_3 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_3.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial label
        fr_4 = Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_4.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_i_info = Label(fr_4, text="Initial condition (T0) Cas9-positive read files", font=controller.label_font)
        lbl_i_info.pack(padx=5, pady=5, side=LEFT)

        # frame for final label
        fr_5 = Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_5.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_f_info = Label(fr_5, text="Final condition (Tf) Cas9-positive read files", font=controller.label_font)
        lbl_f_info.pack(padx=5, pady=5, side=LEFT)

        # container frame for list boxes
        list_container = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        list_container.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial listbox
        initial_list_container = Frame(list_container, relief=controller.frame_relief,
                                       borderwidth=controller.frame_borderwidth)
        initial_list_container.pack(fill=BOTH, padx=0, pady=0, side=LEFT, expand=True)
        # frame for initial listbox buttons
        initial_btn_container = Frame(initial_list_container)
        initial_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # initial listbox
        self.lst_initial_files = Listbox(initial_list_container, height="8", width="45", selectmode=EXTENDED)

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
                                                                                                          frame_name="LoadFiles",
                                                                                                          tar_name="lst_initial_files"))
        btn_i_add.pack(side=TOP, padx=0, pady=0)
        # remove
        btn_i_rem = Button(initial_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="LoadFiles",
                                                                           tar_name="lst_initial_files"))
        btn_i_rem.pack(side=TOP, padx=0, pady=0)
        # clear
        btn_i_clear = Button(initial_btn_container, text="Clear", command=lambda: self.lst_initial_files.delete(0, END))
        btn_i_clear.pack(side=TOP, padx=0, pady=0)

        # frame for initial listbox
        final_list_container = Frame(list_container, relief=controller.frame_relief,
                                     borderwidth=controller.frame_borderwidth)
        final_list_container.pack(fill=BOTH, padx=0, pady=0, side=LEFT, expand=True)
        # frame for initial listbox buttons
        final_btn_container = Frame(final_list_container)
        final_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # final listbox
        self.lst_final_files = Listbox(final_list_container, height="8", width="45", selectmode=EXTENDED)

        # y scroll for initial listbox
        yscrl_lst_final_files = Scrollbar(final_list_container, orient=VERTICAL)
        yscrl_lst_final_files.pack(side=RIGHT, fill=Y)
        yscrl_lst_final_files.config(command=self.lst_final_files.yview)
        self.lst_final_files.config(yscrollcommand=yscrl_lst_final_files.set)
        # x scroll for initial listbox
        xscrl_lst_final_files = Scrollbar(final_list_container, orient=HORIZONTAL)
        xscrl_lst_final_files.pack(side=BOTTOM, fill=X)
        xscrl_lst_final_files.config(command=self.lst_final_files.xview)
        self.lst_final_files.config(xscrollcommand=xscrl_lst_final_files.set)

        self.lst_final_files.pack(side=LEFT, padx=0, pady=0)

        # add
        btn_f_add = Button(final_btn_container, text="Add",
                           command=lambda: self.controller.browse_file(multi_files="TRUE",
                                                                       frame_name="LoadFiles",
                                                                       tar_name="lst_final_files"))
        btn_f_add.pack(side=TOP, padx=0, pady=0)
        # remove
        btn_f_rem = Button(final_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="LoadFiles",
                                                                           tar_name="lst_final_files"))
        btn_f_rem.pack(side=TOP, padx=0, pady=0)
        # clear
        btn_f_clear = Button(final_btn_container, text="Clear", command=lambda: self.lst_final_files.delete(0, END))
        btn_f_clear.pack(side=TOP, padx=0, pady=0)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("LoadCas9NegFiles_v2"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("NewExpSettings"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


# function to load Cas9- files
class LoadCas9NegFiles_v2(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = Frame(self, relief=GROOVE, borderwidth=0)
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = Label(fr_header_1, text="Step 3: Load Cas9-negative files", font=controller.title_font)
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_load_files_trans_800x175.png")
        lbl_exp_design_image = Label(fr_banner, image=img_banner)
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # info
        fr_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = Label(fr_1, text="Select Cas9-negative sequencing read files. These are the read files from your"
                                    " cells that do not express Cas9 and only contain the sgRNA pool! ",
                         wraplength="900", justify="left", font=controller.body_font)
        lbl_info.pack(side=LEFT, padx=5, pady=5)

        # container frame for listbox labels
        fr_3 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_3.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial label
        fr_4 = Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_4.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_i_info = Label(fr_4, text="Initial condition (T0) Cas9-negative read files", font=controller.label_font)
        lbl_i_info.pack(padx=5, pady=5, side=LEFT)

        # frame for final label
        fr_5 = Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_5.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_f_info = Label(fr_5, text="Final condition (T0) Cas9-negative read files", font=controller.label_font)
        lbl_f_info.pack(padx=5, pady=5, side=LEFT)

        # container frame for list boxes
        list_container = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        list_container.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial listbox
        initial_list_container = Frame(list_container, relief=controller.frame_relief,
                                       borderwidth=controller.frame_borderwidth)
        initial_list_container.pack(fill=BOTH, padx=5, pady=5, side=LEFT, expand=True)

        # frame for initial listbox buttons
        initial_btn_container = Frame(initial_list_container)
        initial_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # initial listbox
        self.lst_initial_files = Listbox(initial_list_container, height="12", width="45", selectmode=EXTENDED)

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
                                                                                                          frame_name="LoadCas9NegFiles_v2",
                                                                                                          tar_name="lst_initial_files"))
        btn_i_add.pack(side=TOP, padx=0, pady=0)
        # remove
        btn_i_rem = Button(initial_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="LoadCas9NegFiles_v2",
                                                                           tar_name="lst_initial_files"))
        btn_i_rem.pack(side=TOP, padx=0, pady=0)
        # clear
        btn_i_clear = Button(initial_btn_container, text="Clear", command=lambda: self.lst_initial_files.delete(0, END))
        btn_i_clear.pack(side=TOP, padx=0, pady=0)

        # frame for initial listbox
        final_list_container = Frame(list_container, relief=controller.frame_relief,
                                       borderwidth=controller.frame_borderwidth)
        final_list_container.pack(fill=BOTH, padx=5, pady=5, side=LEFT, expand=True)
        # frame for initial listbox buttons
        final_btn_container = Frame(final_list_container)
        final_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # initial listbox
        self.lst_final_files = Listbox(final_list_container, height="12", width="45", selectmode=EXTENDED)

        # y scroll for initial listbox
        yscrl_lst_final_files = Scrollbar(final_list_container, orient=VERTICAL)
        yscrl_lst_final_files.pack(side=RIGHT, fill=Y)
        yscrl_lst_final_files.config(command=self.lst_final_files.yview)
        self.lst_final_files.config(yscrollcommand=yscrl_lst_final_files.set)

        # x scroll for initial listbox
        xscrl_lst_final_files = Scrollbar(final_list_container, orient=HORIZONTAL)
        xscrl_lst_final_files.pack(side=BOTTOM, fill=X)
        xscrl_lst_final_files.config(command=self.lst_final_files.xview)
        self.lst_final_files.config(xscrollcommand=xscrl_lst_final_files.set)

        self.lst_final_files.pack(side=LEFT, padx=0, pady=0)

        # add
        btn_f_add = Button(final_btn_container, text="Add",
                           command=lambda: self.controller.browse_file(multi_files="TRUE",
                                                                       frame_name="LoadCas9NegFiles_v2",
                                                                       tar_name="lst_final_files"))
        btn_f_add.pack(side=TOP, padx=0, pady=0)
        # remove
        btn_f_rem = Button(final_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="LoadCas9NegFiles_v2",
                                                                           tar_name="lst_final_files"))
        btn_f_rem.pack(side=TOP, padx=0, pady=0)
        # clear
        btn_f_clear = Button(final_btn_container, text="Clear", command=lambda: self.lst_final_files.delete(0, END))
        btn_f_clear.pack(side=TOP, padx=0, pady=0)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: self.process_inputs())
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("LoadFiles"))
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
        global_vars.EXPERIMENT_SETTINGS['Library reference file'] = \
            self.controller.frames['LoadFiles'].txt_reflibrary.get()  # library reference file

        # get library read fastq file
        global_vars.LIBRARY_READ_PATH = self.controller.frames['LoadFiles'].txt_libreads.get()

        print("\nLibrary read path: %s" % global_vars.LIBRARY_READ_PATH)

        # get Initial condition Cas9+ fastq files
        global_vars.LIST_INITIAL_PATHS = self.controller.frames['LoadFiles'].lst_initial_files.get(0, END)

        # get Final condition Cas9+ fastq files
        global_vars.LIST_FINAL_PATHS = self.controller.frames['LoadFiles'].lst_final_files.get(0, END)

        # get Initial condition Cas9- fastq files
        global_vars.LIST_CAS9_NEG_I_PATHS = self.controller.frames['LoadCas9NegFiles_v2'].lst_initial_files.get(0, END)

        # get Final condition Cas9- fastq files
        global_vars.LIST_CAS9_NEG_F_PATHS = self.controller.frames['LoadCas9NegFiles_v2'].lst_final_files.get(0, END)

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
        fr_main_1 = Frame(self, relief=GROOVE, borderwidth=0)
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = Label(fr_header_1, text="Summary", font=controller.title_font)
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_exp_settings_Cas9_trans_800x175.png")
        lbl_exp_design_image = Label(fr_banner, image=img_banner)
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # experiment name
        fr_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = Label(fr_1, text="Please confirm all details below before continuing:",
                         wraplength="900", justify="left", font=controller.body_font)
        lbl_info.pack(fill=X, side=LEFT, padx=5, pady=5)

        # summary text box
        fr_2 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        self.txt_summary = Text(fr_2, height="18")
        self.txt_summary.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        scrl_txt_summary = Scrollbar(fr_2)
        scrl_txt_summary.pack(side=RIGHT, fill=Y)
        scrl_txt_summary.config(command=self.txt_summary.yview)
        self.txt_summary.config(yscrollcommand=scrl_txt_summary.set)

        # self.txt_summary.config(state=DISABLED)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("TrimmingReads"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("LoadCas9NegFiles_v2"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

        btn_Back.pack(side=RIGHT, padx=5, pady=5)


class FileExtractor(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        fr_main_1 = Frame(self, relief=GROOVE, borderwidth=0)
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # info
        fr_1 = Frame(fr_main_1)
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = Label(fr_1, text="Select files to extract",
                         wraplength="800", justify="left")
        lbl_info.pack(padx=5, pady=5)

        # container frame for listbox labels
        fr_2 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial listbox
        initial_list_container = Frame(fr_main_1, relief=controller.frame_relief,
                                       borderwidth=controller.frame_borderwidth)
        initial_list_container.pack(fill=BOTH, padx=5, pady=5, side=LEFT, expand=True)
        self.file_list = Listbox(initial_list_container, height="8", width="45", selectmode=EXTENDED)
        self.file_list.pack(side=LEFT, padx=5, pady=5)

        initial_btn_container = Frame(initial_list_container)
        initial_btn_container.pack(fill=BOTH, padx=5, pady=5, side=LEFT)

        # add
        btn_i_add = Button(initial_btn_container, text="Add",
                           command=lambda: self.controller.browse_file(multi_files="TRUE",
                                                                       frame_name="FileExtractor",
                                                                       tar_name="file_list"))
        btn_i_add.pack(side=TOP, padx=5, pady=5)
        # remove
        btn_i_rem = Button(initial_btn_container, text="Remove",
                           command=lambda: self.controller.delete_selected(frame_name="FileExtractor",
                                                                           tar_name="file_list"))
        btn_i_rem.pack(side=TOP, padx=5, pady=5)
        # clear
        btn_i_clear = Button(initial_btn_container, text="Clear", command=lambda: self.file_list.delete(0, END))
        btn_i_clear.pack(side=TOP, padx=5, pady=5)

        # summary text box
        fr_3 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_3.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        self.txt_summary = Text(fr_3, height="10")
        self.txt_summary.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        scrl_txt_summary = Scrollbar(fr_3)
        scrl_txt_summary.pack(side=RIGHT, fill=Y)
        scrl_txt_summary.config(command=self.txt_summary.yview)
        self.txt_summary.config(yscrollcommand=scrl_txt_summary.set)

        fr_4 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_4.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)

        btnExtract = Button(fr_4, text="Extract", command=lambda: self.controller.extract_gzip())
        btnExtract.pack()

        # navigation buttons
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("MainPage"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


class TrimmingReads(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = Frame(self, relief=GROOVE, borderwidth=0)
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = Label(fr_header_1, text="Step 4: Run analysis", font=controller.title_font)
        lbl_heading.pack(padx=5, pady=5, expand=True)

        # page banner
        fr_banner = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_banner.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        img_banner = PhotoImage(file="images/banners/banner_exp_settings_Cas9_trans_800x175.png")
        lbl_exp_design_image = Label(fr_banner, image=img_banner)
        lbl_exp_design_image.image = img_banner
        lbl_exp_design_image.pack(padx=5, pady=5, expand=True)

        # experiment name
        fr_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_1.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = Label(fr_1, text="You are now ready to start analysis! You have two options:"
                                    "\n\n1. If you are starting from raw read files to get read counts, select"
                                    " option 1. TRACS will automatically do the following:"
                                    "\n\n-Trim your reads\n-Build a reference library\n-Perform alignments"
                                    "\n-Generate a read counts\n-Run the TRACS algorithm to determine"
                                    " gene essentiality."
                                    "\n\n2. If you have already generated a read counts file in the correct format"
                                    " and have it saved in the 'readcounts' folder in your experiment directory, select"
                                    " option 2. Using your read count file, TRACS will automatically:"
                                    "\n\n-Run the TRACS algorithm to determine gene essentiality.",
                         wraplength="900", justify="left", font=controller.body_font)
        lbl_info.pack(fill=X, side=LEFT, padx=5, pady=5)

        # summary text box
        fr_2 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        #self.txt_summary = Text(fr_2, height="10")
        #self.txt_summary.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        #scrl_txt_summary = Scrollbar(fr_2)
        #scrl_txt_summary.pack(side=RIGHT, fill=Y)
        #scrl_txt_summary.config(command=self.txt_summary.yview)
        #self.txt_summary.config(yscrollcommand=scrl_txt_summary.set)

        # btnTrim = Button(fr_2, text="Trim Files", command=lambda: self.do_trimming_v3(file_list=global_vars.LIST_INITIAL_PATHS))

        # option 1 button: run the entire pipline
        btnTrim = Button(fr_2, text="1. I have raw reads (FASTQ files). Do everything for me!",
                         command=lambda: self.do_steps_option_1())
        btnTrim.pack(side=TOP)

        # option 2 button: only run the TRACS algorithm
        btnTest = Button(fr_2, text="2. I already have a read counts file (TXT). Just give me the essential genes!",
                         command=lambda: self.do_steps_option_2())
        btnTest.pack(side=TOP)

        # navigation buttons
        #btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("AnalyzeCounts"))
        #btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("SummaryPage"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    def do_steps_option_1(self):

        # call function to convert to fasta library
        self.convert_to_fasta()

        # trim reads
        print("\nTrimming reads...")
        #self.txt_summary.insert(END, "\nTrimming reads...")
        #self.txt_summary.update_idletasks()

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
        #self.txt_summary.insert(END, "\nAligning reads with library index...")
        #self.txt_summary.update_idletasks()

        # call do_bowtie2_alignment function to align reads for all samples
        # Library reads
        self.do_bowtie2_alignment(file_list=global_vars.LIBRARY_READ_PATH,
                                  condition_name="Library",
                                  cas9_type="")

        # all Cas9+ samples
        self.do_bowtie2_alignment(file_list=global_vars.LIST_INITIAL_PATHS,
                                  condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                                  cas9_type=global_vars.CAS9_TYPE_NAMING[0])
        self.do_bowtie2_alignment(file_list=global_vars.LIST_FINAL_PATHS,
                                  condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                                  cas9_type=global_vars.CAS9_TYPE_NAMING[0])

        # all Cas9- samples
        self.do_bowtie2_alignment(file_list=global_vars.LIST_CAS9_NEG_I_PATHS,
                                  condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                                  cas9_type=global_vars.CAS9_TYPE_NAMING[1])
        self.do_bowtie2_alignment(file_list=global_vars.LIST_CAS9_NEG_F_PATHS,
                                  condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                                  cas9_type=global_vars.CAS9_TYPE_NAMING[1])

        # run mageck_count function to get raw read counts
        self.do_mageck_count(file_list_1=global_vars.LIST_INITIAL_PATHS,
                             file_list_2=global_vars.LIST_FINAL_PATHS,
                             condition_name_1=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                             condition_name_2=global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        print("\n\nDone alignments and read count generation!")
        print("\n\nStarting TRACS algorithm...")

        # start do_steps_option_2 function
        self.do_steps_option_2()

    # function to force status box to scroll
    def scroll_down(self):
        self.txt_summary.see("end")
        self.after(1000, self.scroll_down)

    # function to call cutadapt to trim provided read sample files
    def do_trimming_v2(self, file_list, condition_name, cas9_type):

        # update status box
        #summary_box = self.txt_summary
        #summary_box.config(state=NORMAL)

        # change directory
        trimmed_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                           global_vars.FILE_FLAGS['Trimmed dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(trimmed_target_path) is False:
            print("\nCreating directory for trimmed files...")
            print("\n%s" % trimmed_target_path)
            #summary_box.insert(END, "\nCreating directory for trimmed files...")
            os.mkdir(global_vars.FILE_FLAGS['Trimmed dir'], 777)

        os.chdir(trimmed_target_path)

        for i in file_list:
            # update status box
            #summary_box.insert(END, "\nNow trimming file %s: %s" % (file_list.index(i) + 1, i))

            if condition_name == "Library":
                # process for library read file
                this_filename = condition_name + global_vars.FILE_FLAGS['Trimmed read file']
                readfile = global_vars.LIBRARY_READ_PATH
            else:
                # process for initial and final reads files
                # get replicate name
                next_char = self.controller.get_rep_char(file_list.index(i))
                #this_filename = condition_name + "-" + next_char + cas9_type + global_vars.FILE_FLAGS['Trimmed read file']
                this_filename = condition_name + cas9_type + next_char + global_vars.FILE_FLAGS['Trimmed read file']
                readfile = i

            # determine line count in file to overcome multicore ban in cutadapt
            line_count = (os.popen('wc -l \'%s\'' % readfile).read()).split(" ", 1)

            # build command pieces
            cmd_prefix = "head \'%s\' -n %s| " % (readfile, line_count[0])
            cmd_main = "cutadapt -g NAAGGACGAAACANNN...GTTTTAGAGCTAGAAAT -> "
            #cmd_f_out = "\'" + current_file + global_vars.FILE_FLAGS['Trimmed read file'] + "\' "  #old filename method
            cmd_f_out = "\'" + this_filename + "\' "
            cmd_args = "-e 0.1 --cores 8 --trimmed-only"

            cmd = cmd_prefix + cmd_main + cmd_f_out + cmd_args

            print("\n\nSaving as filename: %s" % this_filename)

            # # # DEBUGGING # # #
            print(cmd)

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
                    #summary_box.insert(END, output)
                    #summary_box.update_idletasks()
                    if show_in_console:
                        sys.stdout.write(output)
                        sys.stdout.flush()

            # check if we are doing the steps for the Library read file, in which case only one file is required
            # if true, then break the for loop from continuing
            # otherwise for loop continues to iterate over other replicates (for initial and final conditions)
            if condition_name == "Library":
                break

        #summary_box.config(state=DISABLED)

    # function to convert CSV reference library file to fasta format for use with bowtie2
    def convert_to_fasta(self):

        # update status box
        #summary_box = self.txt_summary
        #summary_box.config(state=NORMAL)

        print("\n>Creating directory for reference library...")
        #summary_box.insert(END, "\n>Creating directory for reference library...")

        # change directory
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        # build command pieces
        cmd_prefix = "awk -F \',\' \'{print \">\"$1\"\\n\"$2}\' "
        cmd_f_in = "\'" + global_vars.EXPERIMENT_SETTINGS['Library reference file'] + "\'"
        cmd_f_out = " > \'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + \
                    global_vars.FILE_FLAGS['Fasta library file'] + "\'"
        cmd = cmd_prefix + cmd_f_in + cmd_f_out

        print("\nConverting reference library file from CSV to FASTA...")
        #summary_box.insert(END, "\nConverting reference library file from CSV to FASTA...")

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
                #summary_box.insert(END, output)
                #summary_box.update_idletasks()
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

        print("\nFinished converting reference library file!")
        #summary_box.insert(END, "\nFinished converting reference library file!")
        #summary_box.config(state=DISABLED)

    # function to generate bowtie2 index for reference library
    def build_index(self):

        # update status box
        #summary_box = self.txt_summary
        #summary_box.config(state=NORMAL)

        print("\nBuilding library index...")
        #summary_box.insert(END, "\nBuilding library index...")

        # change directory
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        global_vars.FILE_FLAGS['Bowtie2 index full path'] = \
            os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                         global_vars.FILE_FLAGS['Bowtie2 index name'] + "/")

        if os.path.exists(global_vars.FILE_FLAGS['Bowtie2 index full path']) is False:
            print("\n\'bowtie2-index\' directory not found; creating it now...")
            #summary_box.insert(END, "\n\'bowtie2-index\' directory not found; creating it now...")
            os.mkdir(global_vars.FILE_FLAGS['Bowtie2 index name'], 777)

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
                #summary_box.insert(END, output)
                #summary_box.update_idletasks()
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

        #summary_box.insert(END, ">\nFinished building index!")
        #summary_box.config(state=DISABLED)

    # function to perform alignments for provided sample read files using bowtie2
    def do_bowtie2_alignment(self, file_list, condition_name, cas9_type):

        # update status box
        #summary_box = self.txt_summary
        #summary_box.config(state=NORMAL)

        # change directory
        alignment_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                             global_vars.FILE_FLAGS['Alignments dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(alignment_target_path) is False:
            print("\nCreating directory for aligned files...")
            print("\n%s" % alignment_target_path)
            #summary_box.insert(END, "\nCreating directory for aligned files...")
            os.mkdir(global_vars.FILE_FLAGS['Alignments dir'], 777)

        trimmed_read_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                         global_vars.FILE_FLAGS['Trimmed dir'])

        os.chdir(alignment_target_path)

        # loop through file list for current condition (Library, Initial, Final conditions)
        for i in file_list:
            # update status box
            print("\nNow aligning file %s: %s" % (file_list.index(i) + 1, i + global_vars.FILE_FLAGS['Trimmed read file']))
            #summary_box.insert(END, ">\nNow aligning file %s: %s" % (file_list.index(i) + 1, i + global_vars.FILE_FLAGS['Trimmed read file']))

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
                #output_file = condition_name + "-" + next_char + cas9_type + global_vars.FILE_FLAGS['Aligned bam file']
                output_file = condition_name + cas9_type + next_char + global_vars.FILE_FLAGS['Aligned bam file']

            # build command pieces
            cmd_prefix = "bowtie2 -x \'%s" % global_vars.FILE_FLAGS['Bowtie2 index full path'] + \
                         global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "-" + global_vars.FILE_FLAGS[
                             'Bowtie2 index name'] + "\' -U "
            cmd_f_in = "\'" + input_file + "\'"
            cmd_f_out = "\'" + output_file + "\'"
            cmd_suffix = " --norc -N 1 -p 64 | samtools view -bS - > %s" % cmd_f_out
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
                    #summary_box.insert(END, output)
                    #summary_box.update_idletasks()
                    if show_in_console:
                        sys.stdout.write(output)
                        sys.stdout.flush()

            # check if we are doing the steps for the Library read file,
            # if true, then break the for loop from continuing because it only has one replicate; no looping required
            if condition_name == "Library":
                break

        #summary_box.config(state=DISABLED)

    # function to run mageck's count function to get read counts from alignments for all samples
    def do_mageck_count(self, file_list_1, file_list_2, condition_name_1, condition_name_2):

        # update status box
        #summary_box = self.txt_summary
        #summary_box.config(state=NORMAL)

        # change directory
        counts_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                             global_vars.FILE_FLAGS['Readcounts dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(counts_target_path) is False:
            print("\nCreating directory for read count files...")
            print("\n%s" % counts_target_path)
            #summary_box.insert(END, "\nCreating directory for read count files...")
            os.mkdir(global_vars.FILE_FLAGS['Readcounts dir'], 777)

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
                print("\nNow processing file %s: %s" % (file_list_1.index(i) + 1, i + global_vars.FILE_FLAGS['Aligned bam file']))
                # update status box
                #summary_box.insert(END, "\nNow processing file %s: %s" % (
                #    file_list_1.index(i) + 1, i + global_vars.FILE_FLAGS['Aligned bam file']))

                # get replicate name
                next_char = self.controller.get_rep_char(file_list_1.index(i))

                input_file = os.path.join(aligned_read_path, condition_name_1 + j + next_char + global_vars.FILE_FLAGS['Aligned bam file'])

                cmd_sample_list = cmd_sample_list + condition_name_1 + j + next_char + ","

                cmd_input_file_list = cmd_input_file_list + "\'" + input_file + "\' "

            # final condition labels and files
            for i in file_list_2:
                print("\nNow processing file %s: %s" % (file_list_2.index(i) + 1, i + global_vars.FILE_FLAGS['Aligned bam file']))
                # update status box
                #summary_box.insert(END, "\nNow processing file %s: %s" % (
                #    file_list_2.index(i) + 1, i + global_vars.FILE_FLAGS['Aligned bam file']))

                # get replicate name
                next_char = self.controller.get_rep_char(file_list_2.index(i))

                input_file = os.path.join(aligned_read_path, condition_name_2 + j + next_char + global_vars.FILE_FLAGS['Aligned bam file'])

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
                #summary_box.insert(END, output)
                #summary_box.update_idletasks()
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

        #summary_box.config(state=DISABLED)

    def do_steps_option_2(self):

        # update summary box display
        #self.txt_summary.insert(END, "\nAnalyzing read counts...")
        #self.txt_summary.update_idletasks()

        # run function to analyze reads with TRACS algorithm
        self.process_reads()

    def scroll_down(self):
        self.txt_summary.see("end")
        self.after(1000, self.scroll_down)

    # function to run TRACS algorithm to calculate Enrichment Scores and Enrichment Ratio from read counts table
    def process_reads(self):

        # update status box
        #summary_box = self.txt_summary
        #summary_box.config(state=NORMAL)

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

        print("\nNumber of replicates: %i" % self.num_replicates)  # DEBUG - print total number of replicates

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

        # # # DEBUGGING save to file
        print("\nSaving to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.raw_counts.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "rawcounts.csv",
                               index=True, sep='\t', float_format='%.16f')

        print("\nDetermining average...")

        # call function to calculate average of log2FC for every gene for each log2FC sample column
        # calculates a single log2FC from all the guides for each gene
        self.avg = pd.DataFrame()  # new dataframe to hold averages
        self.avg_calculator(name = global_vars.COLUMN_NAMES['Library'])
        self.avg_calculator(name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'])
        self.avg_calculator(name = global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        # # # DEBUGGING save to file
        print("\nSaving to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.avg.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".average.csv", index=True,
                        sep='\t',
                        float_format='%.16f')

        # # # ENRICHMENT SCORE # # #
        print("\nCalculating sgw...")

        # determine sgw for each gene
        self.sgw = pd.DataFrame()  # new dataframe to hold sgw
        self.get_sgw()

        print("\nCalculating rank...")
        # determine rank for each sample
        self.get_rank()

        # # # DEBUGGING save to file
        print("\nSaving to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.sgw.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".sgw.csv", index=True,
                        sep='\t',
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

        # # # DEBUGGING save to file
        print("\nSaving to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.ES.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".ES.csv", index=True,
                               sep='\t',
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

        # save to file
        print("\nSaving to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.ES.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".csv",
                       index=True, sep='\t', float_format='%.6f')

        print("\nAll done!")

    # function to calculate log2FC for Library. Log2FC = log2[Initial Cas9- / Library]
    def lib_log2FC_calculator(self):

        # loop through all initial Cas9- replicates
        for col in range(0, int(self.num_replicates)):
            # get name of library column
            library_norm_name = global_vars.COLUMN_NAMES['Library'] + global_vars.COLUMN_NAMES['Normalized']

            # get name of current initial cas9- replicate
            current_col_name = global_vars.EXPERIMENT_SETTINGS['Initial condition name'] + global_vars.COLUMN_NAMES['Cas9 neg'] + self.controller.get_rep_char(col) + global_vars.COLUMN_NAMES['Normalized']
                                #Adh_Cas9neg-A.norm
                                # Adh-A-Cas9-neg
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
    # ie. 6 sgRNA ID log2FC --> reduced to single log2FC per gene
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
    app.resizable(False, False)
    app.mainloop()



