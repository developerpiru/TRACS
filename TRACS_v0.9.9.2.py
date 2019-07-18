import tkinter as tk
from tkinter import font as tkfont
from tkinter import *
from tkinter import filedialog
from tkinter import ttk
from tkinter.ttk import Frame, Label, Button, Style
import os
import subprocess, sys
import pandas as pd
import numpy as np
import global_vars, plugins
import math

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
# - DONE - get log2FC per sample
# - DONE - get sgw per sample
# - DONE - get rank per sample
# - DONE - get Library, Initial, and Final ES per gene
# - DONE - get ES ratio (ER) (log2 fold change of Initial and Final ES) per gene
# TRACS algorithm functions work for n replicates

# # # TO DO # # #
# output R file for web-based viewing tool
# implement multithreading
# move summary_box outputs to separate function call ???? --- stuck in infinite loop with bowtie2 function
# implement error catching

# # # obsolete/deprecated method # # #
# split filename function (plugins.split_filename) working but cutadapt fails to get fullpath


# # # OUTPUTS # # #
#
# output of raw read counts from mageck:
#
# | sgRNA ID | Gene | Cas9+ Initial A | Cas9+ Initial B | Cas9+ Initial C | Cas9+ Final A | Cas9+ Final B | Cas9+ Final C | Cas9- Initial A | Cas9- Initial B | Cas9- Initial C | Cas9- Final A | Cas9- Final B | Cas9- Final C |
#
# sgRNA ID
# Gene name
# n replicates of Cas9+ Initial
# n replicates of Cas9+ Final
# n replicates of Cas9- Initial
# n replicates of Cas9- Final


# Currently disabled file loading in process_inputs() #


# version #
app_version = "0.9.9.2"


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
        self.winfo_toplevel().title("TRACS")

        # set fonts
        self.title_font = tkfont.Font(family='Helvetica', size=18, weight="bold")
        self.label_font = tkfont.Font(family='Helvetica', size=10, weight="bold")
        self.label_font_small = tkfont.Font(family='Helvetica', size=8, weight="bold")
        self.body_font = tkfont.Font(family='Helvetica', size=12, weight="bold")

        # main container
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # define page frames
        self.frames = {}
        # create frames for each page
        for F in (MainPage, NewExpSettings, LoadFiles, LoadCas9NegFiles_v2, SummaryPage, TrimmingReads, FileExtractor, AnalyzeCounts):
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

    def browse_folder(self, frame_name, tar_name):
        target_obj = getattr(self.frames[frame_name], tar_name)
        folder = filedialog.askdirectory(initialdir=os.getcwd(), title='Select folder...')
        target_obj.delete(0, tk.END)
        target_obj.insert(tk.END, folder)

    def browse_file(self, multi_files, frame_name, tar_name):
        target_obj = getattr(self.frames[frame_name], tar_name)
        if multi_files == "TRUE":
            files = filedialog.askopenfilenames(title='Select files...',
                                                filetypes=(("FASTQ files", "*.fastq"),
                                                           ("GZipped Fastq", "*fastq.gz"), ("All files", "*")))
            for item in files:
                #target_obj.insert(END, plugins.split_filename(item))      #split filename
                target_obj.insert(END, item)
        else:
            file = filedialog.askopenfilename(title='Select file...',
                                              filetypes=(("FASTQ files", "*.fastq"),
                                                         ("GZipped Fastq", "*fastq.gz"), ("CSV", "*.csv"),
                                                         ("All files", "*")))
            target_obj.delete(0, END)
            #target_obj.insert(END, plugins.split_filename(file))      #split filename
            target_obj.insert(END, file)

    def delete_selected(self, frame_name, tar_name):
        target_obj = getattr(self.frames[frame_name], tar_name)
        selected = target_obj.curselection()
        for i in selected[::-1]:
            target_obj.delete(i)

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

    def get_rep_names(self, x):
        sep = global_vars.REPLICATE_SELECTED[0]
        next_char = chr(ord(global_vars.REPLICATE_SELECTED[1]) + int(x)).upper()
        return sep, next_char

    def extract_gzip(self):
        FILE_PATH_LIST = self.frames['FileExtractor'].file_list.get(0, END)
        summary_box = self.frames['FileExtractor'].txt_summary
        cmd = "gunzip -k"
        for i in FILE_PATH_LIST:
            cmd = cmd + " %s" % i

        summary_box.insert(END, "\n>" + cmd)
        cmd_output = os.popen(cmd).read()
        summary_box.insert(END, cmd_output)
        #os.system(cmd)

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
        lbl_info = Label(fr_1, text="TRACS",
                         wraplength="800", justify="left", font=controller.title_font)
        lbl_info.pack(padx=5, pady=5)

        fr_2 = Frame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=X, side=TOP, padx=5, pady=5, anchor=N, expand=True)
        lbl_info = Label(fr_2, text="Toolset for Ranked Analysis of CRISPR Screens",
                         wraplength="800", justify="left", font=controller.body_font)
        lbl_info.pack(side=TOP, padx=5, pady=5)
        lbl_version = Label(fr_2, text="Version " + app_version,
                            wraplength="800", justify="left", font=controller.body_font)
        lbl_version.pack(side=TOP, padx=5, pady=5)

        btnNewExp = Button(frame, text="New Experiment", command=lambda: controller.show_frame("NewExpSettings"))
        btnNewExp.pack(padx=5, pady=5)

        btnFileExtractor = Button(frame, text="Extract Files", command=lambda: controller.show_frame("FileExtractor"))
        btnFileExtractor.pack(padx=5, pady=5)

        self.pack(fill=BOTH, expand=True)

        btnBack = Button(self, text="Exit", command=lambda: app.destroy())
        btnBack.pack(side=RIGHT, padx=5, pady=5)


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
        lbl_heading = Label(fr_header_1, text="Step 1: New Experiment Settings", font=controller.title_font)
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
        self.txt_name_initial.bind("<Key>", self.update_labels)

        # name of final group
        fr_5 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_5.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_name_final = Label(fr_5, text="Name of final condition (Tf)", width="40", font=controller.label_font)
        lbl_name_final.pack(side=LEFT, padx=5, pady=5)
        self.txt_name_final = Entry(fr_5)
        self.txt_name_final.pack(fill=X, padx=5, pady=5, expand=True)
        self.txt_name_final.bind("<Key>", self.update_labels)

        # replicate naming scheme
        fr_6 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_6.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)
        lbl_replicates = Label(fr_6, text="Replicate labeling scheme", width="40", font=controller.label_font)
        lbl_replicates.pack(side=LEFT, padx=5, pady=5)

        # dropdown box
        self.opt_repnaming = ttk.Combobox(fr_6, state="readonly",  width="8",
                                          values=["-1, -2, -3", "_1, _2, _3", "-A, -B, -C", "_A, _B, _C"])
        self.opt_repnaming.pack(side=LEFT, padx=5, pady=5)
        self.opt_repnaming.current(0)
        self.lbl_rep_sample = Label(fr_6, text="Initial-1 and Final-1", width="40", font=controller.label_font_small)
        self.lbl_rep_sample.pack(fill=X, padx=5, pady=5, anchor=N, expand=True)

        self.opt_repnaming.bind("<<ComboboxSelected>>", self.update_labels)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: self.controller.show_frame("LoadFiles"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("MainPage"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    def update_labels(self, event):
        suffix = self.opt_repnaming.get()

        if self.txt_name_initial.get() == "":
            initial_name = "Initial"
        else:
            initial_name = self.txt_name_initial.get()

        if self.txt_name_final.get() == "":
            final_name = "Final"
        else:
            final_name = self.txt_name_final.get()

        if suffix == "-1, -2, -3":
            text = "%s-1, %s-2" % (initial_name, initial_name) + " and %s-1, %s-2" % (final_name, final_name)
        elif suffix == "_1, _2, _3":
            text = "%s_1, %s_2" % (initial_name, initial_name) + " and %s_1, %s_2" % (final_name, final_name)
        elif suffix == "-A, -B, -C":
            text = "%s-A, %s-B" % (initial_name, initial_name) + " and %s-A, %s-B" % (final_name, final_name)
        elif suffix == "_A, _B, _C":
            text = "%s_A, %s_B" % (initial_name, initial_name) + " and %s_A, %s_B" % (final_name, final_name)

        self.lbl_rep_sample.config(text=text)
        self.lbl_rep_sample.update_idletasks()


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
        lbl_heading = Label(fr_header_1, text="Step 2: Load Cas9+ Files", font=controller.title_font)
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
        lbl_reflibrary = Label(fr_2, text="Library reference file", width="40", font=controller.label_font)
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
        lbl_libreads = Label(fr_6, text="Library reads file", width="40", font=controller.label_font)
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
        lbl_i_info = Label(fr_4, text="Initial condition reads files", font=controller.label_font)
        lbl_i_info.pack(padx=5, pady=5, side=LEFT)

        # frame for final label
        fr_5 = Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_5.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_f_info = Label(fr_5, text="Final condition reads files", font=controller.label_font)
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
        # move up
        btn_i_up = Button(initial_btn_container, text="Up",
                            command=lambda: self.controller.move_up(frame_name="LoadFiles",
                                                                    tar_name="lst_initial_files"))
        btn_i_up.pack(side=TOP, padx=0, pady=0)
        # move down
        btn_i_down = Button(initial_btn_container, text="Down",
                              command=lambda: self.controller.move_down(frame_name="LoadFiles",
                                                                        tar_name="lst_initial_files"))
        btn_i_down.pack(side=TOP, padx=0, pady=0)

        # frame for initial listbox
        final_list_container = Frame(list_container, relief=controller.frame_relief,
                                       borderwidth=controller.frame_borderwidth)
        final_list_container.pack(fill=BOTH, padx=5, pady=5, side=LEFT, expand=True)
        # frame for initial listbox buttons
        final_btn_container = Frame(final_list_container)
        final_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # initial listbox
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
        # move up
        btn_f_up = Button(final_btn_container, text="Up",
                          command=lambda: self.controller.move_up(frame_name="LoadFiles",
                                                                  tar_name="lst_final_files"))
        btn_f_up.pack(side=TOP, padx=0, pady=0)
        # move down
        btn_f_down = Button(final_btn_container, text="Down",
                            command=lambda: self.controller.move_down(frame_name="LoadFiles",
                                                                      tar_name="lst_final_files"))
        btn_f_down.pack(side=TOP, padx=0, pady=0)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("LoadCas9NegFiles_v2"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("NewExpSettings"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)


# new function to load Cas9- files
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
        lbl_heading = Label(fr_header_1, text="Step 3: Load Cas9- Files", font=controller.title_font)
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
        lbl_info = Label(fr_1, text="Select Cas9- negative control read files",
                         wraplength="800", justify="left", font=controller.body_font)
        lbl_info.pack(side=LEFT, padx=5, pady=5)

        # container frame for listbox labels
        fr_3 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_3.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=True)

        # frame for initial label
        fr_4 = Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_4.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_i_info = Label(fr_4, text="Initial condition read files", font=controller.label_font)
        lbl_i_info.pack(padx=5, pady=5, side=LEFT)

        # frame for final label
        fr_5 = Frame(fr_3, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_5.pack(fill=X, padx=5, pady=5, side=LEFT, expand=True)
        lbl_f_info = Label(fr_5, text="Final condition read files", font=controller.label_font)
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
        # move up
        btn_i_up = Button(initial_btn_container, text="Up",
                            command=lambda: self.controller.move_up(frame_name="LoadCas9NegFiles_v2",
                                                                    tar_name="lst_initial_files"))
        btn_i_up.pack(side=TOP, padx=0, pady=0)
        # move down
        btn_i_down = Button(initial_btn_container, text="Down",
                              command=lambda: self.controller.move_down(frame_name="LoadCas9NegFiles_v2",
                                                                        tar_name="lst_initial_files"))
        btn_i_down.pack(side=TOP, padx=0, pady=0)

        # frame for initial listbox
        final_list_container = Frame(list_container, relief=controller.frame_relief,
                                       borderwidth=controller.frame_borderwidth)
        final_list_container.pack(fill=BOTH, padx=5, pady=5, side=LEFT, expand=True)
        # frame for initial listbox buttons
        final_btn_container = Frame(final_list_container)
        final_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        # initial listbox
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
        # move up
        btn_f_up = Button(final_btn_container, text="Up",
                          command=lambda: self.controller.move_up(frame_name="LoadCas9NegFiles_v2",
                                                                  tar_name="lst_final_files"))
        btn_f_up.pack(side=TOP, padx=0, pady=0)
        # move down
        btn_f_down = Button(final_btn_container, text="Down",
                            command=lambda: self.controller.move_down(frame_name="LoadCas9NegFiles_v2",
                                                                      tar_name="lst_final_files"))
        btn_f_down.pack(side=TOP, padx=0, pady=0)

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
        global_vars.EXPERIMENT_SETTINGS['Library reference file'] = \
            self.controller.frames['LoadFiles'].txt_reflibrary.get()  # number of replicates

        # get selected replicate naming scheme
        global_vars.REPLICATE_SELECTED = self.controller.frames['NewExpSettings'].opt_repnaming.get()

        # get library read fastq file
        global_vars.LIBRARY_READ_PATH = self.controller.frames['LoadFiles'].txt_libreads.get()
        global_vars.EXPERIMENT_SETTINGS['Library read file'] = \
            self.controller.frames['LoadFiles'].txt_libreads.get()

        print("\nLIBRARY_READ_PATH = %s" % global_vars.LIBRARY_READ_PATH)
        print("\nLibrary read file = %s" % global_vars.EXPERIMENT_SETTINGS['Library read file'])

        # get Initial condition Cas9+ fastq files
        global_vars.LIST_INITIAL_PATHS = self.controller.frames['LoadFiles'].lst_initial_files.get(0, END)
        # get Final condition Cas9+ fastq files
        global_vars.LIST_FINAL_PATHS = self.controller.frames['LoadFiles'].lst_final_files.get(0, END)
        # get Initial condition Cas9- fastq files
        global_vars.LIST_CAS9_NEG_I_PATHS = self.controller.frames['LoadCas9NegFiles_v2'].lst_initial_files.get(0, END)
        # get Final condition Cas9- fastq files
        global_vars.LIST_CAS9_NEG_F_PATHS = self.controller.frames['LoadCas9NegFiles_v2'].lst_final_files.get(0, END)

        summary_box = self.controller.frames['SummaryPage'].txt_summary
        summary_box.config(state=NORMAL)
        summary_box.delete(1.0, END)

        # output inputs to summary text box
        for keys, values in global_vars.EXPERIMENT_SETTINGS.items():
            tmp_text = "%s: %s\n" % (keys, values)
            summary_box.insert(END, tmp_text)

        self.print_files(tar_name=global_vars.LIBRARY_READ_PATH,
                         condition="Library", cas9_type="+")
        self.print_files(tar_name=global_vars.LIST_INITIAL_PATHS,
                         condition=global_vars.EXPERIMENT_SETTINGS['Initial condition name'], cas9_type="+")
        self.print_files(tar_name=global_vars.LIST_FINAL_PATHS,
                         condition=global_vars.EXPERIMENT_SETTINGS['Final condition name'], cas9_type="+")
        self.print_files(tar_name=global_vars.LIST_CAS9_NEG_I_PATHS,
                         condition=global_vars.EXPERIMENT_SETTINGS['Initial condition name'], cas9_type="-")
        self.print_files(tar_name=global_vars.LIST_CAS9_NEG_F_PATHS,
                         condition=global_vars.EXPERIMENT_SETTINGS['Final condition name'], cas9_type="-")

        summary_box.config(state=DISABLED)

        # show summary page after processing all inputs
        self.controller.show_frame("SummaryPage")

    def print_files(self, tar_name, condition, cas9_type):
        summary_box = self.controller.frames['SummaryPage'].txt_summary

        if condition == "Library":
            tmp_text = "\n\nInitial Library read file"
        else:
            tmp_text = "\n\nSample files for Cas9%s '%s' condition\n" % (cas9_type, condition)

        summary_box.insert(END, tmp_text)
        for i in tar_name:
            # sep, next_char = self.controller.get_rep_names(file_list.index(i))
            # this_filename = condition + sep + str(next_char) + global_vars.FILE_FLAGS['Trimmed read file']

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
        lbl_heading = Label(fr_header_1, text="Step 4: Summary", font=controller.title_font)
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
                         wraplength="800", justify="left", font=controller.body_font)
        lbl_info.pack(fill=X, side=LEFT, padx=5, pady=5)

        # summary text box
        fr_2 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        self.txt_summary = Text(fr_2, height="10")
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
        lbl_heading = Label(fr_header_1, text="Step 5: Running analyses", font=controller.title_font)
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
        lbl_info = Label(fr_1, text="Please wait as reads are being trimmed...",
                         wraplength="800", justify="left", font=controller.body_font)
        lbl_info.pack(fill=X, side=LEFT, padx=5, pady=5)

        # summary text box
        fr_2 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        self.txt_summary = Text(fr_2, height="10")
        self.txt_summary.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        scrl_txt_summary = Scrollbar(fr_2)
        scrl_txt_summary.pack(side=RIGHT, fill=Y)
        scrl_txt_summary.config(command=self.txt_summary.yview)
        self.txt_summary.config(yscrollcommand=scrl_txt_summary.set)

        # btnTrim = Button(fr_2, text="Trim Files", command=lambda: self.do_trimming_v3(file_list=global_vars.LIST_INITIAL_PATHS))
        btnTrim = Button(fr_2, text="Start",
                         command=lambda: self.do_steps())
        btnTrim.pack(side=RIGHT)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("AnalyzeCounts"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("SummaryPage"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    def do_steps(self):

        # convert to fasta library
        self.convert_to_fasta()

        # trim reads
        print("\nTrimming reads...")
        self.txt_summary.insert(END, "\nTrimming reads...")
        self.txt_summary.update_idletasks()

        # Library reads
        self.do_trimming_v2(file_list=global_vars.LIBRARY_READ_PATH,
                            condition_name="Library",
                            cas9_type="")

        # Cas9+
        self.do_trimming_v2(file_list=global_vars.LIST_INITIAL_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[0])
        self.do_trimming_v2(file_list=global_vars.LIST_FINAL_PATHS,
                            condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                            cas9_type=global_vars.CAS9_TYPE_NAMING[0])

        # Cas9-
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
        self.txt_summary.insert(END, "\nAligning reads with library index...")
        self.txt_summary.update_idletasks()

        # Library reads
        self.do_bowtie2_alignment(file_list=global_vars.LIBRARY_READ_PATH,
                                  condition_name="Library",
                                  cas9_type="")

        # Cas9+
        self.do_bowtie2_alignment(file_list=global_vars.LIST_INITIAL_PATHS,
                                  condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                                  cas9_type=global_vars.CAS9_TYPE_NAMING[0])
        self.do_bowtie2_alignment(file_list=global_vars.LIST_FINAL_PATHS,
                                  condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                                  cas9_type=global_vars.CAS9_TYPE_NAMING[0])

        # Cas9-
        self.do_bowtie2_alignment(file_list=global_vars.LIST_CAS9_NEG_I_PATHS,
                                  condition_name=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                                  cas9_type=global_vars.CAS9_TYPE_NAMING[1])
        self.do_bowtie2_alignment(file_list=global_vars.LIST_CAS9_NEG_F_PATHS,
                                  condition_name=global_vars.EXPERIMENT_SETTINGS['Final condition name'],
                                  cas9_type=global_vars.CAS9_TYPE_NAMING[1])

        # run mageck count to get raw read counts
        self.do_mageck_count(file_list_1=global_vars.LIST_INITIAL_PATHS,
                             file_list_2=global_vars.LIST_FINAL_PATHS,
                             condition_name_1=global_vars.EXPERIMENT_SETTINGS['Initial condition name'],
                             condition_name_2=global_vars.EXPERIMENT_SETTINGS['Final condition name'])

        print("\n\nDone alignments and read count generation")

    def print_to_summary_box(self, output):
        summary_box = self.txt_summary
        summary_box.config(state=NORMAL)
        summary_box.insert(END, output)
        summary_box.config(state=DISABLED)
        summary_box.update_idletasks()
        summary_box.see("end")
        #self.scroll_down()

    def scroll_down(self):
        self.txt_summary.see("end")
        self.after(1000, self.scroll_down)

    def do_trimming_v2(self, file_list, condition_name, cas9_type):
        summary_box = self.txt_summary
        summary_box.config(state=NORMAL)

        # change directory
        trimmed_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                           global_vars.FILE_FLAGS['Trimmed dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(trimmed_target_path) is False:
            print("\nCreating directory for aligned files...")
            print("\n%s" % trimmed_target_path)
            summary_box.insert(END, "\nCreating directory for aligned files...")
            os.mkdir(global_vars.FILE_FLAGS['Trimmed dir'], 777)

        os.chdir(trimmed_target_path)

        for i in file_list:
            # update status box
            summary_box.insert(END, "\nNow trimming file %s: %s" % (file_list.index(i) + 1, i))
            # command pieces

            current_file = plugins.split_filename(i)

            print("\nPRINT condition_name = %s" % condition_name)

            if condition_name == "Library":
                # process for library read file
                this_filename = condition_name + global_vars.FILE_FLAGS['Trimmed read file']
                readfile = global_vars.EXPERIMENT_SETTINGS['Library read file']
            else:
                # process for initial and final reads files
                # get replicate name
                sep, next_char = self.controller.get_rep_names(file_list.index(i))
                this_filename = condition_name + sep + str(next_char) + cas9_type + global_vars.FILE_FLAGS['Trimmed read file']
                readfile = i

            # determine line count in file to overcome multicore ban in cutadapt
            line_count = (os.popen('wc -l \'%s\'' % readfile).read()).split(" ", 1)

            cmd_prefix = "head \'%s\' -n %s| " % (readfile, line_count[0])
            cmd_main = "cutadapt -g NAAGGACGAAACANNN...GTTTTAGAGCTAGAAAT -> "
            #cmd_f_out = "\'" + current_file + global_vars.FILE_FLAGS['Trimmed read file'] + "\' "  #old filename method
            cmd_f_out = "\'" + this_filename + "\' "
            cmd_args = "-e 0.1 --cores 8 --trimmed-only"

            cmd = cmd_prefix + cmd_main + cmd_f_out + cmd_args

            print("\n\nSaving as filename: %s" % this_filename)

            print(cmd)

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
                    summary_box.insert(END, output)
                    summary_box.update_idletasks()
                    if show_in_console:
                        sys.stdout.write(output)
                        sys.stdout.flush()

            # check if we are doing the steps for the Library read file,
            # in which case only one file is required
            # if true, then break the for loop from continuing
            # otherwise for loop iterates over each char of global_vars.LIBRARY_READ_PATH
            if condition_name == "Library":
                break

        summary_box.config(state=DISABLED)

    def convert_to_fasta(self):
        # awk -F ',' '{print ">"$1"\n"$2}' LibraryAB.csv > LibraryAB_fasta.fa
        summary_box = self.txt_summary
        summary_box.config(state=NORMAL)

        print("\n>Creating directory for reference library...")
        summary_box.insert(END, "\n>Creating directory for reference library...")

        # change directory
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        # command pieces
        cmd_prefix = "awk -F \',\' \'{print \">\"$1\"\\n\"$2}\' "
        cmd_f_in = "\'" + global_vars.EXPERIMENT_SETTINGS['Library reference file'] + "\'"
        cmd_f_out = " > \'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + \
                    global_vars.FILE_FLAGS['Fasta library file'] + "\'"
        cmd = cmd_prefix + cmd_f_in + cmd_f_out

        print("\nConverting reference library file from CSV to FASTA...")
        summary_box.insert(END, "\nConverting reference library file from CSV to FASTA...")

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
                summary_box.insert(END, output)
                summary_box.update_idletasks()
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

        print("\nFinished converting reference library file!")
        summary_box.insert(END, "\nFinished converting reference library file!")
        summary_box.config(state=DISABLED)

    def build_index(self):
        summary_box = self.txt_summary
        summary_box.config(state=NORMAL)

        print("\nBuilding library index...")
        summary_box.insert(END, "\nBuilding library index...")

        # change directory
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        global_vars.FILE_FLAGS['Bowtie2 index full path'] = \
            os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                         global_vars.FILE_FLAGS['Bowtie2 index name'] + "/")

        if os.path.exists(global_vars.FILE_FLAGS['Bowtie2 index full path']) is False:
            print("\n\'bowtie2-index\' directory not found; creating it now...")
            summary_box.insert(END, "\n\'bowtie2-index\' directory not found; creating it now...")
            os.mkdir(global_vars.FILE_FLAGS['Bowtie2 index name'], 777)

        os.chdir(global_vars.FILE_FLAGS['Bowtie2 index name'])

        # command pieces
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
                summary_box.insert(END, output)
                summary_box.update_idletasks()
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

        summary_box.insert(END, ">\nFinished building index!")
        summary_box.config(state=DISABLED)

    def do_bowtie2_alignment(self, file_list, condition_name, cas9_type):
        summary_box = self.txt_summary
        summary_box.config(state=NORMAL)

        # change directory
        alignment_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                             global_vars.FILE_FLAGS['Alignments dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(alignment_target_path) is False:
            print("\nCreating directory for aligned files...")
            print("\n%s" % alignment_target_path)
            summary_box.insert(END, "\nCreating directory for aligned files...")
            os.mkdir(global_vars.FILE_FLAGS['Alignments dir'], 777)

        trimmed_read_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                         global_vars.FILE_FLAGS['Trimmed dir']) #  + "/"

        os.chdir(alignment_target_path)

        for i in file_list:
            # update status box
            summary_box.insert(END, ">\nNow aligning file %s: %s" % (file_list.index(i) + 1, i + global_vars.FILE_FLAGS['Trimmed read file']))

            if condition_name == "Library":
                # process for library read file
                input_file = os.path.join(trimmed_read_path, condition_name + global_vars.FILE_FLAGS['Trimmed read file'])
                output_file = condition_name + global_vars.FILE_FLAGS['Aligned bam file']
            else:
                # process for initial and final read files
                # get replicate name
                sep, next_char = self.controller.get_rep_names(file_list.index(i))
                input_file = os.path.join(trimmed_read_path, condition_name + sep + str(next_char) + cas9_type + global_vars.FILE_FLAGS['Trimmed read file'])
                output_file = condition_name + sep + str(next_char) + cas9_type + global_vars.FILE_FLAGS['Aligned bam file']

            # command pieces
            cmd_prefix = "bowtie2 -x \'%s" % global_vars.FILE_FLAGS['Bowtie2 index full path'] + \
                         global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "-" + global_vars.FILE_FLAGS[
                             'Bowtie2 index name'] + "\' -U "
            # old method
            # cmd_f_in = "\'../" + input_file + global_vars.FILE_FLAGS['Trimmed read file'] + "\'"
            # cmd_f_out = "\'" + current_file + global_vars.FILE_FLAGS['Aligned bam file'] + "\'"

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
                    summary_box.insert(END, output)
                    summary_box.update_idletasks()
                    if show_in_console:
                        sys.stdout.write(output)
                        sys.stdout.flush()

            # check if we are doing the steps for the Library read file,
            # if true, then break the for loop from continuing
            if condition_name == "Library":
                break

        summary_box.config(state=DISABLED)

    def do_mageck_count(self, file_list_1, file_list_2, condition_name_1, condition_name_2):
        summary_box = self.txt_summary
        summary_box.config(state=NORMAL)

        # change directory
        counts_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                             global_vars.FILE_FLAGS['Readcounts dir'] + "/")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])

        if os.path.exists(counts_target_path) is False:
            print("\nCreating directory for read count files...")
            print("\n%s" % counts_target_path)
            summary_box.insert(END, "\nCreating directory for read count files...")
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

        # loop for Cas9+ and Cas-
        for j in global_vars.CAS9_TYPE_NAMING:
            # initial condition labels and files

            print("\n *** j = %s" % j)
            print("\n *** global_vars.CAS9_TYPE_NAMING.index(j) = %s" % global_vars.CAS9_TYPE_NAMING.index(j))
            print("\n *** length of variable = %s" % len(global_vars.CAS9_TYPE_NAMING))

            for i in file_list_1:
                # update status box
                summary_box.insert(END, "\nNow processing file %s: %s" % (
                    file_list_1.index(i) + 1, i + global_vars.FILE_FLAGS['Aligned bam file']))

                # get replicate name
                sep, next_char = self.controller.get_rep_names(file_list_1.index(i))

                input_file = os.path.join(aligned_read_path, condition_name_1 + sep + str(next_char) + j + global_vars.FILE_FLAGS['Aligned bam file'])

                cmd_sample_list = cmd_sample_list + condition_name_1 + sep + str(next_char) + j + ","

                cmd_input_file_list = cmd_input_file_list + "\'" + input_file + "\' "

            # final condition labels and files
            for i in file_list_2:
                # update status box
                summary_box.insert(END, "\nNow processing file %s: %s" % (
                    file_list_2.index(i) + 1, i + global_vars.FILE_FLAGS['Aligned bam file']))

                # get replicate name
                sep, next_char = self.controller.get_rep_names(file_list_2.index(i))

                input_file = os.path.join(aligned_read_path, condition_name_2 + sep + str(next_char) + j + global_vars.FILE_FLAGS['Aligned bam file'])

                cmd_sample_list = cmd_sample_list + condition_name_2 + sep + str(next_char) + j + ","

                cmd_input_file_list = cmd_input_file_list + "\'" + input_file + "\' "

        # command pieces
        cmd_prefix = "mageck count -l \'" + global_vars.EXPERIMENT_SETTINGS['Library reference file'] + "\'"
        cmd_namtag = " -n \'" + global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "\'"
        cmd_sampletags = " --sample-label \'" + cmd_sample_list[:-1] + "\'"
        cmd_readfiles = " --fastq " + cmd_input_file_list
        cmd = cmd_prefix + cmd_namtag + cmd_sampletags + cmd_readfiles

        print("\n\n**Mageck command output:\n %s" % cmd)

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
                summary_box.insert(END, output)
                summary_box.update_idletasks()
                if show_in_console:
                    sys.stdout.write(output)
                    sys.stdout.flush()

        summary_box.config(state=DISABLED)


class AnalyzeCounts(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        # main frame for page body content
        fr_main_1 = Frame(self, relief=GROOVE, borderwidth=0)
        fr_main_1.pack(fill=BOTH, padx=5, pady=5, expand=True)

        # heading
        fr_header_1 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_header_1.pack(fill=BOTH, padx=5, pady=5, anchor=N)
        lbl_heading = Label(fr_header_1, text="Step 6: Analyze Read counts", font=controller.title_font)
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
        lbl_info = Label(fr_1, text="Please wait while read counts are analyzed...",
                         wraplength="800", justify="left", font=controller.body_font)
        lbl_info.pack(fill=X, side=LEFT, padx=5, pady=5)

        # summary text box
        fr_2 = Frame(fr_main_1, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth)
        fr_2.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        self.txt_summary = Text(fr_2, height="10")
        self.txt_summary.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        scrl_txt_summary = Scrollbar(fr_2)
        scrl_txt_summary.pack(side=RIGHT, fill=Y)
        scrl_txt_summary.config(command=self.txt_summary.yview)
        self.txt_summary.config(yscrollcommand=scrl_txt_summary.set)

        # btnTrim = Button(fr_2, text="Trim Files", command=lambda: self.do_trimming_v3(file_list=global_vars.LIST_INITIAL_PATHS))
        btnTest = Button(fr_2, text="Read counts",
                         command=lambda: self.do_steps())
        btnTest.pack(side=RIGHT)

        # navigation buttons
        btn_Next = Button(self, text="Next", command=lambda: controller.show_frame("AnalyzeCounts"))
        btn_Next.pack(side=RIGHT)
        btn_Back = Button(self, text="Back", command=lambda: controller.show_frame("TrimmingReads"))
        btn_Back.pack(side=RIGHT, padx=5, pady=5)
        btn_Back.pack(side=RIGHT, padx=5, pady=5)

    def do_steps(self):

        # update summary box display
        self.txt_summary.insert(END, "\nAnalyzing read counts...")
        self.txt_summary.update_idletasks()

        # run function to analyze reads with TRACS algorithm
        self.process_reads()

    def print_to_summary_box(self, output):
        summary_box = self.txt_summary
        summary_box.config(state=NORMAL)
        summary_box.insert(END, output)
        summary_box.config(state=DISABLED)
        summary_box.update_idletasks()
        summary_box.see("end")
        # self.scroll_down()

    def scroll_down(self):
        self.txt_summary.see("end")
        self.after(1000, self.scroll_down)

    def process_reads(self):
        summary_box = self.txt_summary
        summary_box.config(state=NORMAL)

        # change directory
        reads_target_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'],
                                           global_vars.FILE_FLAGS['Readcounts dir'] + "/")

        os.chdir(reads_target_path)

        readsfile = global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".count.txt"
        #readsfile = "DemoJuly15.count.txt"

        # Table format:
        # sgRNA ID
        # Gene name
        # Library
        # n replicates of Cas9+ Initial
        # n replicates of Cas9+ Final
        # n replicates of Cas9- Initial
        # n replicates of Cas9- Final

        # load counts into pandas dataframe
        self.raw_counts = pd.read_csv(readsfile, delimiter='\t')
        # set gene names as index
        self.raw_counts.index = self.raw_counts['Gene']
        # remove gene and sgRNA columns
        del self.raw_counts['Gene'], self.raw_counts['sgRNA']
        self.raw_counts = self.raw_counts.sort_index()
        print(self.raw_counts)

        # # # NORMALIZING READS # # #

        # determine number of columns (which is equal to number of samples, including library)
        num_samples = len(self.raw_counts.columns)
        print("\nNumber of samples: %i" % num_samples)  # variable not used anymore

        # add 1 to each read count to prevent divide by 0 errors when calculating log2FC
        # self.raw_counts = self.raw_counts + 1
        self.raw_counts += 1

        print("\nNormalizing read counts...")

        # get sum of each column (total reads)
        total_reads = self.raw_counts.sum(axis=0)
        print(total_reads)

        # normalize by total reads for each column
        i = 0
        for col in self.raw_counts.columns:
            name = col + ".norm"
            print("\nCurrent column: %s and sum of column is %i" % (name, total_reads[i]))
            self.raw_counts[name] = (self.raw_counts[col] * 1.0) / total_reads[i]
            # format for correct digits
            self.raw_counts[name] = self.raw_counts[name].map(lambda x: '%.16f' % x)
            self.raw_counts[name] = self.raw_counts[name].astype(float)
            i += 1

        print(self.raw_counts)
        print(self.raw_counts.dtypes)

        # # # CALCULATE LOG2 FOLD CHANGE # # #
        # # the log2FC is Cas9+ / Cas9- = Enrichment Score

        # Table format:
        # sgRNA ID
        # Gene name
        # Library
        # n replicates of Cas9+ Initial
        # n replicates of Cas9+ Final
        # n replicates of Cas9- Initial
        # n replicates of Cas9- Final

        # determine number of replicates for Cas9+ Initial condition
        cas9p_intial_reps = len(global_vars.LIST_INITIAL_PATHS)
        # determine number of replicates for Cas9+ Final condition
        cas9p_final_reps = len(global_vars.LIST_FINAL_PATHS)
        # determine number of replicates for Cas9- Initial condition
        cas9n_intial_reps = len(global_vars.LIST_CAS9_NEG_I_PATHS)
        # determine number of replicates for Cas9- Final condition
        cas9n_final_reps = len(global_vars.LIST_CAS9_NEG_F_PATHS)

        # determine number of replicates
        self.num_replicates = len(global_vars.LIST_INITIAL_PATHS)

        # for debugging
        # print("\n# of initial cas9+: %i, should be 3" % cas9p_intial_reps)
        # print("\n# of final cas9+: %i, should be 3" % cas9p_final_reps)
        # print("\n# of initial cas9-: %i, should be 3" % cas9n_intial_reps)
        # print("\n# of final cas9-: %i, should be 3" % cas9n_final_reps)

        # determine total number of samples (all Cas9+ + all Cas9- columns)
        total_sample_columns = cas9p_intial_reps + cas9p_final_reps + cas9n_intial_reps + cas9n_final_reps

        # determine starting positions of the normalized reads column for each sample
        # Normalized Library column
        library_norm = total_sample_columns + 1

        # Normalized Cas9+ Initial start column
        i_p_norm_start = library_norm + 1  # library is only one column, so next column is Cas9+ initial

        # Normalized Cas9+ Final start column
        f_p_norm_start = i_p_norm_start + self.num_replicates

        # Normalized Cas9- Initial start column
        i_n_norm_start = f_p_norm_start + self.num_replicates

        # Normalized Cas9- Final start column
        f_n_norm_start = i_n_norm_start + self.num_replicates

        print("\nnorm library column starts at %i, should be 13" % library_norm)
        print("\nnorm initial Cas9+ column starts at %i, should be 14" % i_p_norm_start)
        print("\nnorm final Cas9+ column starts at %i, should be 17" % f_p_norm_start)
        print("\nnorm initial Cas9- column starts at %i, should be 20" % i_n_norm_start)
        print("\nnorm final Cas9- column starts at %i, should be 23" % f_n_norm_start)

        self.column_names = self.raw_counts.columns

        print("\nCalculating log2 fold changes...")

        # call function to calculate log2FC for library-Cas9-initial
        self.liblog2FC_calculator(initial=library_norm, final=i_n_norm_start)

        # call function to calculate log2FC for initial condition
        self.log2FC_calculator(initial=i_n_norm_start, final=i_p_norm_start)

        # call function to calculate log2FC for final condition
        self.log2FC_calculator(initial=f_n_norm_start, final=f_p_norm_start)

        # for debugging
        print(self.raw_counts)

        # # # GENE SCORE # # #

        print("\nCalculating gene scores...")

        # determine number of sgRNAs for each gene
        self.num_sgRNAs = pd.DataFrame()
        self.num_sgRNAs['num sgRNA'] = self.raw_counts.groupby(self.raw_counts.index)[self.column_names[0]].count()

        # determine where log2FC columns begin in raw_counts table
        # columns = int(len(self.raw_counts.columns)) - (len(global_vars.LIST_INITIAL_PATHS)*2)
        # start from library log2FC until Cas9n log2FC
        start_col = f_n_norm_start + self.num_replicates
        end_col = start_col + 3*self.num_replicates - 1

        print("\nStart column for average is %i" % start_col)
        print("\nEnd column for average is %i" % end_col)

        print("\n length of rawcounts %i " % len(self.raw_counts.columns))

        # # # DEBUGGING save to file
        print("\nSaving to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.raw_counts.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + "rawcounts.csv", index=True, sep='\t',
                       float_format='%.16f')

        print("\nDetermining average...")

        # call function to calculate average of log2FC for every gene for each log2FC sample column
        self.avg = pd.DataFrame()  # new dataframe to hold averages
        self.average_calculator(start=start_col, end=end_col)

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

        # insert sgRNA number into ES dataframe
        self.ES['Num.sgRNA'] = self.num_sgRNAs['num sgRNA']

        # self.get_enrichment_score()
        # print(self.ES)

        # new ES method
        # call function to get Library.ES
        self.get_ES(start=3*self.num_replicates)

        # call function to get Initial.ES
        self.get_ES(start=4 * self.num_replicates)

        # call function to get Final.ES
        self.get_ES(start=5 * self.num_replicates)

        # # # DEBUGGING save to file
        print("\nSaving to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.ES.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".ES.csv", index=True,
                               sep='\t',
                               float_format='%.6f')

        print("\nCalculating enrichment ratio...")
        # determine enrichment ratio (ER) -- log2 ratio of ES ( log2(final ES / initial ES)
        self.ES['Enrichment Ratio'] = np.log2(self.ES[self.ES.columns[3]].astype(float) / self.ES[self.ES.columns[2]].astype(float))
        self.ES['Enrichment Ratio'] = self.ES['Enrichment Ratio'].map(lambda x: '%.6f' % x)

        # save to file
        print("\nSaving to file ...")
        os.chdir(global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
        self.ES.to_csv(global_vars.EXPERIMENT_SETTINGS['Experiment name'] + ".csv", index=True, sep='\t', float_format='%.6f')

        print("\nAll done!")

    def liblog2FC_calculator(self, initial, final):

        for col in range(0, self.num_replicates):
            current_f_col = col + final
            new_col_name = self.column_names[current_f_col] + ".Lib.log2FC"

            self.raw_counts[new_col_name] = np.log2(self.raw_counts[self.raw_counts.columns[current_f_col]].astype(float) / self.raw_counts[self.raw_counts.columns[initial]].astype(float))

    def log2FC_calculator(self, initial, final):

        for col in range(0, self.num_replicates):
            current_f_col = col + final
            current_i_col = col + initial
            new_col_name = self.column_names[current_f_col] + ".log2FC"

            self.raw_counts[new_col_name] = np.log2(self.raw_counts[self.raw_counts.columns[current_f_col]].astype(float) / self.raw_counts[self.raw_counts.columns[current_i_col]].astype(float))

    def average_calculator(self, start, end):

        column_names = self.raw_counts.columns

        for col in range(0, 3*self.num_replicates):
            current_col = col + start
            print("\n current current_col = %i" % current_col)
            current_name = column_names[current_col]
            new_col_name = current_name[:-12] + ".avg"  # remove 12 characters from the end (removes ".norm.log2FC")

            # average log2FC from every sgRNA per gene
            self.avg[new_col_name] = self.raw_counts.groupby(self.raw_counts.index)[current_name].mean()

    def get_sgw(self):

        column_names = self.avg.columns

        for col in range(0, len(self.avg.columns)):
            current_name = column_names[col]
            new_col_name = current_name[:-4] + ".sgw"  # remove 4 characters from the end (removes ".avg")

            # determine sgw: sgRNA weighted log2FC (log2FC * # of sgRNAs)
            self.sgw[new_col_name] = self.avg[current_name].astype(float) * self.num_sgRNAs['num sgRNA'].astype(float)

    def get_rank(self):

        column_names = self.sgw.columns

        for col in range(0, len(self.sgw.columns)):
            current_name = column_names[col]
            new_col_name = current_name[:-4] + ".rank"  # remove 4 characters from the end (removes ".sgw")

            # sort values by ascending order in new_col_name
            self.sgw = self.sgw.sort_values(current_name)

            # determine rank of gene (Enrichment Score)
            # increment the rank from 1 to x
            self.sgw[new_col_name] = pd.RangeIndex(stop=self.sgw.shape[0]) + 1

    def get_ES(self, start):

        column_names = self.sgw.columns
        new_col_name = column_names[start]
        new_col_name = new_col_name[:-5] + ".ES"

        self.ES[new_col_name] = self.sgw.iloc[:, int(start):int(start+self.num_replicates)].mean(axis=1) / self.num_sgRNAs['num sgRNA'].astype(float)

    def get_enrichment_score(self):     # update to any number of replicates

        column_names = self.sgw.columns
        i_start_column = int(len(self.sgw.columns)) / 2
        f_start_column = i_start_column + (i_start_column / 2)

        initial_name = column_names[i_start_column]
        initial_name = initial_name[:-7] + ".ES"

        final_name = column_names[f_start_column]
        final_name = final_name[:-7] + ".ES"

        # determine enrichment scores
        self.ES[global_vars.EXPERIMENT_SETTINGS['Initial condition name']] = self.sgw.iloc[:, int(i_start_column):int(f_start_column)].sum(axis=1) / self.num_sgRNAs['num sgRNA'].astype(float)
        self.ES[global_vars.EXPERIMENT_SETTINGS['Final condition name']] = self.sgw.iloc[:, int(f_start_column):int(len(self.sgw.columns))].sum(axis=1) / self.num_sgRNAs['num sgRNA'].astype(float)

        # determine enrichment ratio (ER) -- log2 ratio of ES ( log2(final ES / initial ES)
        self.ES['Enrichment Ratio'] = np.log2(self.ES[global_vars.EXPERIMENT_SETTINGS['Final condition name']] / self.ES[global_vars.EXPERIMENT_SETTINGS['Initial condition name']])
        self.ES['Enrichment Ratio'] = self.ES['Enrichment Ratio'].map(lambda x: '%.6f' % x)


if __name__ == "__main__":
    app = TRACSApp()
    app.geometry("835x600+300+300")
    # app.resizable(False, False)
    app.resizable(True, True)
    app.mainloop()



