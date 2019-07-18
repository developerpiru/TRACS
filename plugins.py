from tkinter import *
import os
import global_vars
import subprocess, sys
import ntpath


def test_func(frame_name, target):
    tar_obj = getattr(frame_name, target)
    tar_obj.insert(END, "\nthis is from the test func")
    print("This is from the function")


def split_filename(fullpath):
    head, tail = ntpath.split(fullpath)
    return tail or ntpath.basename(head)


