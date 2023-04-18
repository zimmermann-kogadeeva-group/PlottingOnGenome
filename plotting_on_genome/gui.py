#!/usr/bin/env python

import PySimpleGUI as sg

from .main import run_pipeline


def main():

    # All the stuff inside your window.
    names = [
        [sg.Text("Search term")],
        [sg.Text("Email")]
    ]

    inputs = [
        [sg.InputText(key="search_term")],
        [sg.InputText(key="email")]
    ]

    layout = [
        [sg.Column(names), sg.Column(inputs)], 
        [sg.Text("Sequences file"), sg.Input(), sg.FilesBrowse(key="seq_file")], 
        [sg.Text("Output folder"), sg.Input(), sg.FolderBrowse(key="output_prefix")],
        [sg.Button('OK', key="OK"), sg.Button("Cancel")]
    ]

    # Create the Window
    window = sg.Window('PlottingOnGenome', layout)
    # Event Loop to process "events" and get the "vals" of the inputs
    while True:
        event, vals = window.read()
        # if user closes window or clicks cancel
        if event in (sg.WIN_CLOSED, 'Cancel'):
            break
        if event == "OK":
            run_pipeline(vals["seq_file"], 
                         vals["search_term"], 
                         vals["email"], 
                         vals["output_prefix"])

            break
    window.close()
    

if __name__ == "__main__":
    main()
