#!/usr/bin/env python

import PySimpleGUI as sg

from main import run_pipeline


def main():

    # All the stuff inside your window.
    names = [
        [sg.Text("sequence file")],
        [sg.Text("search term")],
        [sg.Text("email")],
        [sg.Text("output prefix")],
        [sg.Text("images prefix (optional)")]
    ]

    inputs = [
        [sg.InputText(key="seq_file",)],
        [sg.InputText(key="search_term")],
        [sg.InputText(key="email")],
        [sg.InputText(key="output_prefix")],
        [sg.InputText(key="images_prefix")]
    ]

    layout = [
        [sg.Column(names), sg.Column(inputs)], 
        [sg.Button('OK', key="OK"), sg.Button("Cancel")]
    ]

    # Create the Window
    window = sg.Window('PlottingOnGenome', layout)
    # Event Loop to process "events" and get the "vals" of the inputs
    while True:
        event, vals = window.read()
        # if user closes window or clicks cancel
        # print()
        if event in (sg.WIN_CLOSED, 'Cancel'):
            break
        if event == "OK":
            if not vals["images_prefix"]: vals["images_prefix"] = None
            run_pipeline(vals["seq_file"], 
                         vals["search_term"], 
                         vals["email"], 
                         vals["output_prefix"], 
                         vals["images_prefix"])

            break
    window.close()
    

if __name__ == "__main__":
    main()
