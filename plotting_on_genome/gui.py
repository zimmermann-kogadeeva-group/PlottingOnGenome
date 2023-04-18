#!/usr/bin/env python

import PySimpleGUI as sg

from .main import run_pipeline


def main():
    # All the stuff inside your window.
    names = [
        [sg.Text("Search term", pad=(5, 7))],
        [sg.Text("Email", pad=(5, 7))],
        [sg.Text("Sequences file", pad=(5, 10))],
        [sg.Text("Output folder")],
    ]

    inputs = [
        [sg.InputText(key="search_term", pad=(5, 10))],
        [sg.InputText(key="email")],
        [sg.Input(), sg.FilesBrowse(key="seq_file")],
        [sg.Input(), sg.FolderBrowse(key="output_prefix")],
    ]

    layout = [
        [sg.Column(names), sg.Column(inputs)],
        [sg.Button("OK", key="OK"), sg.Button("Cancel")],
    ]

    # Create the Window
    window = sg.Window("PlottingOnGenome", layout, element_padding=(5, 5))
    # Event Loop to process "events" and get the "vals" of the inputs
    while True:
        event, vals = window.read()
        # if user closes window or clicks cancel
        if event in (sg.WIN_CLOSED, "Cancel"):
            break
        if event == "OK":
            run_pipeline(
                vals["seq_file"],
                vals["search_term"],
                vals["email"],
                vals["output_prefix"],
            )

            break
    window.close()


if __name__ == "__main__":
    main()
