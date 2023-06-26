#!/usr/bin/env python

import json
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import PySimpleGUI as sg

from .main import Pipeline


matplotlib.use("agg")
settings_cache = Path.home() / ".plotting_on_genome_settings.json"


def load_settings():
    # Default settings in case of first use of the package
    settings = {
        "email": None,
        "search_term": None,
        "fwd_suffix": "_F",
        "rev_suffix": "_R",
        "seq_file": None,
        "output_prefix": None,
    }
    if settings_cache.exists():
        with open(settings_cache, "r") as fh:
            settings = json.load(fh)
    return settings


def save_settings(settings):
    # The two line below make sure that values in settings dictionary
    # are correct, as seq_file and output_prefix can be None when
    # previous input gets used due to how InputText and FolderBrowse
    # object interact in PySimpleGUI
    settings["seq_file"] = settings[0]
    settings["output_prefix"] = settings[1]
    with open(settings_cache, "w") as fh:
        json.dump(settings, fh)


def get_layout():
    # Load previous user input
    settings = load_settings()

    # Set up the layout of the GUI
    names = [
        [sg.Text("Search term", pad=(5, 7))],
        [sg.Text("Email", pad=(5, 7))],
        [sg.Text("Forward suffix", pad=(5, 7))],
        [sg.Text("Reverse suffix", pad=(5, 7))],
        [sg.Text("Sequences file", pad=(5, 10))],
        [sg.Text("Output folder")],
    ]

    inputs = [
        [
            sg.InputText(
                key="search_term", default_text=settings["search_term"], pad=(5, 7)
            )
        ],
        [sg.InputText(key="email", default_text=settings["email"], pad=(5, 7))],
        [
            sg.InputText(
                key="fwd_suffix", default_text=settings["fwd_suffix"], pad=(5, 7)
            )
        ],
        [
            sg.InputText(
                key="rev_suffix", default_text=settings["rev_suffix"], pad=(5, 7)
            )
        ],
        [sg.Input(default_text=settings["seq_file"]), sg.FilesBrowse(key="seq_file")],
        [
            sg.Input(default_text=settings["output_prefix"]),
            sg.FolderBrowse(key="output_prefix"),
        ],
    ]

    layout = [
        [sg.Column(names), sg.Column(inputs)],
        [sg.Button("OK", key="OK"), sg.Button("Cancel")],
    ]
    return layout


def save_figures(user_input):
    pipeline = Pipeline(
        user_input["seq_file"],
        user_input["search_term"],
        user_input["email"],
        user_input["output_prefix"],
        fwd_suffix=user_input["fwd_suffix"],
        rev_suffix=user_input["rev_suffix"],
    )

    for seq_id in pipeline.seq_ids:
        for i, insert in enumerate(pipeline.get_inserts(seq_id)):
            fig, axs = plt.subplots(1, 2, figsize=(10, 8))
            pipeline.plot_insert(insert, axs=axs)
            fig.savefig(pipeline.work_dir / f"{seq_id}_hit{i}.png")
            plt.close()

    fig, ax = plt.subplots(figsize=(10, 30))
    pipeline.plot_all_db_seqs(ax=ax)
    fig.savefig(pipeline.work_dir / "genome_plot.png")


def main():
    # Create the Window
    window = sg.Window("PlottingOnGenome", get_layout(), element_padding=(5, 5))

    # Event Loop to process "events" and get the "user_input" of the inputs
    while True:
        event, user_input = window.read()
        # if user closes window or clicks cancel
        if event in (sg.WIN_CLOSED, "Cancel"):
            break
        if event == "OK":
            # Save the settings for next use
            # TODO: popup to let the user know that the package is running
            save_settings(user_input)
            try:
                save_figures(user_input)
            except Exception as e:
                sg.Popup(e)
            break

    window.close()


if __name__ == "__main__":
    main()
