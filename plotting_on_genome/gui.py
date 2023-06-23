#!/usr/bin/env python

import json
import matplotlib.pyplot as plt
from pathlib import Path
import PySimpleGUI as sg

from .main import Pipeline


settings_cache = Path.home() / ".plotting_on_genome_settings.json"


def load_settings():
    settings = {
        "email": None,
        "search_term": None,
        "seq_file": None,
        "output_prefix": None,
    }
    if settings_cache.exists():
        with open(settings_cache, "r") as fh:
            settings = json.load(fh)
    return settings


def save_settings(settings):
    settings[0] = settings["seq_file"]
    settings[1] = settings["output_prefix"]
    with open(settings_cache, "w") as fh:
        json.dump(settings, fh)


def main():
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
        [sg.InputText(key="fwd_suffix", default_text="_F", pad=(5, 7))],
        [sg.InputText(key="rev_suffix", default_text="_R", pad=(5, 7))],
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

    # Create the Window
    window = sg.Window("PlottingOnGenome", layout, element_padding=(5, 5))
    # Event Loop to process "events" and get the "vals" of the inputs
    while True:
        event, vals = window.read()
        # if user closes window or clicks cancel
        if event in (sg.WIN_CLOSED, "Cancel"):
            break
        if event == "OK":
            # The two line below make sure that values in settings dictionary
            # are correct, as seq_file and output_prefix can be None when
            # previous input gets used due to how InputText and FolderBrowse
            # object interact in PySimpleGUI
            vals["seq_file"] = vals[0]
            vals["output_prefix"] = vals[1]

            # Save the settings for next use
            # TODO: popup to let the user know that the package is running
            save_settings(vals)
            try:
                pipeline = Pipeline(
                    vals["seq_file"],
                    vals["search_term"],
                    vals["email"],
                    vals["output_prefix"],
                    fwd_suffix=vals["fwd_suffix"],
                    rev_suffix=vals["rev_suffix"],
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

            except Exception as e:
                sg.Popup(e)

            break
    window.close()


if __name__ == "__main__":
    main()
