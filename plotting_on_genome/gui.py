#!/usr/bin/env python

import json
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import PySimpleGUI as sg
from threading import Thread

from .main import Pipeline


matplotlib.use("agg")
settings_cache = Path.home() / ".plotting_on_genome_settings.json"


class PropagatingThread(Thread):
    def run(self):
        self.exc = None
        try:
            if hasattr(self, "_Thread__target"):
                # Thread uses name mangling prior to Python 3.
                self.ret = self._Thread__target(
                    *self._Thread__args, **self._Thread__kwargs
                )
            else:
                self.ret = self._target(*self._args, **self._kwargs)
        except BaseException as e:
            self.exc = e

    def join(self, timeout=None):
        super(PropagatingThread, self).join(timeout)
        if self.exc:
            raise self.exc
        return self.ret


def load_settings():
    # Default settings in case of first use of the package
    settings = {}
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
        [sg.Text("Filter threshold", pad=(5, 7))],
        [sg.Text("Retmax (Max. # of records from NCBI)", pad=(5, 7))],
        [sg.Text("Sequences file", pad=(5, 10))],
        [sg.Text("Output folder")],
    ]

    inputs = [
        [
            sg.InputText(
                key="search_term", default_text=settings.get("search_term"), pad=(5, 7)
            )
        ],
        [sg.InputText(key="email", default_text=settings.get("email"), pad=(5, 7))],
        [
            sg.InputText(
                key="fwd_suffix",
                default_text=settings.get("fwd_suffix", "_F"),
                pad=(5, 7),
            )
        ],
        [
            sg.InputText(
                key="rev_suffix",
                default_text=settings.get("rev_suffix", "_R"),
                pad=(5, 7),
            )
        ],
        [
            sg.InputText(
                key="filter", default_text=settings.get("filter", 0.5), pad=(5, 7)
            )
        ],
        [
            sg.InputText(
                key="retmax", default_text=settings.get("retmax", 200), pad=(5, 7)
            )
        ],
        [
            sg.Input(default_text=settings.get("seq_file")),
            sg.FilesBrowse(key="seq_file"),
        ],
        [
            sg.Input(default_text=settings.get("output_prefix")),
            sg.FolderBrowse(key="output_prefix"),
        ],
    ]

    layout = [
        [sg.Column(names), sg.Column(inputs)],
        [sg.Button("OK", key="OK"), sg.Button("Cancel")],
    ]
    return layout


def save_figures(user_input):
    # Create a new Pipeline object
    pipeline = Pipeline(
        user_input["seq_file"],
        user_input["search_term"],
        user_input["email"],
        user_input["output_prefix"],
        retmax=int(user_input["retmax"]),
        fwd_suffix=user_input["fwd_suffix"],
        rev_suffix=user_input["rev_suffix"],
    )

    # Create linear plots of inserts with annotations
    for seq_id in pipeline.seq_ids:
        for i, insert in enumerate(
            pipeline.get_inserts(
                seq_id, output="matched", filter_threshold=user_input["filter"]
            )
        ):
            fig, axs = plt.subplots(2, 1, figsize=(10, 8))
            pipeline.plot_insert(insert, axs=axs)
            fig.savefig(pipeline.work_dir / f"{seq_id}_hit{i}.png")
            plt.close()

    # Create a plot of genome / all contigs as circular plot with inserts
    # layered on top
    fig, ax = plt.subplots(figsize=(10, 30))
    pipeline.plot_all_db_seqs(
        output="matched", filter_threshold=user_input["filter"], ax=ax
    )
    fig.savefig(pipeline.work_dir / "genome_plot.png")

    fig, axs = plt.subplots(1, 3, figsize=(12, 5))
    pipeline.plot_insert_dists(
        output="matched", filter_threshold=user_input["filter"], axs=axs
    )
    fig.savefig(pipeline.work_dir / "insert_length_dist.png")

    pipeline.to_dataframe(
        output="matched", filter_threshold=user_input["filter"]
    ).to_csv(pipeline.work_dir / "inserts.csv", index=False)


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
            user_input["filter"] = float(user_input.get("filter"))
            # Save the settings for next use
            save_settings(user_input)
            try:
                thread = PropagatingThread(target=save_figures, args=[user_input])
                thread.start()
                sg.popup_non_blocking("Please wait ...")
                thread.join()
            except Exception as e:
                sg.Popup(e)
            break

    window.close()


if __name__ == "__main__":
    main()
