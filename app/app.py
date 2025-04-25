#!/usr/bin/env python3
import argparse

import streamlit as st
from page1 import get_main_inputs, run_pipeline
from page2 import show_results

# st.set_page_config(layout="wide")

if "stage" not in st.session_state:
    st.session_state.stage = 0

if "search_term_count" not in st.session_state:
    st.session_state.search_term_count = 1


def submit(*args):
    with st.spinner("Processing..."):
        res = {name: run_pipeline(**param_set) for name, param_set in args[0].items()}
        if res is not None:
            st.session_state.results = res
            st.session_state.stage = 1


def main():
    parser = argparse.ArgumentParser("PlottingOnGenome")
    parser.add_argument("--workdir", action="store_true")
    args = parser.parse_args()

    col1, col2 = st.columns((6, 1), gap="large", vertical_alignment="bottom")
    with col1:
        st.title("PlottingOnGenome")
    with col2:
        if st.session_state.stage == 1:
            if st.button("Reset"):
                st.session_state.stage = 0
                st.session_state.search_term_count = 1
                st.rerun()

    if st.session_state.stage == 0:
        all_inputs = get_main_inputs(args.workdir)
        st.button(
            "Submit", on_click=submit, args=[all_inputs], use_container_width=True
        )

    if st.session_state.stage >= 1:
        show_results()


if __name__ == "__main__":
    main()
