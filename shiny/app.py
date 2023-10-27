
from pathlib import Path
from typing import Dict

import duckdb
from termite.vis import data

from shiny import App, reactive, render, ui


def get_databases() -> Dict[str,str]:
    """Return a list of databases."""
    return {
        str(x): x.name
        for x in Path('/Users/u0089478/data/termite').glob('*.db')
        if '8' not in x.name
        }


def get_experiments(db) -> Dict[int, str]:
    """Return a dictionary of experiment: descriptions."""
    return (data.get_experiments(db=db)
            .set_index('experiment_id')['full']
            .to_dict())


def get_datasets(expid: int, db) -> Dict[int, str]:
    if not str(expid).strip():
        return {}

    return (data.get_datasets(expid, db=db)
            .set_index('dataset_id')['full']
            .to_dict())


app_ui = ui.page_fluid(
    ui.navset_tab_card(
        ui.nav(
            "Table Builder",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_selectize(
                        "db", "Database", get_databases(),
                        ),
                    ui.input_selectize(
                        "exp", "Experiment", {},
                        ),
                    ui.input_selectize(
                        "dset", "Dataset", {},
                        ),
                ),
                ui.panel_main(
                    ui.h2("??")
                ),
            ),
        ),
        ui.nav("Experiments", "Experiments"),
        ui.nav("Test",
               [
                   ui.h2("Hello Shiny!"),
                   ui.input_slider("n", "N", 0, 100, 20),
                   ui.output_text_verbatim("txt")
                   ]),


        # create gap ----
        ui.nav_spacer(),

        # right hand side ----
        ui.nav_control(
            ui.a("Python", href="https://python.org", target="_blank")
        ),
    ),
)


def server(input, output, session):

    @reactive.Calc
    def db() -> duckdb.DuckDBPyConnection:
        return duckdb.connect(input.db(), read_only=True)

    @output
    @render.text
    def txt() -> str:
        return f"n*2 is {input.n() * 2}"


    @reactive.Effect()
    def _exp() -> None:
        xp = get_experiments(db=db())
        ui.update_select("exp", choices=xp)

    @reactive.Effect()
    def _ds() -> None:
        exp_id = input.exp()
        ds = get_datasets(exp_id, db=db())
        ui.update_select("dset", choices=ds)

app = App(app_ui, server)
