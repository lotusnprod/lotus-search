# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class PlotlyDashKetcher(Component):
    """A PlotlyDashKetcher component.


    Keyword arguments:

    - id (string; optional):
        Unique ID to identify this component in Dash callbacks.

    - buttonLabel (string; optional)

    - molecule (string; optional)"""

    _children_props = []
    _base_nodes = ["children"]
    _namespace = "plotly_dash_ketcher"
    _type = "PlotlyDashKetcher"

    @_explicitize_args
    def __init__(
        self,
        id=Component.UNDEFINED,
        molecule=Component.UNDEFINED,
        buttonLabel=Component.UNDEFINED,
        **kwargs,
    ):
        self._prop_names = ["id", "buttonLabel", "molecule"]
        self._valid_wildcard_attributes = []
        self.available_properties = ["id", "buttonLabel", "molecule"]
        self.available_wildcard_properties = []
        _explicit_args = kwargs.pop("_explicit_args")
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs and excess named props
        args = {k: _locals[k] for k in _explicit_args}

        super().__init__(**args)
