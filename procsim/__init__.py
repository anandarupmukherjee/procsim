"""procsim: a steady-state process simulator for chemical and mineral processing flowsheets."""
from .flowsheet import Flowsheet, solve, solve_json, catalogue_json, FlowsheetError
from .stream import Stream, Context
from .units import UNIT_TYPES
from .components import COMPONENTS, MINERALS

__version__ = "0.3.0"
