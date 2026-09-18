"""Parser interfaces and parser plugin registration."""

from quantum_macaroni.parsers.base import (
    ElectronicStructureParser,
    ParserResult,
    available_parsers,
    get_parser,
    register_parser,
)
from quantum_macaroni.parsers.crystal_outp import CrystalOutpParser
from quantum_macaroni.parsers.fleur_outxml import FleurOutxmlParser

DEFAULT_PARSER = FleurOutxmlParser()
register_parser(DEFAULT_PARSER)
register_parser(CrystalOutpParser())

__all__ = [
    "ElectronicStructureParser",
    "ParserResult",
    "FleurOutxmlParser",
    "CrystalOutpParser",
    "DEFAULT_PARSER",
    "register_parser",
    "get_parser",
    "available_parsers",
]
