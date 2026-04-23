"""Parser interfaces and parser plugin registration."""

from quantum_macaroni.parsers.base import (
    ElectronicStructureParser,
    ParserResult,
    available_parsers,
    get_parser,
    register_parser,
)
from quantum_macaroni.parsers.fleur_outxml import FleurOutxmlParser

DEFAULT_PARSER = FleurOutxmlParser()
register_parser(DEFAULT_PARSER)

__all__ = [
    "ElectronicStructureParser",
    "ParserResult",
    "FleurOutxmlParser",
    "DEFAULT_PARSER",
    "register_parser",
    "get_parser",
    "available_parsers",
]
