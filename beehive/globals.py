
import logging

from dataclasses import dataclass, field
from pathlib import Path
from scanpy import AnnData


lg = logging.getLogger(__name__)

@dataclass
class GLOBALS:
    adata: AnnData = field(default_factory=AnnData)
    h5adfile: Path = field(default_factory=Path)
    init: bool = False

    
globals = GLOBALS()
