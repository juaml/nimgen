"""Provide tests for nimgen/web.py."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import tempfile

from nimgen import web


def test_run_webgestalt():
    """Test run_webgestalt."""

    with tempfile.TemporaryDirectory() as tmp:
        genes = os.path.join(tmp, "genes.txt")
        os.system(f"touch {genes}")
        web.run_webgestalt(genes)
