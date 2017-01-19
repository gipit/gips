"""Unit tests for core functions, such as those found in data.core."""

import sys

import pytest
import mock

import gips

def t_version_override(mocker):
    """Test gips.__init__.detect_version() for correct override of __version__."""
    env = mocker.patch.object(gips.os, 'environ')
    # os.environ.get is called by libs as well as detect_version(); fortunately no harm seems to
    # come from giving them bad results.

    # no override requested
    env.get.side_effect = lambda key, default=None: default # key not found
    version_a = gips.detect_version()

    # override requested
    env.get.side_effect = lambda key, default=None: 'fancy-new-version'
    version_b = gips.detect_version()

    env.get.assert_has_calls([ # assert two identical calls
        mock.call('GIPS_OVERRIDE_VERSION', gips.version.__version__) for _ in range(2)
    ])
    assert (version_a, version_b) == (gips.version.__version__, 'fancy-new-version')
