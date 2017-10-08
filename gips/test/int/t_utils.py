"""Integration tests for code in gips.utils."""

import os

import pytest

from gips import utils


class TestException(Exception):
    """Custom exception to avoid masking any real one below."""


test_content = "Have you ever retired a human by mistake?"

@pytest.mark.parametrize('should_raise', (False, True))
def t_make_temp_dir(should_raise):
    """Test basic functionality of utils.make_temp_dir().

    Confirm a file in the temp dir is usable and the dir is deleted even
    if an exception is raised.
    """
    # make a temp dir and write to a file in it
    did_raise = False
    try:
        with utils.make_temp_dir() as temp_dir:
            fname = temp_dir + '/testfile'
            with open(fname, 'a') as f:
                f.write(test_content)
            with open(fname, 'r') as f:
                read_content = f.read()
            if should_raise:
                raise TestException('Like our owl?')
    except TestException:
        did_raise = True

    # confirm writing worked, and the tempdir is destroyed on context exit, and an exception,
    # if raised, was permitted to bubble up
    assert (test_content == read_content and not os.path.exists(temp_dir)
            and (did_raise if should_raise else True))
