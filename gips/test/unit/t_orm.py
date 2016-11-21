import os, fnmatch, datetime

import pytest

from gips.inventory import orm

"""Test inventory.orm.* for correctness."""

@pytest.mark.parametrize("env_var, expected", (
    ('true',        True),
    ('TrUe',        True),
    ('TRUE',        True),
    ('false',       False),
    ('99agmavfa95', False),
    ('1',           True),
    ('1.0',         True),
    ('1.1',         True),
    ('-1.0',        True),
    ('-1.1',        True),
    ('0',           False),
    ('-0',          False),
    ('0.0',         False),
    ('',            False),
    (' ',           False),
))
def t_use_orm(mocker, env_var, expected):
    """Test various inputs for the env var."""
    mocker.patch.dict('gips.inventory.orm.os.environ', GIPS_ORM=env_var)
    assert expected == orm.use_orm()

def t_use_orm_unset_var(mocker):
    """Unset env var should result in ORM use."""
    mocker.patch.dict('gips.inventory.orm.os.environ', clear=True)
    assert orm.use_orm()

@pytest.mark.parametrize('setup_complete, use_orm, expected', (
    (False, True,  True),
    (False, False, False),
    (True,  True,  False),
    (True,  False, False),
))
def t_setup(mocker, setup_complete, use_orm, expected):
    """Test orm.setup() to confirm it only runs setup at the right times."""
    mocker.patch('gips.inventory.orm.setup_complete', setup_complete)
    m_use_orm = mocker.patch('gips.inventory.orm.use_orm')
    m_use_orm.return_value = use_orm
    m_django_setup = mocker.patch('gips.inventory.orm.django.setup')
    orm.setup()
    # confirm it ran/didn't run django.setup as expected, and confirm 
    # that the global var got set in any case
    assert expected == m_django_setup.called and orm.setup_complete
