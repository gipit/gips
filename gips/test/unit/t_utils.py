"""Unit tests for code found in gips.utils."""

import sys
import datetime

import pytest

from gips import utils
from gips import inventory


def t_remove_files(mocker):
    """remove_files should remove the permutations of its arguments."""
    filenames = ['a.hdf', 'b.hdf']
    extensions = ['.index', '.aux.xml']
    removals = ('a.hdf', 'a.hdf.index', 'a.hdf.aux.xml', 'b.hdf', 'b.hdf.index', 'b.hdf.aux.xml')
    m_isfile = mocker.patch.object(utils.os.path, 'isfile')
    m_isfile.return_value = True
    m_os_remove = mocker.patch.object(utils.os, 'remove')
    utils.remove_files(filenames, extensions)
    [m_isfile.assert_any_call(r) for r in removals]
    [m_os_remove.assert_any_call(r) for r in removals]
    assert m_os_remove.call_count == 6


def t_remove_files_no_ext(mocker):
    """remove_files should work correctly when `extensions` is defaulted."""
    filenames = ['a.hdf', 'b.hdf']
    removals = ('a.hdf', 'b.hdf')
    m_isfile = mocker.patch.object(utils.os.path, 'isfile')
    m_isfile.return_value = True
    m_os_remove = mocker.patch.object(utils.os, 'remove')
    utils.remove_files(filenames)
    [m_os_remove.assert_any_call(r) for r in removals]
    [m_isfile.assert_any_call(r) for r in removals]
    assert m_os_remove.call_count == 2


def t_settings_user(mocker):
    """gips.settings should load user settings first."""
    mocker.patch.object(utils.os.path, 'isfile').return_value = True
    mocker.patch.object(utils.os.path, 'expanduser').return_value = 'whatever'
    m_load_source = mocker.patch.object(utils.imp, 'load_source')
    fake_settings = m_load_source.return_value # a MagicMock
    assert utils.settings() == fake_settings


def t_settings_global(mocker):
    """gips.settings should fall back on gips.settings when user settings fail."""
    # force into the second clause
    mocker.patch.object(utils.os.path, 'isfile').return_value = False
    # fake out `import gips.settings` with mocks and trickery:
    fake_gips = mocker.Mock()
    fake_settings = fake_gips.settings
    sys.modules['gips'] = fake_gips
    sys.modules['gips.settings'] = fake_settings
    assert utils.settings() == fake_settings


def t_open_vector_error_handling(mocker):
    m_settings = mocker.patch.object(utils, 'settings')
    m_settings.return_value.DATABASES = {'fakedbsetting': {
        'NAME': 'fake-db',
        'USER': 'fake-user',
        'PASSWORD': 'fake-password',
        'HOST': 'fake-oka',
        'PORT': '5432',
    }}
    m_GeoVector = mocker.patch.object(utils, 'GeoVector')
    fname = 'fakedbsetting:bar'
    utils.open_vector(fname)
    m_GeoVector.return_value.SetPrimaryKey.assert_called_once_with("")


@pytest.mark.parametrize('dt, variant, expected', (
    # standard julian day
    (datetime.datetime(1950, 6, 1, 17, 0), None, 2433434.208333),
    (datetime.datetime(2005, 10, 10, 3, 15), None, 2453653.635417),
    (datetime.datetime(2016, 2, 29, 11, 59), None, 2457447.999306), # leap year
    # modified julian day
    (datetime.datetime(1950, 6, 1, 17, 0), 'modified',   33433.708333),
    (datetime.datetime(2005, 10, 10, 3, 15), 'modified', 53653.135417),
    (datetime.datetime(2016, 2, 29, 11, 59), 'modified', 57447.499306), # leap year
    # CNES julian day (based on 1950/1/1, 0h)
    (datetime.datetime(1950, 6, 1, 17, 0), 'cnes',   151.708333),
    (datetime.datetime(2005, 10, 10, 3, 15), 'cnes', 20371.135417),
    (datetime.datetime(2016, 2, 29, 11, 59), 'cnes', 24165.499306), # leap year
))
def t_julian_date(dt, variant, expected):
    """Test the gloriously pure function, julian_date.

    Expectations established by USNO's website:
    http://aa.usno.navy.mil/data/docs/JulianDate.php
    """
    # rounding chosen to match USNO calculator page ---------v
    assert expected == round(utils.julian_date(dt, variant), 6)


@pytest.yield_fixture
def restore_error_handler():
    handler = utils.error_handler
    yield
    utils.set_error_handler(handler)


def t_set_error_handler(restore_error_handler):
    """Test setting of the GIPS error handler."""
    def fake_handler():
        pass
    utils.set_error_handler(fake_handler)
    assert utils.error_handler is fake_handler


@pytest.mark.parametrize('verbosity, print_exc_call_cnt', (
    (utils._traceback_verbosity - 3, 0),
    (utils._traceback_verbosity + 0, 1),
    (utils._traceback_verbosity + 3, 1)))
def t_report_error(mocker, verbosity, print_exc_call_cnt):
    """Test GIPS' general purpose error reporting."""
    mocker.patch.object(utils.gippy.Options, 'Verbose').return_value = verbosity
    m_print_exc = mocker.patch.object(utils.traceback, 'print_exc')
    m_verbose_out = mocker.patch.object(utils, 'verbose_out')

    utils.report_error(Exception('blarg'), 'error message here')

    assert (m_print_exc.call_count == print_exc_call_cnt and
            m_verbose_out.call_count == 1)


# for t_gips_exit
e = Exception('aaah!')
e.msg_prefix = 'oh noze!'

@pytest.mark.parametrize('accum_errors', ([], [e, e, e]))
def t_gips_exit(mocker, accum_errors):
    """Test GIPS' exit function."""
    m_report_error = mocker.patch.object(utils, 'report_error')
    m_sys_exit = mocker.patch.object(utils.sys, 'exit')
    m_sys_exit.side_effect = RuntimeError('real one raises SystemExit')
    mocker.patch.object(utils, '_accumulated_errors', new=accum_errors)

    with pytest.raises(RuntimeError):
        utils.gips_exit()

    # confirm exit status is 0 if no errors, 1 otherwise:
    m_sys_exit.assert_called_once_with(0 if len(accum_errors) == 0 else 1)
    # confirm an error report is generated for each error:
    assert m_report_error.call_count == len(accum_errors)


@pytest.mark.parametrize('driver_string, stop_on_error, setup_orm', (
    (None,    False, False),
    (None,    True,  False),
    ('modis', False, False),
    (None,    False, True)))
def t_gips_script_setup(mocker, driver_string, stop_on_error, setup_orm):
    """Test setup function for GIPS-as-a-CLI-application."""
    mocker.patch.object(utils, '_stop_on_error', new=False)
    m_idc = mocker.patch.object(utils, 'import_data_class')
    m_gips_inv = mocker.Mock()
    m_gips_inv.orm.setup # to force the child mocks to exist
    mocker.patch.dict('sys.modules', {'gips.inventory': m_gips_inv})

    rv = utils.gips_script_setup(driver_string, stop_on_error, setup_orm)

    if driver_string is None:
        m_idc.assert_not_called()
    else:
        m_idc.assert_called_once_with(driver_string)
    assert (utils._stop_on_error == stop_on_error and
            rv is {True: None, False: m_idc.return_value}[driver_string is None] and
            m_gips_inv.orm.setup.call_count == {False: 0, True: 1}[setup_orm])


def t_lib_error_handler_base_case(mocker):
    """Test error handler for GIPS' library mode:  non-entry into except."""
    m_report_error = mocker.patch.object(utils, 'report_error')

    # implicit assertion is that no exception is raised
    with utils.lib_error_handler():
        pass # stand-in for non-raising code

    m_report_error.assert_not_called()


@pytest.mark.parametrize('continuable, stop_on_error', (
    (False, False),
    (False, True),
    (True,  True)))
def t_lib_error_handler_reraise(mocker, continuable, stop_on_error):
    """Test exception-reraising cases of lib_error_handler."""
    m_report_error = mocker.patch.object(utils, 'report_error')

    mocker.patch.object(utils, '_stop_on_error', new=stop_on_error)
    with pytest.raises(RuntimeError):
        with utils.lib_error_handler('PREFIX', continuable):
            raise RuntimeError("AAAAAH!")
    m_report_error.assert_not_called()


def t_lib_error_handler_continuable_case(mocker):
    """Test error handler for GIPS' library mode."""
    m_report_error = mocker.patch.object(utils, 'report_error')
    mocker.patch.object(utils, '_stop_on_error', new=False)

    rte = RuntimeError("AAAAAH!")
    with utils.lib_error_handler('PREFIX', True):
        raise rte
    m_report_error.assert_called_once_with(rte, 'PREFIX')


def t_cli_error_handler_base_case(mocker):
    """Test error handler for GIPS' CLI mode:  non-entry into except."""
    m_report_error = mocker.patch.object(utils, 'report_error')

    # implicit assertion is that no exception is raised
    with utils.cli_error_handler():
        pass # stand-in for non-raising code

    m_report_error.assert_not_called()

@pytest.mark.parametrize('continuable, stop_on_error', (
    (False, False),
    (False, True),
    (True,  True)))
def t_cli_error_handler_reraise(mocker, continuable, stop_on_error):
    """Test exception-reraising cases of cli_error_handler."""
    m_report_error = mocker.patch.object(utils, 'report_error')
    m_gips_exit = mocker.patch.object(utils, 'gips_exit')

    mocker.patch.object(utils, '_stop_on_error', new=stop_on_error)
    with utils.cli_error_handler('PREFIX', continuable):
        raise RuntimeError("AAAAAH!")
    assert 0 == m_report_error.call_count and 1 == m_gips_exit.call_count


def t_cli_error_handler_continuable_case(mocker):
    """Test error handler for GIPS' CLI mode."""
    m_report_error = mocker.patch.object(utils, 'report_error')
    mocker.patch.object(utils, '_stop_on_error', new=False)

    rte = RuntimeError("AAAAAH!")
    with utils.cli_error_handler('PREFIX', True):
        raise rte
    m_report_error.assert_called_once_with(rte, 'PREFIX')
