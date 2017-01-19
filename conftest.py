import logging, os, shutil
import envoy

from _pytest.assertion.util import _compare_eq_dict, _diff_text
from gips.test.sys.util import GipsProcResult, set_constants

# configure logging (can config in command-line & config file, see below)
log_format = '%(levelname)-8s %(filename)s:%(lineno)d: %(message)s'
root_logger = logging.getLogger()
root_streamer = logging.StreamHandler()
root_streamer.setFormatter(logging.Formatter(log_format))
root_logger.addHandler(root_streamer)

_log = logging.getLogger(__name__)

# pytest_* functions are hooks automatically detected by pytest
def pytest_addoption(parser):
    """Add custom options & settings to py.test.

    These set up the data repo & configure log level."""
    help_str = ("Set log level to one of:  'debug', 'info', 'warning', "
                "'error', or 'critical'.  Default is 'warning'.")
    parser.addoption("--log-level", dest='log_level', help=help_str)
    parser.addoption("--ll",        dest='log_level', help='Alias for --log-level')

    help_str = "Set up a test data repo & download data for test purposes."
    parser.addoption("--setup-repo", action="store_true", help=help_str)

    help_str = ("Replace the data repo with a new empty one before testing begins.  "
                "Implies --setup-repo.")
    parser.addoption("--clear-repo", action="store_true", help=help_str)

    help_str = ("The directory housing the data repo for testing purposes.  "
                "MUST match GIPS' configured REPOS setting.")
    parser.addini('data-repo', help=help_str)

    parser.addini('output-dir',
                  help="The directory housing output files from test runs.")

    parser.addoption("--slow", action="store_true", help="Do not skip @slow tests.")

    help_str = ("Do not skip @src_altering tests, which are tests that may "
                "alter or remove source data in the repo.")
    parser.addoption("--src-altering", action="store_true", help=help_str)


def pytest_configure(config):
    """Process user config & command-line options."""
    raw_level = config.getoption("log_level")
    level = ('warning' if raw_level is None else raw_level).upper()
    root_logger.setLevel(level)

    dr = str(config.getini('data-repo'))
    if not dr:
        raise ValueError("No value specified for 'data-repo' in pytest.ini")
    else:
        _log.debug("value detected for data-repo: " + dr)

    set_constants(config)

    if config.getoption("slow"):
        _log.debug("--slow detected; will not skip @slow tests")

    if config.getoption("clear_repo"):
        _log.debug("--clear-repo detected; trashing and rebuilding data repo")
        path = config.getini('data-repo')
        os.path.lexists(path) and shutil.rmtree(path)
        setup_data_repo()
        config.option.setup_repo = True
    elif config.getoption("setup_repo"):
        _log.debug("--setup-repo detected; setting up data repo")
        setup_data_repo()


def setup_data_repo():
    """Construct the data repo if it is absent."""
    # confirm the user's done basic config
    # TODO do these checks every run if fast enough
    gcp = envoy.run("gips_config print")
    if gcp.status_code != 0:
        raise RuntimeError("config check via `gips_config print` failed",
                           gcp.std_out, gcp.std_err, gcp)
    _log.debug("`gips_config print` check succeeded.")

    # set up data root if it isn't there already
    gcp = envoy.run("gips_config env")
    if gcp.status_code != 0:
        raise RuntimeError("data root setup via `gips_config env` failed",
                           gcp.std_out, gcp.std_err, gcp)
    _log.debug("`gips_config env` succeeded; data repo (possibly) created")


def pytest_assertrepr_compare(config, op, left, right):
    """When asserting equality between process results, show detailed differences."""
    checks = (op == '==', isinstance(left, GipsProcResult), isinstance(right, GipsProcResult))
    if not all(checks):
        return

    def header_and_indent(header, lines):
        if lines:
            return [header + ':'] + ['  ' + line for line in lines]
        else:
            return [header + ':  matches']

    verbose = config.getoption('verbose')

    output = ['GipsProcResult == GipsProcResult:']

    oper = {True: '==', False: '!='}[left.exit_status == right.exit_status]
    output += ['exit_status:  {} {} {}'.format(left.exit_status, oper, right.exit_status)]

    for s in ('stdout', 'stderr'):
        l_compare = getattr(left, 'compare_' + s)
        r_compare = getattr(right, 'compare_' + s)
        if l_compare and r_compare:
            (l_stream, r_stream) = (getattr(left, s), getattr(right, s))
            # TODO text diff breaks due to terminal control sequences (I
            # think it's escaping them wrong or not at all).
            output += header_and_indent(s, _diff_text(l_stream, r_stream, verbose))
        else:
            output += [s + ':  ignored']

    updated_diff = _compare_eq_dict(left.updated, right.updated, verbose)
    deleted_diff = _compare_eq_dict(left.deleted, right.deleted, verbose)
    created_diff = _compare_eq_dict(left.created, right.created, verbose)
    output += header_and_indent('updated', updated_diff)
    output += header_and_indent('deleted', deleted_diff)
    output += header_and_indent('created', created_diff)

    return output
