GIPS Automated Testing
======================
See the end of this file for a quick-start section, after configuring pytest:

Configuration
-------------
Create a file named `pytest.ini` in the project directory, and use it to set
the `data-repo` and `output-dir` settings to full paths to the respective
directories.

```
[pytest]
# GIPS' configured data repository:
data-repo = /home/your-user-name-here/src/gips/data-repo
# a directory of your choice, used for output from, e.g., gips_project
output-dir = /home/your-user-name-here/src/gips/testout
# These should remain as they are (placing them in a config file is a TODO):
python_functions =  t_
python_classes =    T_
python_files =      t_*.py
```

Library Dependencies
--------------------
GIPS automated testing uses a few libraries, mainly pytest and mock, and a
library for gluing them together (pytest-mock).  These are installed by
`install.sh`.  Further reading:

* Pytest:  http://docs.pytest.org/en/latest/index.html
* Mock:  http://www.voidspace.org.uk/python/mock/
* pytest-mock:  https://pypi.python.org/pypi/pytest-mock

Test selection
--------------
Running tests is straightforward, but deselect the system tests to save time:

```
py.test 	  	# runs all tests, even slow system tests
py.test -k unit # only runs super-fast unit tests
============================= test session starts ==============================
platform linux2 -- Python 2.7.11+, pytest-2.9.2, py-1.4.31, pluggy-0.3.1
rootdir: /home/tolson/src/gips, inifile: pytest.ini
plugins: mock-1.1, cov-2.2.1, timeout-1.0.0
collected 23 items

gips/test/unit/t_core.py .
gips/test/unit/t_modis_fetch.py ...........

======================= 11 tests deselected by '-kunit' ========================
=================== 12 passed, 11 deselected in 0.21 seconds ===================
```

Select specific tests per usual for pytest:

https://pytest.org/latest/usage.html#specifying-tests-selecting-tests

You can also select tests based on mark (marks are tags for tests):

https://pytest.org/latest/example/markers.html

And reading about fixtures is especially recommended as the GIPS test suite
uses them substantially:

http://pytest.org/latest/fixture.html

**Important caveat:** If you specify the same test file multiple times on the
command line (say, to specify multiple tests in that file), any module-scoped
test fixtures will run once for each time the file is listed, which may not be
what was expected.  This is because pytest sees each item on the command line
as a separate run of the module, even though really it's the same file listed
mulitple times.  Instead, avoid the problem by using `-k` to specify tests.

Writing Unit Tests with Mocking
===============================
Automated testing can be difficult due to the code's interaction with the
outside.  Some interactions may be undesirable for testing purposes, such as
access to an RNG or the network.  Others may be important to observe for
testing purposes, but it may be difficult to do so.  For instance, if a
function makes a library call, it may be useful to know what the arguments
were, but difficult to measure without special tools.

One such special tool is python's mock library.  `t_auth_settings` in
`gips/test/unit/t_modis_fetch.py` is a good example:  Due to `fetch_mocks`, no
real I/O is permitted, and instead, every outside interaction in
`modis.modisAsset.fetch` is carefully monitored.  The goal of the test is to
confirm that `fetch` only asks for auth credentials when circumstances warrant,
so calls to `get_setting` are intercepted and checked accordingly.

Mocking is a complex subject; read about it here:

http://www.voidspace.org.uk/python/mock/

Running System Tests
====================
Each test runs `gips_something` commands, usually just one, as subprocesses.
The test harness captures the command's exit status, standard output, and
standard error.  In addition, the filesystem is observed before and after the
test, and any created, deleted, or updated files are noted.  These observations
are then compared with known correct values to see if the test ought to pass.
Finally any cleanup actions occur, which may include removing created files to
leave the filesystem in a pristine state.

System tests require a little setup, unfortunately:  Edit `settings.py` and set
`OVERRIDE_VERSION` to the version the system tests expect.  You can check this
by looking in `gips/test/sys/expected/modis.py` and checking on the version
that is expected to be output.  You should also have a correct `pytest.ini`
file (see the top-level README.md).  Back up any data files you want to retain
in your data repo directory; some system tests need to operate destructively on
it.

It's not normal to retain state in between test runs, but the pattern of GIPS
operations require it for efficiency, mostly downloading data files via
`gips_inventory --fetch`.  Once a certain set of these files are in place, they
are not downloaded again to improve performance.  Use options to control this
behavior:

```
py.test                 # run tests without fetching data files first; may fail
py.test --setup-repo    # fetch data files then run tests
py.test                 # now this will work since data files are retained
py.test --clear-repo    # delete the entire repo, then fetch data files,
                        # then run tests
```

Also some specially-marked tests are skipped unless command-line options are
given; this is for performance and safety reasons.  The repo-altering tests,
mostly tests of `gips_inventory --fetch`, will always be skipped without the
`--src-altering` option.  Some other tests may not run unless `--slow` is
given.  Note also that the source-altering tests may leave the repo in a state
that requires `--setup-repo` or `--clear-repo` for the other tests to pass.

Writing System Tests:  Testing a New Driver
===========================================
Each driver's tests are split in half; here are the files for modis:

* `gips/test/sys/t_modis.py`:  The tests themselves
* `gips/test/sys/expected/modis.py`:  A file of known good outcomes, one per
  test.

As a quick demonstration, edit a new test file, `gips/test/sys/t_$DRIVER.py`,
and add a test to it:

```
from .util import *
def t_example(repo_env, expected):
    actual = repo_env.run('echo "hello world!"')
    assert expected == actual
```

And an output expectation to a new expectations file,
`gips/test/sys/expected/$DRIVER.py`:

```
t_example = { 'stdout': """hello world!\n""" }
```

Now run the test (and only this test, to save time), and you should observe it
pass:

```
$ py.test -k t_example
============================= test session starts ==============================
platform linux2 -- Python 2.7.11+, pytest-2.9.2, py-1.4.31, pluggy-0.3.1
rootdir: /home/tolson/src/gips, inifile: pytest.ini
plugins: mock-1.1, cov-2.2.1, timeout-1.0.0
collected 23 items

gips/test/sys/t_$DRIVER.py .

===================== 22 tests deselected by '-kt_example' =====================
=================== 1 passed, 22 deselected in 0.42 seconds ====================
```

In practice, copying the patterns in `t_modis.py` for testing a new driver can
be considered adequate.  Also examine any deprecated shell script in
`gips/test`, as these are organized by driver and can supply needed test cases.
For a driver's source-altering tests, edit a single pair of files regardless of
driver:

* `sys/t_repo_src_alteration.py`: For the tests
* `sys/expected/repo_src_alteration.py`: For the expectations.

Working with Environment Fixtures and the Expectations File
-----------------------------------------------------------
For GIPS, system tests means observing command-line processes as they execute
and comparing their output to known-good values.  These outputs only take two
forms:  Standard streams and altering files in the filesystem (including
creating and deleting files).  Some tests are expected to modify the data repo,
while others will modify a project directory.  Choose an appropriate pytest
fixture accordingly, using existing driver tests as a guide.  Pytest fixtures
are available to make this easier; see `gips/test/sys/util.py` for details.

The `expected` fixture is somewhat special:  It provides a convenient way to
store bulky data needed for tests in a separate file, without needing to
explicitly access the file itself nor configure anything.  See
`gips/test/sys/expected/modis.py` for the way this works in practice:  The
`expected` fixture looks for an object with the same name as the currently
running test.  It then loads that dictionary into a `GipsProcResult` object,
which is convenient for making assertions in these system tests.

The known-good values stored in `expected/` are usually captured from initial
test runs.  For convenience, newly developed tests may be run with
`py.test -s --log-level debug`; parts of the output can be cutpasted into an
`expected/` file to establish these needed known-good values.

Caveats & Tips
==============
If you want to view logger output or `print` statements, or if you want to run
the debugger (`import pdb; pdb.set_trace()`) you _must_ specify `py.test -s`.
This is due to pytest's somewhat draconian defaults regarding the brevity of
test run output by default.

Quick Start
===========
```
py.test -k unit # to run unit tests (after activating your virtualenv):
py.test         # run both the system & unit tests (be prepared to wait; note
                # also that additional configuration may be necessary).
```
