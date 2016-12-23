"""Unit tests for data.modis.modis.modisAsset.fetch."""

import sys
import datetime
import urllib2


import pytest
import mock

from ...data.modis import modis

dt = datetime.datetime

# probably test separately; it preumaturely returns: ('MCD12Q1', 'h12v04', dt(2012, 12, 1, 0, 0)),

# everything for NH's tile for Dec 1 2012 that fails:
http_404_params = (
    (('MOD11A2', 'h12v04', dt(2012, 12, 1, 0, 0)),                  # call
        'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.005/2012.12.01'),  # passed to urlopen

    (('MCD43A2', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'http://e4ftl01.cr.usgs.gov/MOTA/MCD43A2.006/2012.12.01'),

    (('MOD09Q1', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'http://e4ftl01.cr.usgs.gov/MOLT/MOD09Q1.005//2012.12.01'),

    (('MCD43A4', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'http://e4ftl01.cr.usgs.gov/MOTA/MCD43A4.006/2012.12.01'),
)


model_404 = (
    '<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 2.0//EN">\n',
    '<html><head>\n',
    '<title>404 Not Found</title>\n',
    '</head><body>\n',
    '<h1>Not Found</h1>\n',
    # shouldn't matter what it says here vvv
    '<p>The requested URL /MOLT/MOD11A2.005/2012.12.01 was not found on this server.</p>\n',
    '</body></html>\n',
)


@pytest.fixture
def fetch_mocks(mocker):
    """Mock out all I/O used in modisAsset.fetch; used for unit tests below."""
    urlopen = mocker.patch.object(modis.urllib, 'urlopen')

    # called in a loop for items in the listing that match criteria
    get = mocker.patch.object(modis.requests, 'get')
    response = get.return_value
    open = mocker.patch.object(modis, 'open')
    file = open.return_value.__enter__.return_value # context manager
    return (urlopen, get, response, open, file)


@pytest.mark.parametrize("call, expected", http_404_params)
def t_no_http_matching_listings(fetch_mocks, call, expected):
    """Unit test for modisAsset.fetch for assets that 404."""
    (urlopen, get, response, open, file) = fetch_mocks

    # give the listing from which names of asset files are extracted
    urlopen.return_value.readlines.return_value = model_404

    modis.modisAsset.fetch(*call)

    # assertions:
    # It should skip the I/O code except for fetching the directory listing
    [f.assert_not_called() for f in (response.iter_content, get, open, file.write)]
    urlopen.return_value.readlines.assert_called_once_with()
    urlopen.assert_called_once_with(expected)

# VERY truncated snippet of an actual listing file
MYD11A1_listing = [
    '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">\n',
    '<html>\n',
    ' <head>\n',
    '  <title>Index of /MOLA/MYD11A1.005/2012.12.01</title>\n',
    ' </head>\n',
    ' <body>\n',
    '<pre><img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MYD11A1.A2012336.h12v03.005.2012341040539.hdf">'
        'MYD11A1.A2012336.h12v03.005.2012341040539.hdf</a>          2012-12-05 22:15  2.4M  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MYD11A1.A2012336.h12v03.005.2012341040539.hdf.xml">'
        'MYD11A1.A2012336.h12v03.005.2012341040539.hdf.xml</a>      2012-12-06 02:22   11K  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MYD11A1.A2012336.h12v04.005.2012341040543.hdf">'
        'MYD11A1.A2012336.h12v04.005.2012341040543.hdf</a>          2012-12-05 22:15  1.2M  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MYD11A1.A2012336.h12v04.005.2012341040543.hdf.xml">'
        'MYD11A1.A2012336.h12v04.005.2012341040543.hdf.xml</a>      2012-12-06 02:22  9.5K  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MYD11A1.A2012336.h12v05.005.2012341040552.hdf">'
        'MYD11A1.A2012336.h12v05.005.2012341040552.hdf</a>          2012-12-05 22:15  382K  \n',
    '<hr></pre>\n',
    '</body></html>\n'
]


MOD11A1_listing = [
    '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2 Final//EN">\n',
    '<html>\n',
    ' <head>\n',
    '  <title>Index of /MOLT/MOD11A1.005/2012.12.01</title>\n',
    ' </head>\n',
    ' <body>\n',
    '<pre><img src="/icons/blank.gif" alt="Icon "> <a href="?C=N;O=D">Name</a>',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MOD11A1.A2012336.h12v02.005.2012339180538.hdf.xml">'
        'MOD11A1.A2012336.h12v02.005.2012339180538.hdf.xml</a>      2012-12-05 02:13   14K  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MOD11A1.A2012336.h12v03.005.2012339180912.hdf">'
        'MOD11A1.A2012336.h12v03.005.2012339180912.hdf</a>          2012-12-04 12:29  3.6M  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MOD11A1.A2012336.h12v03.005.2012339180912.hdf.xml">'
        'MOD11A1.A2012336.h12v03.005.2012339180912.hdf.xml</a>      2012-12-05 02:13   10K  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MOD11A1.A2012336.h12v04.005.2012339180517.hdf">'
        'MOD11A1.A2012336.h12v04.005.2012339180517.hdf</a>          2012-12-04 12:29  1.8M  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MOD11A1.A2012336.h12v04.005.2012339180517.hdf.xml">'
        'MOD11A1.A2012336.h12v04.005.2012339180517.hdf.xml</a>      2012-12-05 02:13  9.0K  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MOD11A1.A2012336.h12v05.005.2012339042544.hdf">'
        'MOD11A1.A2012336.h12v05.005.2012339042544.hdf</a>          2012-12-03 23:44  385K  \n',
    '<img src="/icons/unknown.gif" alt="[   ]">'
        '<a href="MOD11A1.A2012336.h12v05.005.2012339042544.hdf.xml">'
        'MOD11A1.A2012336.h12v05.005.2012339042544.hdf.xml</a>      2012-12-04 02:12  9.4K  \n',
    '<hr></pre>\n',
    '</body></html>\n'
]


# everything for NH's tile for Dec 1 2012 that succeeds, including both inputs and outputs
http_200_params = (
    (('MYD11A1', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'http://e4ftl01.cr.usgs.gov/MOLA/MYD11A1.005/2012.12.01',
        MYD11A1_listing, # HTML index page
        'MYD11A1.A2012336.h12v04.005.2012341040543.hdf'),

    (('MOD11A1', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/2012.12.01',
        MOD11A1_listing, # HTML index page
        'MOD11A1.A2012336.h12v04.005.2012339180517.hdf'),
)


@pytest.mark.parametrize("call, listing_url, listing, asset_fn", http_200_params)
def t_http_matching_listings(mocker, fetch_mocks, call, listing_url, listing, asset_fn):
    """Query http server, extract asset URL, then download it."""
    (urlopen, get, response, open, file) = fetch_mocks

    # rely on the real one because mocking was too painful
    user = modis.modisAsset.Repository.get_setting('username')
    passwd = modis.modisAsset.Repository.get_setting('password')

    # give the listing from which names of asset files are extracted
    urlopen.return_value.readlines.return_value = listing
    content = ("If you think you understand, ",
                    "you don't.  ",
               "If you think you don't understand, ",
                    "you still don't.")
    response.iter_content.return_value = content
    response.status_code = 200

    modis.modisAsset.fetch(*call)

    # These assertions are too involved, possibly, because fetch() does too much
    # listing assertions:  urllib.urlopen(mainurl).readlines()
    urlopen.assert_called_once_with(listing_url)
    urlopen.return_value.readlines.assert_called_once_with()
    # request assertions:  response = request.get(...) && response.iter_content()
    get.assert_called_once_with(listing_url + '/' + asset_fn, auth=(user, passwd), timeout=10)
    response.iter_content.assert_called_once_with()
    # file write assertions:  open(...) as fd && fd.write(...)
    assert open.call_args[0][0].endswith(asset_fn) # did we open the right filename?
    file.write.assert_has_calls([mock.call(c) for c in content])


MOD10A1_listing = ['total 184064\r\n',  # drastically foreshortened ftp listing
    'lrwxrwxrwx  1 7372 90      97 Nov 23  2015 '
        'BROWSE.MOD10A1.A2012336.h12v02.005.2012339213026.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.04/'
        'BROWSE.MOD10A1.A2012336.h12v02.005.2012339213026.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Nov 23  2015 '
        'BROWSE.MOD10A1.A2012336.h12v03.005.2012339212814.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.04/BROWSE.MOD10A1.A2012336.h12v03.005.2012339212814.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Nov 23  2015 '
        'BROWSE.MOD10A1.A2012336.h12v04.005.2012339213007.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.04/BROWSE.MOD10A1.A2012336.h12v04.005.2012339213007.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Nov 23  2015 '
        'BROWSE.MOD10A1.A2012336.h12v05.005.2012339213007.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.04/BROWSE.MOD10A1.A2012336.h12v05.005.2012339213007.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Nov 23  2015 '
        'BROWSE.MOD10A1.A2012336.h12v07.005.2012339212842.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.04/BROWSE.MOD10A1.A2012336.h12v07.005.2012339212842.1.jpg\r\n',
    '-rw-r--r--  1 7372 90    8711 Jan 15  2013 '
        'MOD10A1.A2012336.h12v02.005.2012339213026.hdf.xml\r\n',
    '-rw-r--r--  1 7372 90 1137281 Dec  4  2012 MOD10A1.A2012336.h12v03.005.2012339212814.hdf\r\n',
    '-rw-r--r--  1 7372 90    8710 Jan 15  2013 '
        'MOD10A1.A2012336.h12v03.005.2012339212814.hdf.xml\r\n',
    '-rw-r--r--  1 7372 90  844320 Dec  4  2012 MOD10A1.A2012336.h12v04.005.2012339213007.hdf\r\n',
    '-rw-r--r--  1 7372 90    7935 Jan 15  2013 '
        'MOD10A1.A2012336.h12v04.005.2012339213007.hdf.xml\r\n',
    '-rw-r--r--  1 7372 90  131621 Dec  4  2012 MOD10A1.A2012336.h12v05.005.2012339213007.hdf\r\n',
]


MYD10A1_listing = ['total 184704\r\n', # drastically foreshortened ftp listing
    'lrwxrwxrwx  1 7372 90      97 Dec  5  2015 '
        'BROWSE.MYD10A1.A2012336.h12v02.005.2012340031948.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.05/BROWSE.MYD10A1.A2012336.h12v02.005.2012340031948.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Dec  5  2015 '
        'BROWSE.MYD10A1.A2012336.h12v03.005.2012340032022.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.05/BROWSE.MYD10A1.A2012336.h12v03.005.2012340032022.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Dec  5  2015 '
        'BROWSE.MYD10A1.A2012336.h12v04.005.2012340031954.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.05/BROWSE.MYD10A1.A2012336.h12v04.005.2012340031954.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Dec  5  2015 '
        'BROWSE.MYD10A1.A2012336.h12v05.005.2012340032147.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.05/BROWSE.MYD10A1.A2012336.h12v05.005.2012340032147.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Dec  5  2015 '
        'BROWSE.MYD10A1.A2012336.h12v02.005.2012340031948.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.05/BROWSE.MYD10A1.A2012336.h12v02.005.2012340031948.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Dec  5  2015 '
        'BROWSE.MYD10A1.A2012336.h12v03.005.2012340032022.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.05/BROWSE.MYD10A1.A2012336.h12v03.005.2012340032022.1.jpg\r\n',
    'lrwxrwxrwx  1 7372 90      97 Dec  5  2015 '
        'BROWSE.MYD10A1.A2012336.h12v04.005.2012340031954.1.jpg -> '
        '../../../../DP0/BRWS/Browse.001/2012.12.05/BROWSE.MYD10A1.A2012336.h12v04.005.2012340031954.1.jpg\r\n',
    '-rw-r--r--  1 7372 90    8304 Jan 10  2013 MYD10A1.A2012336.h12v03.005.2012340032022.hdf.xml\r\n',
    '-rw-r--r--  1 7372 90  262727 Dec  4  2012 MYD10A1.A2012336.h12v04.005.2012340031954.hdf\r\n',
    '-rw-r--r--  1 7372 90    7916 Jan 10  2013 MYD10A1.A2012336.h12v04.005.2012340031954.hdf.xml\r\n',
    '-rw-r--r--  1 7372 90  164381 Dec  4  2012 MYD10A1.A2012336.h12v05.005.2012340032147.hdf\r\n',
    '-rw-r--r--  1 7372 90    7920 Jan 10  2013 MYD10A1.A2012336.h12v05.005.2012340032147.hdf.xml\r\n',
]


@pytest.mark.parametrize("call, listing_url, listing, asset_fn", (
    (('MYD10A1', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'ftp://n5eil01u.ecs.nsidc.org/SAN/MOSA/MYD10A1.005/2012.12.01',
        MYD10A1_listing, # FTP index page
        'MYD10A1.A2012336.h12v04.005.2012340031954.hdf'),
    (('MOD10A1', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'ftp://n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A1.005/2012.12.01',
        MOD10A1_listing, # FTP index page
        'MOD10A1.A2012336.h12v04.005.2012339213007.hdf'),
    )
)
def t_ftp_matching_listings(mocker, fetch_mocks, call, listing_url, listing, asset_fn):
    """Query ftp server, extract asset URL, then download it."""
    (urlopen, get, _, open, file) = fetch_mocks

    # give the listing from which names of asset files are extracted
    urlopen.return_value.readlines.return_value = listing
    # content returned by ftp download
    content = ("If you think you understand, you don't.  \n"
               "If you think you don't understand, you still don't.")
    urlopen2 = mocker.patch.object(modis.urllib2, 'urlopen')
    connection = urlopen2.return_value
    connection.read.return_value = content

    modis.modisAsset.fetch(*call)

    # listing assertions:  urllib.urlopen(mainurl).readlines()
    urlopen.assert_called_once_with(listing_url)
    urlopen.return_value.readlines.assert_called_once_with()
    # test the way the url was opened
    get.assert_not_called() # should use urllib2 instead
    urlopen2.assert_called_once_with(listing_url + '/' + asset_fn)
    # file write assertions:  open(...) as fd && fd.write(...)
    assert open.call_args[0][0].endswith(asset_fn) # did we open the right filename?
    file.write.assert_called_once_with(content)


def t_auth_settings(mocker, fetch_mocks):
    """Confirm auth settings are only queried when the asset warrants it."""
    gs_mock = mocker.patch.object(modis.modisAsset.Repository, 'get_setting')

    call = ('MYD11A1', 'h12v04', dt(2012, 12, 1, 0, 0))
    modis.modisAsset.fetch(*call)

    call = ('MOD10A2', 'h12v04', dt(2012, 12, 1, 0, 0))
    modis.modisAsset.fetch(*call)

    assert all([
        gs_mock.call_count == 2,
        gs_mock.mock_calls[0] == (('username',), {}),
        gs_mock.mock_calls[1] == (('password',), {}),
    ])


@pytest.mark.parametrize("call, listing_url, listing, asset_fn", http_200_params)
def t_auth_error(fetch_mocks, mocker, call, listing_url, listing, asset_fn):
    """Confirm auth failures are handled gracefully."""
    (urlopen, get, response, open, file) = fetch_mocks
    print_mock = mocker.patch.object(modis, 'print') # to confirm user notification is produced

    urlopen.return_value.readlines.return_value = listing # coerce the loop to start

    # rig requests.get() to fail
    get.return_value.status_code = 401
    get.return_value.reason = 'Unauthorized'
    get.return_value.text = 'HTTP Basic: Access denied.\n'

    modis.modisAsset.fetch(*call)

    assert all([
        not open.called, # Confirm short-circuiting by showing open was never called
        print_mock.call_count == 1,
        # making sure the errors go to stderr important; precise wording probably isn't
        print_mock.mock_calls[0][-1] == dict(file=sys.stderr),
    ])
