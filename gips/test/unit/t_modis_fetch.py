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
        'https://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.006/2012.12.01'), # passed to urlopen

    (('MCD43A2', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'https://e4ftl01.cr.usgs.gov/MOTA/MCD43A2.006/2012.12.01'),

    (('MOD09Q1', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'https://e4ftl01.cr.usgs.gov/MOLT/MOD09Q1.006//2012.12.01'),

    (('MCD43A4', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'https://e4ftl01.cr.usgs.gov/MOTA/MCD43A4.006/2012.12.01'),
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

    # called once to get a listing of files then again when one is chosen to download
    get = mocker.patch.object(modis.requests, 'get')
    listing = mocker.Mock()
    listing.status_code = 200
    content = mocker.Mock()
    content.status_code = 200
    get.side_effect = [listing, content]

    open = mocker.patch.object(modis, 'open')
    file = open.return_value.__enter__.return_value # context manager
    return (get, listing, content, open, file)


@pytest.mark.parametrize("call, expected", http_404_params)
def t_no_http_matching_listings(fetch_mocks, call, expected):
    """Unit test for modisAsset.fetch for assets that 404."""
    # setup & mocks
    (get, listing, content, open, file) = fetch_mocks
    listing.iter_lines.return_value = model_404 # . . . just content that doesn't match

    # call
    modis.modisAsset.fetch(*call)

    # assertions
    assert expected in get.call_args[0]
    get.assert_called_once()
    listing.iter_lines.assert_called_once_with()
    # It should skip the I/O code except for fetching the directory listing
    [f.assert_not_called() for f in (open, file.write)]


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
        'https://e4ftl01.cr.usgs.gov/MOLA/MYD11A1.006/2012.12.01',
        MYD11A1_listing, # HTML index page
        'MYD11A1.A2012336.h12v04.005.2012341040543.hdf'),

    (('MOD11A1', 'h12v04', dt(2012, 12, 1, 0, 0)),
        'https://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.006/2012.12.01',
        MOD11A1_listing, # HTML index page
        'MOD11A1.A2012336.h12v04.005.2012339180517.hdf'),
)


@pytest.mark.parametrize("call, listing_url, listing_html, asset_fn", http_200_params)
def t_http_matching_listings(mocker, fetch_mocks, call, listing_url, listing_html, asset_fn):
    """Query http server, extract asset URL, then download it."""
    #(urlopen, get, response, open, file) = fetch_mocks
    (get, listing, content, open, file) = fetch_mocks

    listing.iter_lines.return_value = listing_html

    content_data = ("If you think you understand, ",
                    "you don't.  ",
               "If you think you don't understand, ",
                    "you still don't.")
    content.iter_content.return_value = content_data

    # rely on the real one because mocking was too painful
    user = modis.modisAsset.Repository.get_setting('username')
    passwd = modis.modisAsset.Repository.get_setting('password')

    modis.modisAsset.fetch(*call)

    # assertions
    assert listing_url in get.call_args_list[0][0] and 2 == get.call_count
    listing.iter_lines.assert_called_once_with()
    # request assertions:  response = request.get(...) && response.iter_content()
    get.assert_called_with(listing_url + '/' + asset_fn, auth=(user, passwd), timeout=30)
    content.iter_content.assert_called_once_with()
    # file write assertions:  open(...) as fd && fd.write(...)
    assert open.call_args[0][0].endswith(asset_fn) # did we open the right filename?
    file.write.assert_has_calls([mock.call(c) for c in content_data])


@pytest.mark.parametrize("call, listing_url, listing_html, asset_fn", http_200_params)
def t_auth_error(fetch_mocks, mocker, call, listing_url, listing_html, asset_fn):
    """Confirm auth failures are handled gracefully."""
    (get, listing, content, open, file) = fetch_mocks
    print_mock = mocker.patch.object(modis, 'print') # to confirm user notification is produced

    listing.iter_lines.return_value = listing_html

    content.status_code = 401 # rig requests.get() to fail for content fetch
    content.reason = 'Unauthorized'
    content.text = 'HTTP Basic: Access denied.\n'

    modis.modisAsset.fetch(*call)

    assert all([
        not open.called, # Confirm short-circuiting by showing open was never called
        print_mock.call_count == 1,
        # making sure the errors go to stderr important; precise wording probably isn't
        print_mock.mock_calls[0][-1] == dict(file=sys.stderr),
    ])
