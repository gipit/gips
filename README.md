# GIPS

See http://gipit.github.io/gips/ for documentation, but know it is not
necessarily current.

## Installation

After cloning this git repo & changing to its directory, run `install.sh`,
which only officially supports recent versions of Ubuntu.  It will use `sudo`
to install system packages, and may ask for authentication accordingly.  It
runs apt-get, which may prompt you for confirmation of its actions.  
Finally, the script will show you how to set GIPS system settings with `gips_config`.  

### MODIS configuration note
If you wish to work with modis data, edit `~/.gips/settings.py` (if 
that file doesn't exist, run `gips_config user` to create it) and
alter the `'modis'` entry in the `REPOS` dict by adding the following
to the end of that file:

```
REPOS['modis'].update(
    {
        'username': 'YOUR-EARTHDATA-USERNAME',
        'password': 'YOUR-EARTHDATA-PASSWORD', # Use os.environ to avoid saving
                                           # a password in a file.
    }
)
```

This username and password must match your Earthdata credentials; for more
information on obtaining credentials:

https://lpdaac.usgs.gov/about/news_archive/important_changes_lp_daac_data_access

For information on setting up automated testing, see `gips/test/README.md`.

## Authors and Contributors
The following have been authors or contributers to GIPS

    Matthew Hanson, matt.a.hanson@gmail.com
    Bobby Braswell (Rob), rbraswell@appliedgeosolutions.com
    Ian Cooke, icooke@appliedgeosolutions.com
    Tom Olson, tolson@appliedgeosolutions.com

## License (GPL v2)

    Copyright (C) 2014 Applied Geosolutions, oss@appliedgeosolutions.com

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program. If not, see <http://www.gnu.org/licenses/>
