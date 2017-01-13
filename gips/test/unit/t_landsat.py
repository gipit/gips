import pytest

from ...data.landsat import landsat

def actual(la):
    # TODO: la.visbands, la.lwbands, and possibly la.meta
    return (la.filename, la.asset, la.version, la.sensor,
            la.date.year, la.date.timetuple().tm_yday)

dn_fn = '/data-root/landsat/tiles/012030/2015352/LC80120302015352LGN00.tar.gz'
sr_fn = '/data-root/landsat/tiles/026035/2015162/LC80260352015162-SC20160317174932.tar.gz'

@pytest.mark.parametrize("fn, expected", [
    (dn_fn, (dn_fn, 'DN', 0, 'LC8',   2015, 352)),
    (sr_fn, (sr_fn, 'SR', 1, 'LC8SR', 2015, 162)),
])
def t_landsatAsset_constructor(fn, expected):
    la = landsat.landsatAsset(fn)
    assert expected == actual(la)
