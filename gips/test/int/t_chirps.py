import os
import shutil
import datetime

import pytest

from gips.inventory import dbinv
from gips.inventory import orm
from gips.data import chirps

@pytest.fixture
def scene_dir_setup():
    asset_bn = 'africa-daily-chirps-v2.0.1997.07.01.tif.gz'
    product_bn = 'africa_1997182_chirps_precip.tif'
    tile_fp = os.path.join(chirps.chirpsRepository.path('tiles'), 'africa')
    scene_fp = os.path.join(tile_fp, '1997182')

    target_asset_fp = os.path.join(scene_fp, asset_bn)
    target_product_fp = os.path.join(scene_fp, product_bn)
    link_target = '/vsigzip/' + target_asset_fp

    # technically a race condition here, but for integration tests, probably close enough
    if os.path.lexists(target_product_fp):
        raise IOError('file obstructing test, aborting: ' + target_product_fp)

    orm.setup()
    # setup by moving the asset into the archive, but only if needed
    use_fake_asset = not os.path.exists(target_asset_fp)
    if use_fake_asset:
        source_asset_fp = os.path.join(os.path.dirname(__file__), 'data', asset_bn)
        made_dirs = [d for d in tile_fp, scene_fp if not os.path.exists(d)]
        [os.mkdir(d) for d in made_dirs]
        shutil.copy(source_asset_fp, target_asset_fp)
        dbinv.rectify_assets(chirps.chirpsAsset)

    yield target_product_fp, link_target

    # cleanup (using a fake DB means don't have to rectify post-delete)
    if os.path.lexists(target_product_fp):
        os.remove(target_product_fp)
    if use_fake_asset:
        os.remove(target_asset_fp)
        [os.rmdir(d) for d in made_dirs]

@pytest.mark.django_db
def t_chirps_product_symlink(mocker, scene_dir_setup):
    target_product_fp, link_target = scene_dir_setup

    date = datetime.date(1997, 7, 1)
    cd = chirps.chirpsData('africa', date)

    cd.needed_products = mocker.Mock()
    m_np = cd.needed_products.return_value
    m_np.requested = {'precip': ['precip']}
    m_np.__len__ = lambda self: 1
    cd.AddFile = mocker.Mock() # to confirm correct call (and save rectification)

    # not a good sign when mocking is more convenient than finding out how to pass arguments
    cd.process(products='dont-care')

    assert (
        os.path.islink(target_product_fp)
        and os.readlink(target_product_fp) == link_target
        and cd.AddFile.call_count == 1
        and cd.AddFile.call_args == mocker.call('chirps', 'precip', target_product_fp)
    )
