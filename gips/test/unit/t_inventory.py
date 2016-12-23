"""Unit tests for inventory/*, mainly the Inventory classes."""

import os
import datetime

import pytest
import mock
from django.forms.models import model_to_dict

from .data import asset_filenames, expected_assets, expected_products

from gips.inventory import DataInventory, dbinv
from gips.data.modis.modis import modisData, modisAsset
from gips.core import SpatialExtent, TemporalExtent


@pytest.fixture
def di_init_mocks(mocker):
    """Appropriate mocks for DataInventory.__init__."""
    modisData.fetch = mocker.Mock() # behavior not super important
    modisData.Asset.archive = mocker.Mock() # has to return list of modisAssets
    return (modisData.fetch,
            modisData.Asset.archive,
            mock.patch('gips.inventory.dbinv.add_asset'))


@pytest.mark.django_db
def t_data_inventory_db_save(di_init_mocks):
    """Confirm DataInventory() saves to the inventory DB on fetch."""
    from gips.inventory.dbinv import models
    (m_fetch, m_archive, m_add_asset) = di_init_mocks
    # need some filenames . . .
    assets = [modisAsset(fn) for fn in asset_filenames]
    # normally _archivefile would set this attrib:
    [setattr(asset, 'archived_filename', asset.filename) for asset in assets]
    m_archive.return_value = assets

    # specify spatial & temporal
    tiles = ['h12v04', 'h13v05', 'h12v05', 'h13v04']
    se = SpatialExtent(modisData, tiles, 0.0, 0.0)
    te = TemporalExtent('2012-12-01,2012-12-03')
    # can ignore rest of **kwargs in constructor due to mocking
    di = DataInventory(modisData, se, te, fetch=True)

    # confirm DB has correct entries in it now
    rows = [model_to_dict(ao) for ao in models.Asset.objects.all()]
    [a.pop('id') for a in rows] # don't care what the keys are
    actual = {r['asset']: r for r in rows} # enforce an order so whims of hashing doesn't ruin test

    assert expected_assets == actual


@pytest.mark.django_db
def t_data_inventory_load():
    """Confirm DataInventory() loads assets and products from the inventory DB.

    In particular check the resulting hierarchy of objects for correctness."""

    # first load DB with data - start with values that should make it through to assertions
    for asset_fields in expected_assets.values():
        dbinv.add_asset(**asset_fields)
    for product_fields in expected_products.values():
        dbinv.add_product(**product_fields)

    # load more data, intentionally putting some 'misses' into the DB:
    # wrong driver
    dbinv.add_asset(**dict(expected_assets['MYD11A1'], driver='prism'))
    dbinv.add_product(**dict(expected_products['quality'], driver='merra'))
    # wrong date
    wrong_date = datetime.datetime(2012, 12, 4)
    dbinv.add_asset(**dict(expected_assets['MYD11A1'], date=wrong_date))
    dbinv.add_product(**dict(expected_products['quality'], date=wrong_date))
    # wrong tile
    dbinv.add_asset(**dict(expected_assets['MYD11A1'], tile='h12v03'))
    dbinv.add_product(**dict(expected_products['quality'], tile='h12v03'))
    # wrong product (DataInventory(products=foo))
    dbinv.add_product(**dict(expected_products['quality'], product='clouds'))

    # instantiate the class under test
    tiles = ['h12v04', 'h12v05', 'h13v04', 'h13v05']
    products = ['fsnow', 'temp8td', 'quality', 'landcover']
    se = SpatialExtent(modisData, tiles, 0.0, 0.0)
    te = TemporalExtent('2012-12-01,2012-12-03')
    di = DataInventory(modisData, se, te, products)

    # confirm DataInventory object is correct
    assert len(di.data) == 3
    for tiles_obj in di.data.values():
        assert len(tiles_obj.tiles) == 1
        data_obj = tiles_obj.tiles['h12v04']
        assert len(data_obj.assets) == 1

        # compare the asset object
        asset_obj = data_obj.assets.values()[0]
        ea = expected_assets[asset_obj.asset]
        for field in ('asset', 'date', 'sensor', 'tile'):
            assert ea[field] == getattr(asset_obj, field)
        assert ea['name'] == asset_obj.filename

        # compare the stored product
        assert len(data_obj.filenames) == 1
        ((sensor, product), fname) = (data_obj.filenames.items())[0]
        ep = expected_products[product]
        assert (ep['sensor'] == sensor and
                ep['product'] == product and
                ep['name'] == fname)
