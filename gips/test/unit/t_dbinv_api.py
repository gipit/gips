import os, fnmatch, datetime

import pytest
import django.db
from django.forms.models import model_to_dict

from .data import asset_filenames, product_filenames, expected_assets, expected_products

from gips.inventory.dbinv import models
from gips.inventory.dbinv import (rectify_assets, rectify_products, list_tiles, add_asset,
        asset_search)
from gips.inventory import dbinv
from gips.data.modis import modisAsset, modisData


@pytest.mark.django_db
def t_rectify_assets(mocker):
    # construct plausible file listing & mock it into glob outcome
    path = modisAsset.Repository.data_path()
    rubbish_filenames = [os.path.join(path, fn) for fn in (
        # landsat is wrong obviously; shouldn't show up in results
        '012030/2015352/012030_2015352_LC8_rad-toa.tif',
        '012030/2015352/012030_2015352_LC8_bqashadow.tif',
        '012030/2015352/012030_2015352_LC8_acca.tif',
        '012030/2015352/LC80120302015352LGN00_MTL.txt',
        # erroneous assets
        'h12v04/2012338/A2003266.h12v09.005.hdf',  # missing generation stamp per issue #79
        'h12v04/2012338/bork.hdf',                 # not enough period-delimited substrings
        'h12v04/2012338/_2004001._02_06.051.2014287173613.hdf', # try to trip it up with bad chars
    )]
    all_filenames = product_filenames + rubbish_filenames + asset_filenames # should skip products
    # simulate iglob() - match against artificial filenames
    mock_iglob = mocker.patch('gips.inventory.dbinv.glob.iglob')
    mock_iglob.side_effect = lambda pat: [fn for fn in all_filenames if fnmatch.fnmatchcase(fn, pat)]

    # items that are stale (in the DB but no longer in the filesystem) should be deleted; test this
    # by preloading stale data into the DB; these should be deleted during rectify_assets() call:
    for fn in ('h12v04/2012338/MCD43A4.A2012338.h12v04.006.2016112020013.hdf',
               'h12v04/2012337/MCD43A4.A2012337.h12v04.006.2016112013509.hdf',
               'h12v04/2012336/MOD10A1.A2012336.h12v04.005.2012339213007.hdf'):
        full_fn = os.path.join(path, fn)
        ma = modisAsset(full_fn)
        asset = models.Asset(
            asset =ma.asset,
            sensor=ma.sensor,
            tile  =ma.tile,
            date  =ma.date,
            name  =full_fn,
            driver='modis',
        )
        asset.save()

    # run the function under test
    rectify_assets(modisAsset)

    # load data for inspection
    rows = [model_to_dict(ao) for ao in models.Asset.objects.all()]
    [a.pop('id') for a in rows] # don't care what the keys are
    actual = {r['asset']: r for r in rows} # enforce an order so whims of hashing doesn't ruin test

    assert len(expected_assets) == len(rows) and expected_assets == actual


@pytest.mark.django_db
def t_rectify_products(mocker):
    # construct plausible file listing & mock it into glob outcome
    path = modisAsset.Repository.data_path()
    rubbish_filenames = [os.path.join(path, fn) for fn in (
        '012030/2015352/LC80120302015352LGN00_MTL.txt',         # not a tiff
        'h12v04/2012337/h12v04_2012337_MOD_temp8td.not-a-tif',  # also not a tiff
        'h12v04/2012336/h12v04_20_12_336_MCD_quality.tif',      # too many _-delimited tokens
        'h12v04/2012338/h12v04_2012338_MCD.tif',                # not enough _-delimited tokens
    )]
    all_filenames = rubbish_filenames + asset_filenames + product_filenames
    # simulate iglob() - match against artificial filenames
    mock_iglob = mocker.patch('gips.inventory.dbinv.glob.iglob')
    mock_iglob.side_effect = lambda pat: [fn for fn in all_filenames if fnmatch.fnmatchcase(fn, pat)]

    # handle staleness of DB entries similarly to t_rectify_assets:
    for fn in ('h09v09/2012336/h09v09_2012336_MCD_quality.tif', # just product_filenames
               'h09v09/2012337/h09v09_2012337_MOD_temp8td.tif', # but with a different tile
               'h09v09/2012338/h09v09_2012338_MCD_fsnow.tif'):
        full_fn = os.path.join(path, fn)
        md = modisData(search=False)
        md.ParseAndAddFiles([full_fn])
        # sanity check, not part of the test
        assert len(md.filenames) == 1 and md.name.lower() == 'modis'
        (sensor, product) = md.filenames.keys()[0] # only one so it must be the right one
        product = models.Product(
            product=product,
            sensor =sensor,
            tile   ='h09v09',
            date   =md.date,
            name   =full_fn,
            driver ='modis',
        )
        product.save()

    # run the function under test
    rectify_products(modisData)

    # load data for inspection
    rows = [model_to_dict(po) for po in models.Product.objects.all()]
    [a.pop('id') for a in rows] # don't care what the keys are
    actual = {r['product']: r for r in rows} # enforce an order so whims of hashing doesn't ruin test

    assert len(expected_products) == len(rows) and expected_products == actual


@pytest.fixture
def basic_asset_db(db):
    # This data isn't entirely valid but is correct enough for simple tests.  Also unicode isn't super
    # relevant here; it's an artifact of cutpasting
    path_prefix = modisAsset.Repository.data_path()
    assets = [
        { # basic row in table
            'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'asset':  u'MCD43A2',
            'date':   datetime.date(2012, 12, 1),
            'driver': u'modis',
            'sensor': u'MCD',
            'tile':   u'h12v04'
        },
        { # duplicate tile string to confirm nonrepetition in outcome
            'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'asset':  u'MCD43A2',
            'date':   datetime.date(2012, 12, 2),
            'driver': u'modis',
            'sensor': u'MCD',
            'tile':   u'h12v04'
        },
        { # second tile to confirm multiple tiles will be returned
            'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'asset':  u'MCD43A2',
            'date':   datetime.date(2012, 12, 2),
            'driver': u'modis',
            'sensor': u'MCD',
            'tile':   u'h13v05'
        },
        { # different driver to confirm no cross-driver pollution of results
            'name':   path_prefix + '/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf',
            'asset':  u'MCD43A2',
            'date':   datetime.date(2012, 12, 1),
            'driver': u'landsat',
            'sensor': u'MCD',
            'tile':   u'trolololo'
        },
    ]
    [models.Asset(**f).save() for f in assets]
    return assets


@pytest.mark.django_db
def t_list_tiles(basic_asset_db):
    """Test gips.inventory.dbinv for correctness using a database fixture."""
    actual = list(list_tiles('modis'))
    expected = [u'h12v04', u'h13v05']
    assert expected == actual


@pytest.mark.django_db
def t_list_dates(basic_asset_db):
    """Test gips.inventory.dbinv for correctness using a database fixture."""
    actual = list(dbinv.list_dates('modis', 'h12v04'))
    expected = [datetime.date(2012, 12, 1), datetime.date(2012, 12, 2)]
    assert expected == actual


@pytest.mark.django_db
def t_asset_search(basic_asset_db):
    """Confirm that dbinv.asset_search works for querying assets."""
    actual = [model_to_dict(a) for a in asset_search(driver='modis', tile='h12v04')]
    [a.pop('id') for a in actual] # don't care what the keys are
    expected = [d for d in basic_asset_db if d['tile'] == 'h12v04']
    assert expected == actual


@pytest.mark.django_db
def t_product_search():
    """Test dbinv.product_search."""
    # set up a few items in the DB
    base_vals = {
        'product': u'some-product',
        'sensor':   u'some-sensor',
        'tile':     u'some-tile',
        'date':     datetime.date(2099, 12, 31),
        'name':     u'/some/file/name.xtn',
        'driver':   u'some-driver',
    }
    for n in range(3):
        vals = dict(base_vals)
        vals['tile'] = 'tile-' + str(n)
        models.Product(**vals).save()

    # perform the query & extract results for inspection
    actual = model_to_dict(dbinv.product_search(driver='some-driver', tile='tile-2')[0])

    # set up for assertions
    actual.pop('id') # don't care what the key is
    expected = dict(base_vals)
    expected['tile'] = 'tile-2'
    assert expected == actual


@pytest.mark.django_db
@pytest.mark.parametrize('mtype, call', (('asset',   dbinv.add_asset),
                                         ('product', dbinv.add_product),
                                         ('asset',   dbinv.update_or_add_asset),
                                         ('product', dbinv.update_or_add_product)))
def t_insert_new_model(mtype, call):
    """Confirm that several functions in the API save new content correctly."""
    values = {
        mtype:    u'some-' + mtype,
        'sensor': u'some-sensor',
        'tile':   u'some-tile',
        'date':   datetime.date(2099, 12, 31),
        'name':   u'/some/file/name.xtn',
        'driver': u'some-driver',
    }
    expected = dict(values)
    a = call(**values)
    returned_actual = model_to_dict(a)
    model = {'asset': models.Asset, 'product': models.Product}[mtype]
    queried_actual = model_to_dict(model.objects.get())
    expected['id'] = queried_actual['id'] # intentional small deviation from ideal test practice
    assert expected == returned_actual == queried_actual and model.objects.count() == 1


@pytest.mark.django_db
@pytest.mark.parametrize('mtype, call', (('asset',   dbinv.update_or_add_asset),
                                         ('product', dbinv.update_or_add_product)))
def t_update_existing_model(mtype, call):
    """Test create_or_add_asset when it updates an existing model."""
    primary_key = 1
    values = {
        # unique-together fields need to be the same
        mtype:    u'some-' + mtype,
        'driver': u'landsat',
        'date':   datetime.date(2012, 12, 1),
        'tile':   u'trolololo',
        # fields that can vary
        'sensor': u'some-cool-sensor',
        'name':   u'original-file-name.xtn',
    }

    # exercise code under test
    call(**values) # call once to set up
    values['name'] = 'new-file-name.xtn' # alter a value to force an update
    returned_actual = model_to_dict(call(**values)) # run the update
    model = {'asset': models.Asset, 'product': models.Product}[mtype]
    queried_actual = model_to_dict(model.objects.get())

    # perform assertions
    expected = dict(values) # now carries replaced filename
    expected['id'] = queried_actual['id'] # intentional small deviation from ideal test practice
    assert expected == returned_actual == queried_actual and model.objects.count() == 1
