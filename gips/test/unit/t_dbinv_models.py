import pytest
import django.db

from gips.inventory.dbinv import models

@pytest.mark.django_db
def t_asset_field_uniqueness():
    """Confirm uniqueness constraint is obeyed."""
    params = dict(
        driver='modis', asset='MOD10A1', sensor='MOD', tile='h11v02', date='2010-11-14',
        name='/repo/modis/tiles/h11v02/2010318/MOD10A1.A2010318.h11v02.005.2010320140621.hdf')
    a1 = models.Asset(**params)
    a1.save()

    # these two can vary without preventing the uniqueness constraint to be violated
    params['sensor'] = 'xxx'
    params['name'] = '/repo/modis/tiles/h11v02/2010318/xxx.hdf'
    a2 = models.Asset(**params)
    with pytest.raises(django.db.IntegrityError):
        a2.save()


@pytest.mark.django_db
def t_product_field_uniqueness():
    """Confirm uniqueness constraint is obeyed."""
    params = dict(
        driver='modis', product='snow', sensor='MCD', tile='h11v02', date='2010-11-14',
        name='/repo/modis/tiles/h11v02/2010318/h11v02_2010318_MCD_snow.tif')
    a1 = models.Product(**params)
    a1.save()

    # these two can vary without preventing the uniqueness constraint to be violated
    params['sensor'] = 'xxx'
    params['name'] = '/repo/modis/tiles/h11v02/2010318/xxx.tif'
    a2 = models.Product(**params)
    with pytest.raises(django.db.IntegrityError):
        a2.save()
