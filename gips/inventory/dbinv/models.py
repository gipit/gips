from __future__ import unicode_literals

from django.db import models


class Asset(models.Model):
    """Inventory for assets, which are downloaded files from data sources.

    Assets undergo little or no processing.  GIPS generates products
    from assets.
    """
    # max_length chosen somewhat arbitrarily since neither of our expected DB backends care (django
    # ORM is a stickler for no discernible reason).
    driver = models.TextField(db_index=True)   # eg 'modis' or 'landsat'
    asset  = models.TextField(db_index=True)   # 'MYD11A1'
    sensor = models.TextField(db_index=True)   # 'MYD'
    tile   = models.TextField(db_index=True)   # 'h12v04'
    date   = models.DateField(db_index=True)   # of observation, not production
    name   = models.TextField()                # file name including full path

    class Meta:
        # These four columns uniquely identify an asset file
        unique_together = ('driver', 'asset', 'tile', 'date')


class Product(models.Model):
    """Inventory for products, which GIPS generates from Assets.

    Main difference from Assets is to replace the asset field with a
    product field.
    """
    # max_length chosen somewhat arbitrarily since neither of our expected DB backends care (django
    # ORM is a stickler for no discernible reason).
    driver  = models.TextField(db_index=True)   # eg 'modis' or 'landsat'
    product = models.TextField(db_index=True)   # 'fsnow'
    sensor  = models.TextField(db_index=True)   # 'MYD'
    tile    = models.TextField(db_index=True)   # 'h12v04'
    date    = models.DateField(db_index=True)   # of observation, not production
    name    = models.TextField()                # file name including full path

    class Meta:
        # These four columns uniquely identify an asset file
        unique_together = ('driver', 'product', 'tile', 'date')
