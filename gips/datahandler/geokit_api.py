#!/usr/bin/env python

import datetime

from django.forms.models import model_to_dict
from gips.inventory import orm
orm.setup()
from gips.inventory.dbinv.models import DataVariable, Result
from gips.datahandler import api


def get_datacatalog ():
    """Return list of installed DataVariables as dicts"""
    dv_list = []
    dvs = DataVariable.objects.all()
    for dv in dvs:
        d = model_to_dict(dv)
        d['asset'] = eval(d['asset'])
        # xmlrpclib doesn't like datetime.date, convert to datetime.datetime
        d['start_date'] = datetime.datetime.fromordinal(d['start_date'].toordinal())
        dv_list.append(d)
    return dv_list


def submit_job(site, variable, spatial_spec, temporal_spec):
    """
    For now, this is a trivial wrapper on datahandler.api.submit_job.
    In the future, will take more geokit-friendly args and translate to
    whatever gips needs.
    """
    if 'key' not in spatial_spec:
        spatial_spec['key'] = 'shaid'
    return api.submit_job(site, variable, spatial_spec, temporal_spec)


def job_status(job_id):
    job_status, product_status = api.job_status(job_id)

    # for now, get rid of the statuses not applicable to products
    product_status.pop('remote', None)
    product_status.pop('retry', None)
    product_status.pop('post-processing', None)
    product_status.pop('pp-scheduled', None)
    product_status.pop('initializing', None)

    return job_status, product_status
    

def stats_request_results(request_id):
    """Get all the results in the specified result set.

    Input dict:
    job -- Job primary key
    """
    if not request_id.has_key('job'):
        raise KeyError('Job ID is a required identifier')
    qs = Result.objects.filter(job=request_id['job'])
    result_list = []
    for result in qs:
        r = model_to_dict(result)
        # xmlrpclib doesn't like datetime.date, convert to datetime.datetime
        r['date'] = datetime.datetime.fromordinal(r['date'].toordinal())
        result_list.append(r)

    return result_list


def stats_request_results_filter(request_id, filters):
    """Get all result objects that match specified filters.

    Inputs:
    request_id:
        job -- string identifier
    filters:
        a dictionary containing one or more django filters
    """
    if not request_id.has_key('job'):
        raise KeyError('Job ID is a required identifier')
    qs = Result.objects.filter(job=request_id['job'])
    qs.filter(**filters)
    result_list = []
    for result in qs:
        r = model_to_dict(result)
        # xmlrpclib doesn't like datetime.date, convert to datetime.datetime
        r['date'] = datetime.datetime.fromordinal(r['date'].toordinal())
        result_list.append(r)

    return result_list
    #return Result.objects.filter(feature_set).filter(filters)

