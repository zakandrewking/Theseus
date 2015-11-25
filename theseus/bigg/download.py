# -*- coding: utf-8 -*-

try:
    # Python 3
    from urllib.request import urlopen
    get_charset = lambda headers: headers.get_content_charset()
except ImportError:
    # Python 2
    from urllib2 import urlopen
    get_charset = lambda headers: headers.getparam('charset')

import cobra


def _add_url_prefix(path):
    return 'http://bigg.ucsd.edu/api/v2/%s' % path.lstrip('/')

def download_model(model_id):
    """Download a COBRA model from the BiGG Models database.

    TODO warn that the unicode downloaded from bigg_models will not play nice
    with libsbml. Ask Ali about this.

    Arguments
    ---------

    model_id: The ID for a model in the BiGG Models (http://bigg.ucsd.edu)

    """
    url = _add_url_prefix('/models/%s/download' % model_id)
    response = urlopen(url)
    json_str = response.read().decode(get_charset(response.headers) or 'utf-8')
    model = cobra.io.json.from_json(json_str)
    return model
