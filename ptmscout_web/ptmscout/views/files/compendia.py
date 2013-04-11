from pyramid.view import view_config
from pyramid.response import FileResponse
from ptmscout.config import settings, strings
import os
from pyramid.httpexceptions import HTTPNotFound
import pickle

def create_file_entry(request, fn, desc, listing):
    entry = {'link': request.route_url('compendia_download', name=fn),
             'name': fn,
             'desc': desc,
             'size': listing[fn]['size'],
             'contents': '%d proteins, %d modifications' % (listing[fn]['proteins'], listing[fn]['modifications']),
             'date': listing[fn]['date']}
    return entry

@view_config(route_name='compendia', renderer='ptmscout:templates/info/filelist.pt', permission='private')
def compendia_listing(request):
    listing = None
    with open(os.path.join(settings.ptmscout_path, settings.export_file_path, 'listing.pyp'),'r') as listing_file:
        listing = pickle.load(listing_file)

    files = [
                create_file_entry( request, 'everything.tsv', 'All proteins and modifications', listing ),
                create_file_entry( request, 'vertebrata.tsv', 'All vertebrate protein modifications', listing ),
                create_file_entry( request, 'mammalia.tsv', 'All mammalian protein modifications', listing ),
                create_file_entry( request, 'phosphorylation.tsv', 'All species phosphorylations', listing ),
                create_file_entry( request, 'acetylation.tsv', 'All species acetylation', listing ),
                create_file_entry( request, 'methylation.tsv', 'All species methylation', listing ),
                create_file_entry( request, 'ubiquitination.tsv', 'All species ubiquitination', listing ),
                create_file_entry( request, 'glycosylation.tsv', 'All species glycosylation', listing ),
            ]

    return {'files': files,
            'desc': strings.compendia_download_page_desc,
            'pageTitle': strings.compendia_download_page_title}

@view_config(route_name='compendia_download', permission='private')
def compendia_download(context, request):
    fname = request.matchdict['name']
    fpath = os.path.join(settings.ptmscout_path, settings.export_file_path, fname)
    if not os.path.exists(fpath):
        raise HTTPNotFound()

    return FileResponse(fpath, request)
