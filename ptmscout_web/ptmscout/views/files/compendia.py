from pyramid.view import view_config
from pyramid.response import FileResponse
from ptmscout.config import settings, strings
import os, time
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
                create_file_entry( request, 'proteomescout_everything.tsv', 'All proteins and modifications', listing ),
                create_file_entry( request, 'proteomescout_vertebrata.tsv', 'All vertebrate protein modifications', listing ),
                create_file_entry( request, 'proteomescout_mammalia.tsv', 'All mammalian protein modifications', listing ),
                create_file_entry( request, 'proteomescout_phosphorylation.tsv', 'All species phosphorylations', listing ),
                create_file_entry( request, 'proteomescout_acetylation.tsv', 'All species acetylation', listing ),
                create_file_entry( request, 'proteomescout_methylation.tsv', 'All species methylation', listing ),
                create_file_entry( request, 'proteomescout_ubiquitination.tsv', 'All species ubiquitination', listing ),
                create_file_entry( request, 'proteomescout_glycosylation.tsv', 'All species glycosylation', listing ),
            ]

    return {'files': files,
            'desc': strings.compendia_download_page_desc,
            'pageTitle': strings.compendia_download_page_title}

@view_config(route_name='compendia_download', permission='private')
def compendia_download(context, request):
    listing = None
    with open(os.path.join(settings.ptmscout_path, settings.export_file_path, 'listing.pyp'),'r') as listing_file:
        listing = pickle.load(listing_file)

    fname = request.matchdict['name']
    fpath = os.path.join(settings.ptmscout_path, settings.export_file_path, fname)
    if not os.path.exists(fpath):
        raise HTTPNotFound()

    t = time.strptime( listing[fname]['date'] )
    strt = time.strftime("%Y%m%d", t)
    response = FileResponse(fpath, request)
    response.content_type = 'text/tab-separated-values'
    response.content_disposition = 'attachment; filename="%s"' % ('%s_%s.tsv' % (fname[:fname.find('.tsv')], strt))

    return response
