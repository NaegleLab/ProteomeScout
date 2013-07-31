from pyramid.view import view_config
from pyramid.httpexceptions import HTTPForbidden
from ptmscout.database import annotations, user
from ptmscout.config import strings
from ptmscout.utils import webutils, decorators, crypto

@view_config(route_name='share_subset', request_method='POST', permission='private', renderer='ptmscout:templates/info/information.pt')
@decorators.get_experiment('id')
def share_subset(context, request, exp):
    subset_name = webutils.post(request, 'saved-subset-select', None)
    subset = annotations.getSubsetByName(exp.id, subset_name, request.user)

    if subset == None:
        return HTTPForbidden()

    subset.share_token = crypto.randomString(10)
    subset.save()

    share_url = request.route_url('share_subset', id=exp.id) + "?token=" + subset.share_token
    return {
        'pageTitle': strings.share_subsets_page_title % (subset_name),
        'header': strings.share_subsets_page_title % (subset_name),
        'message': strings.share_subsets_token_message % ( subset_name, exp.name, share_url, share_url )
        }

@view_config(route_name='share_subset', request_method='GET', permission='private', renderer='ptmscout:templates/info/information.pt')
@decorators.get_experiment('id')
def view_shared_subset(context, request, exp):
    share_token = webutils.get(request, 'token', None)

    subset = annotations.getSubsetByShareToken(share_token)

    if subset==None:
        return HTTPForbidden()

    owner = user.getUserById(subset.owner_id)

    nsubset = subset.copy()
    nsubset.owner_id = request.user.id
    nsubset.name = "%s: %s" % (owner.username, nsubset.name)
    nsubset.save()

    annotation_set_id = subset.annotation_set_id
    if annotation_set_id is not None and not request.user.canViewAnnotations(annotation_set_id):
        ap = annotations.AnnotationPermission()
        ap.user_id = request.user.id
        ap.annotation_set_id = annotation_set_id
        ap.access_level = 'view'

        request.user.annotation_sets.append(ap)
        request.user.saveUser()

    subset_page = request.route_url('experiment_subset', id=exp.id)
    return {
            'pageTitle': strings.share_subsets_page_title % (subset.name),
            'header': strings.share_subsets_page_title % (subset.name),
            'message': strings.share_subset_success_message % (exp.name, subset_page)
            }
