from ptmscout.utils.url_data import URLBuilder

def experiment_filter(fn):
    def build_filter(request):
        urlfilter = URLBuilder()
        urlfilter.create_field('experiment_id', URLBuilder.INTEGER)
        urlfilter.decode_request(request)
        request.urlfilter = urlfilter

        return fn(request)

    build_filter.__name__ = fn.__name__
    return build_filter
