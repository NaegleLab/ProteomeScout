from ptmscout.utils.url_data import URLBuilder

def experiment_filter(fn):
    def build_filter(*args):
        urlfilter = URLBuilder()
        urlfilter.create_field('experiment_id', URLBuilder.INTEGER)
        request = args[-1]
        urlfilter.decode_request(request)
        request.urlfilter = urlfilter

        return fn(*args)

    build_filter.__name__ = fn.__name__
    return build_filter
