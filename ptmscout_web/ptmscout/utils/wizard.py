class HtmlLiteral(object):
    def __init__(self, html):
        self.html = html
    
    def __html__(self):
        return self.html


class WizardNavigation(object):
    def __init__(self, request):
        self.request = request
        self.page_routes = []
        self.page_urls = {}
        self.page_names = {}
        self.page_enable = {}
        self.current_page = None

    def add_page(self, page_route, page_name, enable, **kwargs):
        self.page_routes.append(page_route)
        self.page_urls[page_route] = self.request.route_url(page_route, **kwargs)
        self.page_names[page_route] = page_name
        self.page_enable[page_route] = enable

    def set_page(self, page_route):
        self.current_page = page_route

    def next_page_url(self):
        idx = self.page_routes.index(self.current_page)
        return self.page_urls[ self.page_routes[idx+1] ]

    def prev_page_url(self):
        idx = self.page_routes.index(self.current_page)
        return self.page_urls[ self.page_routes[idx-1] ]

    def has_next(self):
        idx = self.page_routes.index(self.current_page)
        return idx < len(self.page_routes) - 1

    def has_prev(self):
        idx = self.page_routes.index(self.current_page)
        return idx > 0

    def nav_crumb(self):
        bc_template   = """<ul class="wizard-breadcrumb">%s</ul>"""
        item_template = """<li><a href="%s">%s</a></li>"""
        cur_template  = """<li class="current"><span>%s</span></li>"""
        dis_template  = """<li class="future"><span>%s</span></li>"""
        nav_items = []

        for page_route in self.page_routes:
            if self.current_page == page_route:
                nav_items.append( cur_template % ( self.page_names[page_route] ) )
            elif self.page_enable[page_route]:
                nav_items.append( item_template % ( self.page_urls[page_route], self.page_names[page_route] ) )
            else:
                nav_items.append( dis_template % ( self.page_names[page_route] ) )

        return HtmlLiteral(bc_template % ("\n".join(nav_items)))
